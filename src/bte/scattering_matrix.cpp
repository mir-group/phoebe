#include "scattering_matrix.h"
#include "constants.h"
#include "mpiHelper.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <nlohmann/json.hpp>
#include <numeric> // std::iota
#include <set>
#include <utility>
#include "io.h"
#ifdef HDF5_AVAIL
#include <highfive/H5Easy.hpp>
#endif

ScatteringMatrix::ScatteringMatrix(Context &context_,
                                   StatisticsSweep &statisticsSweep_,
                                   BaseBandStructure &innerBandStructure_,
                                   BaseBandStructure &outerBandStructure_)
    : context(context_), statisticsSweep(statisticsSweep_),
      innerBandStructure(innerBandStructure_),
      outerBandStructure(outerBandStructure_)  {

  internalDiagonal = std::make_shared<VectorBTE>(statisticsSweep, outerBandStructure, 1);

  numStates = int(outerBandStructure.irrStateIterator().size());
  numCalculations = statisticsSweep.getNumCalculations();

  dimensionality_ = int(context.getDimensionality());
  highMemory = context.getScatteringMatrixInMemory();

  // we want to know the state index of acoustic modes at gamma,
  // so that we can set their populations to zero
  // TODO shouldn't we check both bandstructures?
  // Or move this to specific subclasses
  if (outerBandStructure.getParticle().isPhonon()) {
    for (int iBte = 0; iBte < numStates; iBte++) {
      auto iBteIdx = BteIndex(iBte);
      StateIndex isIdx = outerBandStructure.bteToState(iBteIdx);
      double en = outerBandStructure.getEnergy(isIdx);
      // TODO need to standardize this
      if (en < 0.1 / ryToCmm1) { // cutoff at 0.1 cm^-1
        excludeIndices.push_back(iBte);
      }

      Eigen::Vector3d k = outerBandStructure.getWavevector(isIdx);
      if (k.squaredNorm() > 1e-8 && en < 0.) {
        Warning("Found a phonon mode q!=0 with negative energies. "
                "Consider improving the quality of your DFT phonon calculation.\n");
      }
    }
  }

  if (context.getConstantRelaxationTime() > 0.) {
    constantRTA = true;
    return;
  }

  // we only output U and N in RTA, as otherwise we would need to do either:
  // 1) the full scattering matrix calc multiple times
  // 2) store 3 copies of the scattering matrix,
  // both of which are bad options
  if ( context.getOutputUNTimes() ) {

    outputUNTimes = true;
    internalDiagonalNormal = std::make_shared<VectorBTE>(statisticsSweep, outerBandStructure, 1);
    internalDiagonalUmklapp = std::make_shared<VectorBTE>(statisticsSweep, outerBandStructure, 1);

    // warn users that if they run an exact bte solver, the U and N will
    // still only come out in RTA
    if( context.getSolverBTE().size() > 0 ) {
      Warning("You've set outputUNTimes to true in input file, but requested an exact BTE solver.\n"
        "Be aware that U and N separation is only currently implemented in the RTA case.");
    }
  }
}

void ScatteringMatrix::setup() {
  // note: here we want to build the matrix or its diagonal
  // builder is a pure virtual function, which is implemented in subclasses
  // c++ discourages calls to pure virtual functions in the constructor

  if (constantRTA) return; // nothing to construct

  std::vector<VectorBTE> emptyVector;

  // user info about memory
  memoryUsage();

  if (highMemory) {
    if (numCalculations > 1) {
      // note: one could write code around this
      // but the methods are very memory intensive for production runs
      Error("High memory BTE methods (scatteringMatrixInMemory=true) can only work with one "
            "temperature and/or chemical potential in a single run");
    }
    int matSize;
    if (context.getUseSymmetries()) {
      matSize = int(3 * numStates);
    } else {
      matSize = int(numStates);
    }

    // user info about memory
    if (mpi->mpiHead()) {
      double x = pow(matSize, 2) / pow(1024., 3) * sizeof(double);
      std::cout << "Allocating " << x << " GB (" << x / (mpi->getSize())
                << " GB per process) for the scattering matrix.\n"
                << std::endl;
    }

    try {
      // The block size of this matrix can really change the performance
      // of the diagonalization method, and also the parallelization
      // of the scattering rate calculation in the full matrix case!
      // Choosing numBlocks = matSize --> perfectly cyclic matrix
      // Choosing numBlocks = numMPIProcs (the default) --> perfectly block
      //
      // Jenny's note: I've benchmarked several values of block size, and
      // blocksize ~between 10-100 is the best. Hopefully this is universally true,
      // benchmarks were done in the phonon case with a matrix with matSize=~31k, 81 MPI procs
      // 64 ended up giving me slightly better performance than 16 or 176.
      // There may be logic to how to choose these sizes,
      // but for now I cannot find advice on this online.

      int nBlocks = int(matSize/64.);
      // TODO Should we move this to PMatrix constructor instead?
      // seems tiny number of blocks is a problem, but it should
      // only come up for very tiny cases where speed is not an issue,
      // therefore, we default to 4 blocks if this happens.
      // Seems like this is maybe a bug in scalapack?
      if(nBlocks < 4) nBlocks = 4;

      theMatrix = ParallelMatrix<double>(matSize, matSize, 0, 0, nBlocks, nBlocks);

    } catch(std::bad_alloc&) {
      Error("Failed to allocate memory for the scattering matrix.\n"
        "You are likely running out of memory.");
    }

    // linewidth and matrix construction
    builder(internalDiagonal, emptyVector, emptyVector);

  } else {
    // calc linewidths only
    builder(internalDiagonal, emptyVector, emptyVector);
  }
}

// NOTE: this is returning the raw diagonal.
// This includes whatever factors there are if it's Omega or not!
VectorBTE ScatteringMatrix::diagonal() {
  if (constantRTA) {
    VectorBTE diagonal(statisticsSweep, outerBandStructure, 1);
    diagonal.setConst(twoPi / context.getConstantRelaxationTime());
    return diagonal;
  } else {
    return *internalDiagonal;
  }
}

VectorBTE ScatteringMatrix::offDiagonalDot(VectorBTE &inPopulation) {

  if (highMemory) {

    VectorBTE outPopulation(statisticsSweep, outerBandStructure, 3);
    // note: we are assuming that ScatteringMatrix has numCalculations = 1

    auto allLocalStates = theMatrix.getAllLocalStates();
    size_t numAllLocalStates = allLocalStates.size();

#pragma omp parallel
    {
      Eigen::MatrixXd dataPrivate = Eigen::MatrixXd::Zero(3, numStates);

      if (context.getUseSymmetries()) {
#pragma omp for
        for (size_t iTup = 0; iTup < numAllLocalStates; iTup++) {
          auto tup = allLocalStates[iTup];
          int iMat1 = std::get<0>(tup);
          int iMat2 = std::get<1>(tup);
          auto t1 = getSMatrixIndex(iMat1);
          auto t2 = getSMatrixIndex(iMat2);
          int iBte1 = std::get<0>(t1).get();
          int iBte2 = std::get<0>(t2).get();
          if (iBte1 == iBte2)
            continue;
          int i = std::get<1>(t1).get();
          int j = std::get<1>(t2).get();
          dataPrivate(i, iBte1) +=
              theMatrix(iMat1, iMat2) * inPopulation(0, j, iBte2);
        }
      } else {
#pragma omp for
        for (size_t iTup = 0; iTup < numAllLocalStates; iTup++) {
          auto tup = allLocalStates[iTup];
          auto iBte1 = std::get<0>(tup);
          auto iBte2 = std::get<1>(tup);
          if (iBte1 == iBte2)
            continue;
          for (int i : {0, 1, 2}) {
            dataPrivate(i, iBte1) +=
                theMatrix(iBte1, iBte2) * inPopulation(0, i, iBte2);
          }
        }
      }

#pragma omp critical
      for (int iBte1=0; iBte1<numStates; ++iBte1) {
        for (int i : {0, 1, 2}) {
          outPopulation(0, i, iBte1) += dataPrivate(i, iBte1);
        }
      }

    } // end pragma omp parallel

    mpi->allReduceSum(&outPopulation.data);
    return outPopulation;
  } else {
    VectorBTE outPopulation = dot(inPopulation);
#pragma omp parallel for collapse(3) default(none)                             \
    shared(outPopulation, internalDiagonal, inPopulation, numCalculations, numStates)
    for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
      for (int iDim = 0; iDim < 3; iDim++) {
        for (int iBte = 0; iBte < numStates; iBte++) {
          outPopulation(iCalc, iDim, iBte) -= internalDiagonal->operator()(iCalc, 0, iBte) *
                                              inPopulation(iCalc, iDim, iBte);
        }
      }
    }
    return outPopulation;
  }
}

std::vector<VectorBTE>
ScatteringMatrix::offDiagonalDot(std::vector<VectorBTE> &inPopulations) {
  // outPopulation = outPopulation - internalDiagonal * inPopulation;
  std::vector<VectorBTE> outPopulations = dot(inPopulations);
  for (unsigned int iVec = 0; iVec < inPopulations.size(); iVec++) {
#pragma omp parallel for collapse(3)
    for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
      for (int iDim = 0; iDim < 3; iDim++) {
        for (int iBte = 0; iBte < numStates; iBte++) {
          outPopulations[iVec](iCalc, iDim, iBte) -=
              internalDiagonal->operator()(iCalc, 0, iBte) *
              inPopulations[iVec](iCalc, iDim, iBte);
        }
      }
    }
  }
  return outPopulations;
}

VectorBTE ScatteringMatrix::dot(VectorBTE &inPopulation) {
  if (highMemory) {
    VectorBTE outPopulation(statisticsSweep, outerBandStructure, 3);
    // note: we are assuming that ScatteringMatrix has numCalculations = 1

    auto allLocalStates = theMatrix.getAllLocalStates();
    size_t numAllLocalStates = allLocalStates.size();

#pragma omp parallel
    {
      Eigen::MatrixXd dataPrivate = Eigen::MatrixXd::Zero(3, numStates);

      if (context.getUseSymmetries()) {
#pragma omp for
        for (size_t iTup = 0; iTup < numAllLocalStates; iTup++) {
          auto tup = allLocalStates[iTup];
          int iMat1 = std::get<0>(tup);
          int iMat2 = std::get<1>(tup);
          auto t1 = getSMatrixIndex(iMat1);
          auto t2 = getSMatrixIndex(iMat2);
          int iBte1 = std::get<0>(t1).get();
          int iBte2 = std::get<0>(t2).get();

          if (std::find(excludeIndices.begin(), excludeIndices.end(), iBte1) != excludeIndices.end())
            continue;
          if (std::find(excludeIndices.begin(), excludeIndices.end(), iBte2) != excludeIndices.end())
            continue;

          int i = std::get<1>(t1).get();
          int j = std::get<1>(t2).get();
          dataPrivate(i, iBte1) +=
              theMatrix(iMat1, iMat2) * inPopulation(0, j, iBte2);
        }
      } else {
#pragma omp for
        for (size_t iTup = 0; iTup < numAllLocalStates; iTup++) {
          auto tup = allLocalStates[iTup];
          auto iBte1 = std::get<0>(tup);
          auto iBte2 = std::get<1>(tup);

          if (std::find(excludeIndices.begin(), excludeIndices.end(), iBte1) != excludeIndices.end())
            continue;
          if (std::find(excludeIndices.begin(), excludeIndices.end(), iBte2) != excludeIndices.end())
            continue;

          for (int i : {0, 1, 2}) {
            dataPrivate(i, iBte1) +=
                theMatrix(iBte1, iBte2) * inPopulation(0, i, iBte2);
          }
        }
      }
#pragma omp critical
      for (int iBte1=0; iBte1<numStates; ++iBte1) {
        for (int i : {0, 1, 2}) {
          outPopulation(0, i, iBte1) += dataPrivate(i, iBte1);
        }
      }
    }

    mpi->allReduceSum(&outPopulation.data);
    return outPopulation;
  } else {
    VectorBTE outPopulation(statisticsSweep, outerBandStructure, 3);
    outPopulation.data.setZero();
    std::vector<VectorBTE> outPopulations;
    std::vector<VectorBTE> inPopulations;
    inPopulations.push_back(inPopulation);
    outPopulations.push_back(outPopulation);
    builder(nullptr, inPopulations, outPopulations);
    return outPopulations[0];
  }
}

std::vector<VectorBTE>
ScatteringMatrix::dot(std::vector<VectorBTE> &inPopulations) {
  if (highMemory) {
    std::vector<VectorBTE> outPopulations;
    for (auto inPopulation : inPopulations) {
      VectorBTE outPopulation(statisticsSweep, outerBandStructure, 3);
      outPopulation = dot(inPopulation);
      outPopulations.push_back(outPopulation);
    }
    return outPopulations;
  } else {
    std::vector<VectorBTE> outPopulations;
    for (unsigned int i = 0; i < inPopulations.size(); i++) {
      VectorBTE outPopulation(statisticsSweep, outerBandStructure, 3);
      outPopulations.push_back(outPopulation);
    }
    builder(nullptr, inPopulations, outPopulations);
    return outPopulations;
  }
}

// set,unset the scaling of omega = A/sqrt(bose1*bose1+1)/sqrt(bose2*bose2+1)
void ScatteringMatrix::a2Omega() {

  if (!highMemory) {
    Error("a2Omega only works if the matrix is stored in memory");
  }

  if (theMatrix.rows() == 0) {
    Error("The scattering matrix hasn't been built yet");
  }

  if (isMatrixOmega) { // it's already with the scaling of omega
    return;
  }

  int iCalc = 0; // as there can only be one temperature

  auto particle = outerBandStructure.getParticle();
  auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
  double temp = calcStatistics.temperature;
  double chemPot = calcStatistics.chemicalPotential;

  auto allLocalStates = theMatrix.getAllLocalStates();
  size_t numAllLocalStates = allLocalStates.size();
#pragma omp parallel for
  for (size_t iTup=0; iTup<numAllLocalStates; iTup++) {
    auto tup = allLocalStates[iTup];

    int iBte1, iBte2, iMat1, iMat2;
    StateIndex is1Idx(-1), is2Idx(-1);
    if (context.getUseSymmetries()) {
      iMat1 = std::get<0>(tup);
      iMat2 = std::get<1>(tup);
      auto tup1 = getSMatrixIndex(iMat1);
      auto tup2 = getSMatrixIndex(iMat2);
      BteIndex iBte1Idx = std::get<0>(tup1);
      BteIndex iBte2Idx = std::get<0>(tup2);
      iBte1 = iBte1Idx.get();
      iBte2 = iBte2Idx.get();
      is1Idx = outerBandStructure.bteToState(iBte1Idx);
      is2Idx = outerBandStructure.bteToState(iBte2Idx);
    } else {
      iBte1 = std::get<0>(tup);
      iBte2 = std::get<1>(tup);
      iMat1 = iBte1;
      iMat2 = iBte2;
      is1Idx = StateIndex(iBte1);
      is2Idx = StateIndex(iBte2);
    }
    if (std::find(excludeIndices.begin(), excludeIndices.end(), iBte1) !=
        excludeIndices.end())
      continue;
    if (std::find(excludeIndices.begin(), excludeIndices.end(), iBte2) !=
        excludeIndices.end())
      continue;

    double en1 = outerBandStructure.getEnergy(is1Idx);
    double en2 = outerBandStructure.getEnergy(is2Idx);

    // n(n+1) for bosons, n(1-n) for fermions
    double term1 = particle.getPopPopPm1(en1, temp, chemPot);
    double term2 = particle.getPopPopPm1(en2, temp, chemPot);

    if (iBte1 == iBte2) {
      internalDiagonal->operator()(0, 0, iBte1) /= term1;
    }

    theMatrix(iMat1, iMat2) /= sqrt(term1 * term2);
  }
  isMatrixOmega = true;
}

// add a flag to remember if we have A or Omega
// void ScatteringMatrix::omega2A() {
//  if (!highMemory) {
//    Error("a2Omega only works if the matrix is stored in memory");
//  }
//
//  if (theMatrix.rows() == 0) {
//    Error("The scattering matrix hasn't been built yet");
//  }
//
//  if (!isMatrixOmega) { // it's already with the scaling of A
//    return;
//  }
//
//  int iCalc = 0; // as there can only be one temperature
//
//  auto particle = outerBandStructure.getParticle();
//  auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
//  double temp = calcStatistics.temperature;
//  double chemPot = calcStatistics.chemicalPotential;
//
//  for (auto tup : theMatrix.getAllLocalStates()) {
//    auto ind1 = std::get<0>(tup);
//    auto ind2 = std::get<1>(tup);
//    if (std::find(excludeIndices.begin(), excludeIndices.end(), ind1) !=
//        excludeIndices.end())
//      continue;
//    if (std::find(excludeIndices.begin(), excludeIndices.end(), ind2) !=
//        excludeIndices.end())
//      continue;
//
//    double en1 = outerBandStructure.getEnergy(ind1);
//    double en2 = outerBandStructure.getEnergy(ind2);
//
//    // n(n+1) for bosons, n(1-n) for fermions
//    double term1 = particle.getPopPopPm1(en1, temp, chemPot);
//    double term2 = particle.getPopPopPm1(en2, temp, chemPot);
//
//    if (ind1 == ind2) {
//      internalDiagonal(iCalc, 0, ind1) *= term1;
//    }
//
//    theMatrix(ind1, ind2) *= sqrt(term1 * term2);
//  }
//  isMatrixOmega = false;
//}

// to compute the RTA, get the single mode relaxation times
VectorBTE ScatteringMatrix::getTimesFromVectorBTE(VectorBTE &diagonal) {

  // just using the shape of internalDiagonal, will be overwritten
  VectorBTE times(statisticsSweep, outerBandStructure, 1);

  if (constantRTA) {
    times.setConst(context.getConstantRelaxationTime() / twoPi);
  } else {
    if (isMatrixOmega) {
      times = diagonal.reciprocal();
    } else { // A_nu,nu = N(1+-N) / tau  -- for phonon case
      auto particle = outerBandStructure.getParticle();
      #pragma omp parallel for
      for (int iBte = 0; iBte < numStates; iBte++) {
        BteIndex iBteIdx(iBte);
        StateIndex isIdx = outerBandStructure.bteToState(iBteIdx);
        double en = outerBandStructure.getEnergy(isIdx);
        for (int iCalc = 0; iCalc < diagonal.numCalculations; iCalc++) {
          auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
          double temp = calcStatistics.temperature;
          double chemPot = calcStatistics.chemicalPotential;
          // n(n+1) for bosons, n(1-n) for fermions
          double popTerm = particle.getPopPopPm1(en, temp, chemPot);
          times(iCalc, 0, iBte) = popTerm / diagonal(iCalc, 0, iBte);
        }
      }
    }
  }
  return times;
}

// function to decide which times to prepare
// 0 -- regular linewidths/internal diagonal
// 1 -- normal processes
// 2 -- umklapp processes
VectorBTE ScatteringMatrix::getSingleModeTimes(int whichTimes) {
  if(whichTimes == 1) { //normal internal diagonal
    return getTimesFromVectorBTE(*internalDiagonalNormal);
  } else if(whichTimes == 2) { // umklapp internal diagonal
    return getTimesFromVectorBTE(*internalDiagonalUmklapp);
  } else { // regular internal diagonal of scattering matrix
    return getTimesFromVectorBTE(*internalDiagonal);
  }
}

// return the diagonal of the scattering matrix,
// converting first to conventional lifetimes if the matrix is
// symmetrized
VectorBTE ScatteringMatrix::getLinewidths() {

  if (constantRTA) {
    VectorBTE linewidths(statisticsSweep, outerBandStructure, 1);
    linewidths.setConst(twoPi / context.getConstantRelaxationTime());
    linewidths.excludeIndices = excludeIndices;
    return linewidths;
  } else {
    VectorBTE linewidths = *internalDiagonal;
    linewidths.excludeIndices = excludeIndices;
    auto particle = outerBandStructure.getParticle();

    if (isMatrixOmega) {
      return linewidths;

    } else {
      // A_nu,nu = Gamma / N(1+N) for phonons, Omega_nu,nu = Gamma for electrons
      if (particle.isElectron()) {
        Error("Attempting to use a numerically unstable quantity");
        // popTerm could be = 0
      }
      #pragma omp parallel for
      for (int iBte = 0; iBte < numStates; iBte++) {
        auto iBteIdx = BteIndex(iBte);
        StateIndex isIdx = outerBandStructure.bteToState(iBteIdx);
        double en = outerBandStructure.getEnergy(isIdx);
        for (int iCalc = 0; iCalc < internalDiagonal->numCalculations; iCalc++) {
          auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
          double temp = calcStatistics.temperature;
          double chemPot = calcStatistics.chemicalPotential;
          // n(n+1) for bosons, n(1-n) for fermions
          double popTerm = particle.getPopPopPm1(en, temp, chemPot);
          linewidths(iCalc, 0, iBte) =
              internalDiagonal->operator()(iCalc, 0, iBte) / popTerm;
        }
      }
      linewidths.excludeIndices = excludeIndices;
      return linewidths;
    }
  }
}

void ScatteringMatrix::setLinewidths(VectorBTE &linewidths) {

  // note: this function assumes that we have not
  // messed with excludeIndices. Should be fine,
  // because when we add phonon contributions to phonon scattering matrix,
  // they will exclude the same indices so long as numStates is the same

  if(internalDiagonal->numCalculations != linewidths.numCalculations) {
    Error("Attempted setting scattering matrix diagonal with"
        " an incorrect number of calculations.");
  }
  if(internalDiagonal->numStates != linewidths.numStates) {
    Error("Attempted setting scattering matrix diagonal with"
        " an incorrect number of states.");
  }

  if (isMatrixOmega) {
    // A_nu,nu = Gamma / N(1+N) for phonons, A_nu,nu = Gamma for electrons
    // because the get function removed these population factors,
    // the set function should add them back in
    internalDiagonal =  std::make_shared<VectorBTE>(linewidths);
   } else {  // phonons
    auto particle = outerBandStructure.getParticle();
    #pragma omp parallel for
    for (int iBte = 0; iBte < numStates; iBte++) {
      auto iBteIdx = BteIndex(iBte);
      StateIndex isIdx = outerBandStructure.bteToState(iBteIdx);
      double en = outerBandStructure.getEnergy(isIdx);
      for (int iCalc = 0; iCalc < linewidths.numCalculations; iCalc++) {
        auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
        double temp = calcStatistics.temperature;
        double chemPot = calcStatistics.chemicalPotential;
        // n(n+1) for bosons, n(1-n) for fermions
        double popTerm = particle.getPopPopPm1(en, temp, chemPot);
        internalDiagonal->operator()(iCalc, 0, iBte) = linewidths(iCalc, 0, iBte) * popTerm;
      }
    }
  }

  // we also need to add these to the full scattering matrix
  // we do this after setting up internal diagonal so that population
  // factors are already in place as needed

  if (highMemory) { // matrix in memory
    int iCalc = 0;
    if (context.getUseSymmetries()) {
      // numStates is defined in scattering.cpp as # of irrStates
      // from the outer band structure
      for (int iBte = 0; iBte < numStates; iBte++) {
        BteIndex iBteIdx(iBte);
        for (int i : {0, 1, 2}) {
          CartIndex iCart(i);
          int iMati = getSMatrixIndex(iBteIdx, iCart);
          for (int j : {0, 1, 2}) {
            CartIndex jCart(j);
            int iMatj = getSMatrixIndex(iBteIdx, jCart);
            theMatrix(iMati, iMatj) = 0.;
          }
          theMatrix(iMati, iMati) += internalDiagonal->operator()(iCalc, 0, iBte);
        }
      }
    } else {
      for (int is = 0; is < numStates; is++) {
        theMatrix(is, is) = internalDiagonal->operator()(iCalc, 0, is);
      }
    }
  }
}

void ScatteringMatrix::outputToHDF5(const std::string &outFileName) {

  // output the matrix elements owned by this process to file.
  // The only tricky part here is that we need to sort them
  // Appropriately in the file without all reducing. This is nearly
  // impossible to do efficiently because of how HighFive is written,
  // wherein it expects chunks of data -- which we cannot expect to have
  // in a block-cyclic distribution of the matrix. Therefore, we do something
  // terrible and we write single elements at a time.

  #ifndef HDF5_AVAIL
    Error("Need Phoebe compiled with parallel HDF5 to write SMatrix to file.");
  #elif HDF5_SERIAL
    Error("Need Phoebe compiled with parallel HDF5 to write SMatrix to file.");
  #else

   try {

      // Create the file to write to, first removing if it exists already
      std::remove(&outFileName[0]);
      HighFive::FileAccessProps fapl;
      fapl.add(HighFive::MPIOFileAccess<MPI_Comm, MPI_Info>(MPI_COMM_WORLD, MPI_INFO_NULL));
      HighFive::File file(outFileName, HighFive::File::Overwrite, fapl);

      // setup the dataset
      unsigned int globalSize = theMatrix.size();

      // Create the data-space to write gWannier to
      std::vector<size_t> dims(2);
      dims[0] = 1;
      dims[1] = size_t(globalSize);
      HighFive::DataSet dsMatrix= file.createDataSet<double>(
                "/scatteringMatrix", HighFive::DataSpace(dims));

      LoopPrint loopPrint("writing the scattering matrix to HDF5",
             "matrix elements / MPI process", theMatrix.getAllLocalStates().size());

      // now we loop over the matrix elements and write them to file
      for(auto matEl : theMatrix.getAllLocalStates()) {

        loopPrint.update();

        // get matrix row and col indices
        size_t iMat1 = std::get<0>(matEl);
        size_t iMat2 = std::get<1>(matEl);
        // convert to 1d index, used on the underlying dataset
        //size_t iGlobal = theMatrix.global2Local(iMat1, iMat2);
        size_t iGlobal = iMat2*theMatrix.cols() + iMat1;
        // get the actual value of the matrix element
        double value = theMatrix(iMat1,iMat2);

        // write the single element to the dataset
        // here, the offset is the global index, and the count is of course 1
        dsMatrix.select({0, iGlobal}, {1, 1}).write(value);
      }

      loopPrint.close();

      // ensure everything has been written
      file.flush();

    } catch (std::exception &error) {
      // catch and print any HDF5 error
      std::cerr << error.what() << std::endl;
      Error("An HDF5 error occurred while trying to write scattering matrix to HDF5.");
    }
  #endif
}

void ScatteringMatrix::outputToJSON(const std::string &outFileName) {

  if (!mpi->mpiHead())
    return;

  VectorBTE times = getSingleModeTimes();
  // this will remove the factor of pop(pop +/- 1) if it's symmetrized
  VectorBTE tmpLinewidths = getLinewidths();
  std::shared_ptr<VectorBTE> timesN;
  std::shared_ptr<VectorBTE> timesU;
  if(outputUNTimes) {
    timesN = std::make_shared<VectorBTE>(getSingleModeTimes(1));
    timesU = std::make_shared<VectorBTE>(getSingleModeTimes(2));
  }

  std::string particleType;
  auto particle = outerBandStructure.getParticle();
  double energyConversion = energyRyToEv;
  std::string energyUnit = "eV";
  std::string relaxationTimeUnit = "fs";
  // we need an extra factor of two pi, likely because of unit conversion
  // (perhaps h vs hbar)
  double energyToTime = energyRyToFs/twoPi;
  if (particle.isPhonon()) {
    particleType = "phonon";
    energyUnit = "meV";
    energyConversion *= 1000;
    relaxationTimeUnit = "ps"; // phonon times more commonly in ps
    energyToTime *= 1e-3;
    // this is a bit of a hack to deal with phel scattering, where stat sweep
    // has nonzero mu values in spite of it being a phonon case

  } else {
    particleType = "electron";
  }

  // need to store as a vector format with dimensions
  // iCalc, ik. ib, iDim (where iState is unfolded into
  // ik, ib) for the velocities and lifetimes, no dim for energies
  std::vector<std::vector<std::vector<double>>> outTimes;
  std::vector<std::vector<std::vector<double>>> outTimesU;
  std::vector<std::vector<std::vector<double>>> outTimesN;
  std::vector<std::vector<std::vector<double>>> outLinewidths;
  std::vector<std::vector<std::vector<std::vector<double>>>> velocities;
  std::vector<std::vector<std::vector<double>>> energies;
  std::vector<double> temps;
  std::vector<double> chemPots;
  std::vector<double> dopings;

  for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
    auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStatistics.temperature;
    double chemPot = calcStatistics.chemicalPotential;
    double doping = calcStatistics.doping;
    temps.push_back(temp * temperatureAuToSi);
    // this is a bit of a hack to deal with phel scattering, where stat sweep
    // has nonzero mu values in spite of it being a phonon case
    if(particle.isElectron()) {
      chemPots.push_back(chemPot * energyConversion);
    } else {
      chemPots.push_back(0);
    }
    dopings.push_back(doping);

    std::vector<std::vector<double>> wavevectorsT;
    std::vector<std::vector<double>> wavevectorsL;
    std::vector<std::vector<std::vector<double>>> wavevectorsV;
    std::vector<std::vector<double>> wavevectorsE;
    // loop over wavevectors
    for (int ik : outerBandStructure.irrPointsIterator()) {
      auto ikIndex = WavevectorIndex(ik);

      std::vector<double> bandsT;
      std::vector<double> bandsL;
      std::vector<std::vector<double>> bandsV;
      std::vector<double> bandsE;
      // loop over bands here
      // get numBands at this point, in case it's an active band structure
      for (int ib = 0; ib < outerBandStructure.getNumBands(ikIndex); ib++) {
        auto ibIndex = BandIndex(ib);
        int is = outerBandStructure.getIndex(ikIndex, ibIndex);
        StateIndex isIdx(is);
        double ene = outerBandStructure.getEnergy(isIdx);
        auto vel = outerBandStructure.getGroupVelocity(isIdx);
        bandsE.push_back(ene * energyConversion);
        int iBte = int(outerBandStructure.stateToBte(isIdx).get());
        double tau = times(iCalc, 0, iBte);
        bandsT.push_back(tau * energyToTime);
        double linewidth = tmpLinewidths(iCalc, 0, iBte);
        bandsL.push_back(linewidth * energyConversion);

        std::vector<double> iDimsV;
        // loop over dimensions
        for (int iDim : {0, 1, 2}) {
          iDimsV.push_back(vel[iDim] * velocityRyToSi);
        }
        bandsV.push_back(iDimsV);
      }
      wavevectorsT.push_back(bandsT);
      wavevectorsL.push_back(bandsL);
      wavevectorsV.push_back(bandsV);
      wavevectorsE.push_back(bandsE);
    }
    outTimes.push_back(wavevectorsT);
    outLinewidths.push_back(wavevectorsL);
    velocities.push_back(wavevectorsV);
    energies.push_back(wavevectorsE);
  }
  if(outputUNTimes) {
    for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
      std::vector<std::vector<double>> wavevectorsU;
      std::vector<std::vector<double>> wavevectorsN;
      // loop over wavevectors
      for (int ik : outerBandStructure.irrPointsIterator()) {
        auto ikIndex = WavevectorIndex(ik);
        std::vector<double> bandsN;
        std::vector<double> bandsU;
        // loop over bands here
        for (int ib = 0; ib < outerBandStructure.getNumBands(ikIndex); ib++) {
          auto ibIndex = BandIndex(ib);
          int is = outerBandStructure.getIndex(ikIndex, ibIndex);
          StateIndex isIdx(is);
          int iBte = int(outerBandStructure.stateToBte(isIdx).get());
          double tauN = timesN->operator()(iCalc, 0, iBte);
          double tauU = timesU->operator()(iCalc, 0, iBte);
          bandsN.push_back(tauN * energyToTime);
          bandsU.push_back(tauU * energyToTime);
        }
        wavevectorsN.push_back(bandsN);
        wavevectorsU.push_back(bandsU);
      }
      outTimesN.push_back(wavevectorsN);
      outTimesU.push_back(wavevectorsU);
    }
  }

  auto points = outerBandStructure.getPoints();
  std::vector<std::vector<double>> meshCoordinates;
  for (int ik : outerBandStructure.irrPointsIterator()) {
    // save the wavevectors
    auto ikIndex = WavevectorIndex(ik);
    auto coord =
        points.cartesianToCrystal(outerBandStructure.getWavevector(ikIndex));
    meshCoordinates.push_back({coord[0], coord[1], coord[2]});
  }

  // output to json
  nlohmann::json output;
  output["temperatures"] = temps;
  output["temperatureUnit"] = "K";
  output["chemicalPotentials"] = chemPots;
  output["chemicalPotentialUnit"] = "eV";
  if (particle.isElectron()) {
    output["dopingConcentrations"] = dopings;
    output["dopingConcentrationUnit"] =
        "cm$^{-" + std::to_string(context.getDimensionality()) + "}$";
  }
  output["linewidths"] = outLinewidths;
  output["linewidthsUnit"] = energyUnit;
  output["relaxationTimes"] = outTimes;
  if(outputUNTimes) {
    output["normalRelaxationTimes"] = outTimesN;
    output["umklappRelaxationTimes"] = outTimesU;
  }
  output["relaxationTimeUnit"] = relaxationTimeUnit;
  output["velocities"] = velocities;
  output["velocityUnit"] = "m/s";
  output["energies"] = energies;
  output["energyUnit"] = energyUnit;
  output["wavevectorCoordinates"] = meshCoordinates;
  output["coordsType"] = "lattice";
  output["particleType"] = particleType;
  std::ofstream o(outFileName);
  o << std::setw(3) << output << std::endl;
  o.close();
}

// TODO this feels redundant with above function, maybe could be simplified --
// plus, it's missing a bunch of things including energies!
void ScatteringMatrix::relaxonsToJSON(const std::string &outFileName,
                                      const Eigen::VectorXd &eigenvalues) {
  if (!mpi->mpiHead()) {
    return;
  }

  Eigen::VectorXd times = 1. / eigenvalues.array();

  std::string particleType;
  auto particle = outerBandStructure.getParticle();
  if (particle.isPhonon()) {
    particleType = "phonon";
  } else {
    particleType = "electron";
  }

  double energyToTime = energyRyToFs;
  double energyConversion = energyRyToEv;

  std::string energyUnit = "eV";

  // need to store as a vector format with dimensions
  // iCalc, ik. ib, iDim (where iState is unfolded into
  // ik, ib) for the velocities and lifetimes, no dim for energies
  std::vector<std::vector<double>> outTimes;
  std::vector<double> temps;
  std::vector<double> chemPots;

  for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
    auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStatistics.temperature;
    double chemPot = calcStatistics.chemicalPotential;
    temps.push_back(temp * temperatureAuToSi);
    chemPots.push_back(chemPot * energyConversion);

    std::vector<double> thisTime;
    for (auto x : times) {
      thisTime.push_back(x * energyToTime);
    }
    outTimes.push_back(thisTime);
  }

  // output to json
  nlohmann::json output;
  output["temperatures"] = temps;
  output["temperatureUnit"] = "K";
  output["chemicalPotentials"] = chemPots;
  output["relaxationTimes"] = outTimes;
  output["relaxationTimeUnit"] = "fs";
  output["particleType"] = particleType;
  std::ofstream o(outFileName);
  o << std::setw(3) << output << std::endl;
  o.close();
}

std::tuple<Eigen::VectorXd, ParallelMatrix<double>>
ScatteringMatrix::diagonalize(int numEigenvalues) {

  // user info about memory
  {
    memoryUsage();
    double xx;
    double x = 2 * pow(theMatrix.rows(), 2) / pow(1024., 3) * sizeof(xx);
    // 2 because one is for eigenvectors, another is the copy of the matrix
    std::cout << std::setprecision(4);
    if (mpi->mpiHead()) {
      std::cout << "About to allocate " << x
                << " (GB) for the scattering matrix diagonalization.\n"
                << std::endl;
    }
  }
  // diagonalize
  std::tuple<std::vector<double>, ParallelMatrix<double>> tup;
  if(numEigenvalues > numStates) { // not possible
    Error("You have requested to calculate more relaxons eigenvalues"
        " than your number of particle states.");
  }
  if(numEigenvalues > 0 && numEigenvalues != numStates) { // zero is default, calculates all of them
    // calculate some of the eigenvalues
    tup = theMatrix.diagonalize(numEigenvalues);
  } else { // calculate all
    tup = theMatrix.diagonalize();
  }
  auto eigenvalues = std::get<0>(tup);
  auto eigenvectors = std::get<1>(tup);

  // place eigenvalues in an VectorBTE object
  Eigen::VectorXd eigenValues(eigenvalues.size());
  for (int is = 0; is < eigenValues.size(); is++) {
    eigenValues(is) = eigenvalues[is];
  }

  return std::make_tuple(eigenValues, eigenvectors);
}

std::vector<std::tuple<std::vector<int>, int>>
ScatteringMatrix::getIteratorWavevectorPairs(const int &switchCase,
                                             const bool &rowMajor) {

  if (rowMajor) { // case for el-ph scattering
    std::vector<std::tuple<std::vector<int>, int>> pairIterator;

    if (switchCase == 1 || switchCase == 2) { // case for linewidth construction
      // here I parallelize over ik1
      // which is the outer loop on q-points
      std::vector<int> k1Iterator =
          outerBandStructure.parallelIrrPointsIterator();

      // I don't parallelize the inner band structure, the inner loop
      std::vector<int> k2Iterator(innerBandStructure.getNumPoints());
      // populate vector with integers from 0 to numPoints-1
      std::iota(std::begin(k2Iterator), std::end(k2Iterator), 0);

      for (int ik1 : k1Iterator) {
        auto t = std::make_tuple(k2Iterator, ik1);
        pairIterator.push_back(t);
      }
    } else { // case el-ph scattering, matrix in memory

      // here we operate assuming innerBandStructure=outerBandStructure
      // list in form [[0,0],[1,0],[2,0],...]
      std::set<std::pair<int, int>> localPairs;
      std::vector<std::tuple<int, int>> localStates = theMatrix.getAllLocalStates();
      int nlocalStates = localStates.size();
#pragma omp parallel default(none) shared(localPairs, theMatrix, localStates, nlocalStates)
      {
        std::vector<std::pair<int, int>> localPairsPrivate;
        // first unpack Bloch Index and get all wavevector pairs
#pragma omp for nowait
        for(int ilocalState = 0; ilocalState < nlocalStates; ilocalState++){
          std::tuple<int,int> tup0 = localStates[ilocalState];
          auto iMat1 = std::get<0>(tup0);
          auto iMat2 = std::get<1>(tup0);
          auto tup1 = getSMatrixIndex(iMat1);
          auto tup2 = getSMatrixIndex(iMat2);
          BteIndex iBte1 = std::get<0>(tup1);
          BteIndex iBte2 = std::get<0>(tup2);
          // map the index on the irr points of BTE to band structure index
          StateIndex is1 = outerBandStructure.bteToState(iBte1);
          StateIndex is2 = outerBandStructure.bteToState(iBte2);
          auto tuple1 = outerBandStructure.getIndex(is1);
          auto tuple2 = outerBandStructure.getIndex(is2);
          WavevectorIndex ik1Index = std::get<0>(tuple1);
          WavevectorIndex ik2Index = std::get<0>(tuple2);
          int ik1Irr = ik1Index.get();
          int ik2Irr = ik2Index.get();
          for (int ik2 : outerBandStructure.getReducibleStarFromIrreducible(ik2Irr)) {
            std::pair<int, int> xx = std::make_pair(ik1Irr, ik2);
            localPairsPrivate.push_back(xx);
          }
        }
#pragma omp critical
        for (auto y : localPairsPrivate) {
          localPairs.insert(y);
        }
      }

      // find set of q1
      // this because we have duplicates with different band indices.
      std::set<int> q1IndexesSet;
      for (std::pair<int, int> p : localPairs) {
        auto iq1 = p.first;
        q1IndexesSet.insert(iq1);
      }

      // we move the set into a vector
      // vector, unlike set, has well-defined indices for its elements
      std::vector<int> q1Indexes;
      for (int x : q1IndexesSet) {
        q1Indexes.push_back(x);
      }

      for (int iq1 : q1Indexes) {
        std::vector<int> x;
        auto y = std::make_tuple(x, iq1);
        pairIterator.push_back(y);
      }

      for (std::pair<int, int> p : localPairs) {
        int iq1 = p.first;
        int iq2 = p.second;

        auto idx = std::find(q1Indexes.begin(), q1Indexes.end(), iq1);
        if (idx != q1Indexes.end()) {
          auto i = idx - q1Indexes.begin();
          // note: get<> returns a reference to the tuple elements
          std::get<0>(pairIterator[i]).push_back(iq2);
        } else {
          Error("iq1 not found, not supposed to happen");
        }
      }
    }

    // Patch to make pooled interpolation of coupling work
    // we need to make sure each MPI process in the pool calls
    // the calculation of the coupling interpolation the same number of times
    // Hence, we correct pairIterator so that it does have the same size across
    // the pool

    if (mpi->getSize(mpi->intraPoolComm) > 1) {
      // std::vector<std::tuple<std::vector<int>, int>> pairIterator;
      auto myNumK1 = int(pairIterator.size());
      int numK1 = myNumK1;
      mpi->allReduceMax(&numK1, mpi->intraPoolComm);

      while (myNumK1 < numK1) {
        std::vector<int> dummyVec;
        dummyVec.push_back(-1);
        auto tt = std::make_tuple(dummyVec,-1);
        pairIterator.push_back(tt);
        myNumK1++;
      }
    }
    return pairIterator;

  } else { // case for ph_scattering

    if (switchCase == 1 || switchCase == 2) { // case for dot
      // must parallelize over the inner band structure (iq2 in phonons)
      // which is the outer loop on q-points
      size_t a = innerBandStructure.getNumPoints();
      auto q2Iterator = mpi->divideWorkIter(a);
      std::vector<int> q1Iterator = outerBandStructure.irrPointsIterator();

      std::vector<std::tuple<std::vector<int>, int>> pairIterator;
      for (int iq2 : q2Iterator) {
        auto t = std::make_tuple(q1Iterator, iq2);
        pairIterator.push_back(t);
      }

      return pairIterator;

    } else { // case for constructing A matrix

      std::set<std::pair<int, int>> localPairs;
      std::vector<std::tuple<int, int>> localStates = theMatrix.getAllLocalStates();
      int nlocalStates = localStates.size();
#pragma omp parallel default(none) shared(theMatrix, localPairs, localStates, nlocalStates)
      {
        std::vector<std::pair<int, int>> localPairsPrivate;
        // first unpack Bloch Index and get all wavevector pairs
#pragma omp for nowait
        for(int ilocalState = 0; ilocalState < nlocalStates; ilocalState++){
          std::tuple<int,int> tup0 = localStates[ilocalState];
          auto iMat1 = std::get<0>(tup0);
          auto iMat2 = std::get<1>(tup0);
          auto tup1 = getSMatrixIndex(iMat1);
          auto tup2 = getSMatrixIndex(iMat2);
          BteIndex iBte1 = std::get<0>(tup1);
          BteIndex iBte2 = std::get<0>(tup2);
          // map the index on the irreducible points of BTE to band structure
          // index
          StateIndex is1 = outerBandStructure.bteToState(iBte1);
          StateIndex is2 = outerBandStructure.bteToState(iBte2);
          auto tuple1 = outerBandStructure.getIndex(is1);
          auto tuple2 = outerBandStructure.getIndex(is2);
          WavevectorIndex iq1Index = std::get<0>(tuple1);
          WavevectorIndex iq2Index = std::get<0>(tuple2);
          int iq1Irr = iq1Index.get();
          int iq2Irr = iq2Index.get();
          for (int iq2 :
               outerBandStructure.getReducibleStarFromIrreducible(iq2Irr)) {
            std::pair<int, int> xx = std::make_pair(iq1Irr, iq2);
            localPairsPrivate.push_back(xx);
          }
        }
#pragma omp critical
        for (auto y : localPairsPrivate) {
          localPairs.insert(y);
        }
      }

      // find set of q2
      std::set<int> q2IndexesSet;
      for (auto tup : localPairs) {
        auto iq2 = std::get<1>(tup);
        q2IndexesSet.insert(iq2);
      }

      // we move the set into a vector
      // vector, unlike set, has well-defined indices for its elements
      std::vector<int> q2Indexes;
      for (int x : q2IndexesSet) {
        q2Indexes.push_back(x);
      }

      std::vector<std::tuple<std::vector<int>, int>> pairIterator;
      for (int iq2 : q2Indexes) {
        std::vector<int> x;
        auto y = std::make_tuple(x, iq2);
        pairIterator.push_back(y);
      }

      for (std::pair<int, int> p : localPairs) {
        int iq1 = p.first;
        int iq2 = p.second;

        auto idx = std::find(q2Indexes.begin(), q2Indexes.end(), iq2);
        if (idx != q2Indexes.end()) {
          auto i = idx - q2Indexes.begin();
          // note: get<> returns a reference to the tuple elements
          std::get<0>(pairIterator[i]).push_back(iq1);
        } else {
          Error("iq2 not found, not supposed to happen");
        }
      }
      return pairIterator;
    }
  }
}

int ScatteringMatrix::getSMatrixIndex(BteIndex &bteIndex, CartIndex &cartIndex) {
  if (context.getUseSymmetries()) {
    int is = bteIndex.get();
    int alpha = cartIndex.get();
    return compress2Indices(is, alpha, numStates, 3);
  } else {
    return bteIndex.get();
  }
}

std::tuple<BteIndex, CartIndex> ScatteringMatrix::getSMatrixIndex(const int &iMat) {
  if (context.getUseSymmetries()) {
    auto t = decompress2Indices(iMat, numStates, 3);
    return std::make_tuple(BteIndex(std::get<0>(t)), CartIndex(std::get<1>(t)));
  } else {
    return std::make_tuple(BteIndex(iMat), CartIndex(0));
  }
}

void ScatteringMatrix::symmetrize() {
  // note: if the matrix is not stored in memory, it's not trivial to enforce
  // that the matrix is symmetric

  if (context.getUseSymmetries()) {
    Error("Symmetrization of the scattering matrix with symmetries"
          " is not verified");
    // not only, probably it doesn't make sense. To be checked
  }

  if (highMemory) {
    ParallelMatrix<double> newMatrix = theMatrix;
    newMatrix *= 0.;

    for (int iRank = 0; iRank < mpi->getSize(); iRank++) {
      int thisSize = 0;
      if (iRank == mpi->getRank()) {
        thisSize = int(theMatrix.getAllLocalStates().size());
      }
      mpi->allReduceSum(&thisSize);

      std::vector<int> index1(thisSize, 0);
      std::vector<int> index2(thisSize, 0);
      std::vector<int> index3(thisSize, 0);
      std::vector<int> index4(thisSize, 0);

      if (iRank == mpi->getRank()) {
        auto allLocalStates = theMatrix.getAllLocalStates();
        size_t numAllLocalStates = allLocalStates.size();

        for (size_t i=0; i<numAllLocalStates; i++) {
          auto tup = allLocalStates[i];
          int iMat1 = std::get<0>(tup);
          int iMat2 = std::get<1>(tup);
          int jMat1, jMat2;
//          if (context.getUseSymmetries()) {
//            Error("Symmetrization of the scattering matrix with symmetries"
//                  " is not verified");
//            // The matrix may not be symmetric
//
//            auto t1 = getSMatrixIndex(iMat1);
//            auto t2 = getSMatrixIndex(iMat2);
//            BteIndex iBte1 = std::get<0>(t1);
//            BteIndex iBte2 = std::get<0>(t2);
//            CartIndex ii = std::get<1>(t1);
//            CartIndex jj = std::get<1>(t2);
//            jMat1 = getSMatrixIndex(iBte1, ii);
//            jMat2 = getSMatrixIndex(iBte2, jj);
//          } else {
          jMat1 = iMat1;
          jMat2 = iMat2;
//          }

          index1[i] = iMat1;
          index2[i] = iMat2;
          index3[i] = jMat2;
          index4[i] = jMat1;
        }
      }

      mpi->allReduceSum(&index1);
      mpi->allReduceSum(&index2);
      mpi->allReduceSum(&index3);
      mpi->allReduceSum(&index4);

      std::vector<double> values(thisSize, 0.);
#pragma omp parallel for
      for (int i = 0; i < thisSize; i++) {
        values[i] = (theMatrix(index1[i], index2[i]) +
                     theMatrix(index3[i], index4[i])) *
                    0.5;
      }
      mpi->allReduceSum(&values);

#pragma omp parallel for
      for (int i = 0; i < thisSize; i++) {
        newMatrix(index1[i], index2[i]) = values[i];
        newMatrix(index3[i], index4[i]) = values[i];
      }
    }
    theMatrix = newMatrix;
  }
}

void ScatteringMatrix::degeneracyAveragingLinewidths(std::shared_ptr<VectorBTE> linewidth) {
  for (int ik : outerBandStructure.irrPointsIterator()) {
    WavevectorIndex ikIdx(ik);
    Eigen::VectorXd en = outerBandStructure.getEnergies(ikIdx);
    auto numBands = int(en.size());

    Eigen::VectorXd linewidthTmp(numBands);
    linewidthTmp.setZero();

    for (int ib1 = 0; ib1 < numBands; ib1++) {
      double ekk = en(ib1);
      int n = 0;
      double tmp2 = 0.;
      for (int ib2 = 0; ib2 < numBands; ib2++) {
        double ekk2 = en(ib2);

        BandIndex ibIdx(ib2);
        int is = outerBandStructure.getIndex(ikIdx, ibIdx);
        StateIndex isIdx(is);
        BteIndex ind1Idx = outerBandStructure.stateToBte(isIdx);
        int iBte1 = ind1Idx.get();

        if (abs(ekk2 - ekk) < 1.0e-6) {
          n++;
          tmp2 = tmp2 + linewidth->data(0, iBte1);
        }
      }
      linewidthTmp(ib1) = tmp2 / float(n);
    }

    for (int ib1 = 0; ib1 < numBands; ib1++) {
      BandIndex ibIdx(ib1);
      int is = outerBandStructure.getIndex(ikIdx, ibIdx);
      StateIndex isIdx(is);
      BteIndex ind1Idx = outerBandStructure.stateToBte(isIdx);
      int iBte1 = ind1Idx.get();
      linewidth->data(0, iBte1) = linewidthTmp(ib1);
    }
  }
}

void ScatteringMatrix::symmetrizeCoupling(Eigen::Tensor<double,3>& coupling,
                                          const Eigen::VectorXd& energies1,
                                          const Eigen::VectorXd& energies2,
                                          const Eigen::VectorXd& energies3){
  int nb1 = energies1.size();
  int nb2 = energies2.size();
  int nb3 = energies3.size();

  for (int ib1 = 0; ib1 < nb1; ib1++) {
    double e1 = energies1(ib1);
    // determine degeneracy degree
    int degDegree = 0;
    for (int i = ib1; i < nb1; i++) {
      double e2 = energies1(i);
      if (abs(e1 - e2) < 1.0e-6) { // at first iteration, ib1=ib2, degDegree=1
        degDegree++;
      } else {
        break;
      }
    }
    Eigen::MatrixXd tmpCoupling(nb2,nb3);
    tmpCoupling.setZero();
    // now do averaging
    for (int i=0; i<degDegree; ++i) {
      for (int ib2 = 0; ib2 < nb2; ib2++) {
        for (int ib3 = 0; ib3 < nb3; ib3++) {
          tmpCoupling(ib2, ib3) += coupling(ib1 + i, ib2, ib3) / float(degDegree);
        }
      }
    }
    // substitute back
    for (int i=0; i<degDegree; ++i) {
      for (int ib2 = 0; ib2 < nb2; ib2++) {
        for (int ib3 = 0; ib3 < nb3; ib3++) {
          coupling(ib1 + i, ib2, ib3) = tmpCoupling(ib2, ib3);
        }
      }
    }
    // now skip to next iteration
    ib1 += degDegree - 1; // -1 because there's another addition in the loop
  }

  for (int ib2 = 0; ib2 < nb2; ib2++) {
    double e1 = energies2(ib2);
    // determine degeneracy degree
    int degDegree = 0;
    for (int i = ib2; i < nb2; i++) {
      double e2 = energies2(i);
      if (abs(e1 - e2) < 1.0e-6) { // at first iteration, ib1=ib2, degDegree=1
        degDegree++;
      } else {
        break;
      }
    }
    Eigen::MatrixXd tmpCoupling(nb1,nb3);
    tmpCoupling.setZero();
    // now do averaging
    for (int i=0; i<degDegree; ++i) {
      for (int ib1 = 0; ib1 < nb1; ib1++) {
        for (int ib3 = 0; ib3 < nb3; ib3++) {
          tmpCoupling(ib1, ib3) += coupling(ib1, ib2+i, ib3) / float(degDegree);
        }
      }
    }
    // substitute back
    for (int i=0; i<degDegree; ++i) {
      for (int ib1 = 0; ib1 < nb1; ib1++) {
        for (int ib3 = 0; ib3 < nb3; ib3++) {
          coupling(ib1, ib2 + i, ib3) = tmpCoupling(ib1, ib3);
        }
      }
    }
    // now skip to next iteration
    ib2 += degDegree - 1; // -1 because there's another addition in the loop
  }

  for (int ib3 = 0; ib3 < nb3; ib3++) {
    double e1 = energies3(ib3);
    // determine degeneracy degree
    int degDegree = 0;
    for (int i = ib3; i < nb3; i++) {
      double e2 = energies3(i);
      if (abs(e1 - e2) < 1.0e-6) { // at first iteration, ib1=ib2, degDegree=1
        degDegree++;
      } else {
        break;
      }
    }
    Eigen::MatrixXd tmpCoupling(nb1,nb2);
    tmpCoupling.setZero();
    // now do averaging
    for (int i=0; i<degDegree; ++i) {
      for (int ib1 = 0; ib1 < nb1; ib1++) {
        for (int ib2 = 0; ib2 < nb2; ib2++) {
          tmpCoupling(ib1, ib2) += coupling(ib1, ib2, ib3+i) / float(degDegree);
        }
      }
    }
    // substitute back
    for (int i=0; i<degDegree; ++i) {
      for (int ib1 = 0; ib1 < nb1; ib1++) {
        for (int ib2 = 0; ib2 < nb2; ib2++) {
          coupling(ib1, ib2, ib3 + i) = tmpCoupling(ib1, ib2);
        }
      }
    }
    // now skip to next iteration
    ib3 += degDegree - 1; // -1 because there's another addition in the loop
  }
}

// helper to precompute particle occupations
Eigen::MatrixXd ScatteringMatrix::precomputeOccupations(BaseBandStructure &bandStructure) {

  int numCalculations = statisticsSweep.getNumCalculations();
  auto numIrrStates = int(bandStructure.irrStateIterator().size());
  Particle particle = bandStructure.getParticle();

  Eigen::MatrixXd occupations(numCalculations, numIrrStates);
  occupations.setZero();
  std::vector<size_t> iBtes = mpi->divideWorkIter(numIrrStates);
  size_t niBtes = iBtes.size();
#pragma omp parallel for default(none)                                         \
    shared(mpi, particle, occupations, numIrrStates, numCalculations, niBtes, iBtes, bandStructure)
  for(size_t iiBte = 0; iiBte < niBtes; iiBte++){
    int iBte = iBtes[iiBte];
    BteIndex iBteIdx(iBte);
    StateIndex isIdx = bandStructure.bteToState(iBteIdx);
    double energy = bandStructure.getEnergy(isIdx);
    for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      double temp = calcStat.temperature;
      // need to do this because even if it's an el bandstructure
      // you could need ph occupations, and statSweep would have a non-zero mu
      double chemPot = (particle.isPhonon()) ? 0 : calcStat.chemicalPotential;
      occupations(iCalc, iBte) = particle.getPopulation(energy, temp, chemPot);
    }
  }
  mpi->allReduceSum(&occupations);
  return occupations;
}

std::vector<int> ScatteringMatrix::getExcludeIndices(BaseBandStructure& bandStructure) {
  // we want to know the state index of acoustic modes at gamma,
  // so that we can set their populations to zero
  std::vector<int> indices;
  int bsNumStates = int(bandStructure.irrStateIterator().size());
  if (bandStructure.getParticle().isPhonon()) {
    for (int iBte = 0; iBte < bsNumStates; iBte++) {
      auto iBteIdx = BteIndex(iBte);
      StateIndex isIdx = bandStructure.bteToState(iBteIdx);
      double en = bandStructure.getEnergy(isIdx);
      // TODO need to standardize this
      if (en < 0.1 / ryToCmm1) { // cutoff at 0.1 cm^-1
        indices.push_back(iBte);
      }

      Eigen::Vector3d k = bandStructure.getWavevector(isIdx);
      if (k.squaredNorm() > 1e-8 && en < 0.) {
        Warning("Found a phonon mode q!=0 with negative energies."
                "omega = " + std::to_string(en) +
                "k = " + std::to_string(k[0]) + std::to_string(k[1]) + std::to_string(k[2]) +
                "Consider improving the quality of your DFT phonon calculation.\n");
      }
    }
  }
  return indices;
}

/* probably never called on regular Scattering matrix, but this is here just in case */
std::tuple<int,int> ScatteringMatrix::shiftToCoupledIndices(
            const int& iBte1, const int& iBte2, const Particle& p1, const Particle& p2) {

  // we shift the iBte indices into the relevant quadrant before savign things to the
  // scattering matrix in the coupled case
  // ibte1 = row, ibte2 = col

  int iBte1Shift = iBte1;
  int iBte2Shift = iBte2;
  if(isCoupled) {
    if(!p1.isPhonon() && p2.isPhonon()) { // electron-phonon drag
      iBte2Shift += numElStates;
    }
    else if(p1.isPhonon() && !p2.isPhonon()) { // phonon-electron drag
      iBte1Shift += numElStates;
    }
    else if(p1.isPhonon() && p2.isPhonon()) { // phonon-self quadrant
      iBte1Shift += numElStates;
      iBte2Shift += numElStates;
    }
    return std::make_tuple(iBte1Shift, iBte2Shift);
  }
  Error("Developer error: there's no reason to shift regular scattering matrix indices!");
  return std::make_tuple(iBte1Shift, iBte2Shift);
}


