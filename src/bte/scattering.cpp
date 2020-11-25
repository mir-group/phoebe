#include "scattering.h"
#include "constants.h"
#include "mpiHelper.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <numeric> // std::iota
#include <set>

ScatteringMatrix::ScatteringMatrix(Context &context_,
                                   StatisticsSweep &statisticsSweep_,
                                   BaseBandStructure &innerBandStructure_,
                                   BaseBandStructure &outerBandStructure_)
    : context(context_), statisticsSweep(statisticsSweep_),
      innerBandStructure(innerBandStructure_),
      outerBandStructure(outerBandStructure_),
      internalDiagonal(statisticsSweep, outerBandStructure, 1) {
  numStates = outerBandStructure.irrStateIterator().size();
  numPoints = outerBandStructure.irrPointsIterator().size();

  double constantRelaxationTime = context.getConstantRelaxationTime();
  if (constantRelaxationTime > 0.) {
    constantRTA = true;
    return;
  }

  dimensionality_ = context.getDimensionality();

  highMemory = context.getScatteringMatrixInMemory();

  smearing = DeltaFunction::smearingFactory(context, innerBandStructure);

  if ( // innerBandStructure != outerBandStructure &&
      smearing->getType() == smearing->tetrahedron) {
    Error e("Tetrahedron smearing for transport untested and thus blocked");
    // not for linewidths. Although this should be double-checked
  }

  numCalcs = statisticsSweep.getNumCalcs();

  // we want to know the state index of acoustic modes at gamma,
  // so that we can set their populations to zero
  if (outerBandStructure.getParticle().isPhonon()) {
    for (long ibte = 0; ibte < numStates; ibte++) {
      auto ibteIdx = BteIndex(ibte);
      long is = outerBandStructure.bteToState(ibteIdx).get();
      double en = outerBandStructure.getEnergy(is);
      if (en < 0.1 / ryToCmm1) { // cutoff at 0.1 cm^-1
        excludeIndeces.push_back(ibte);
      }
    }
  }
}

ScatteringMatrix::~ScatteringMatrix() { delete smearing; }

// copy constructor
ScatteringMatrix::ScatteringMatrix(const ScatteringMatrix &that)
    : context(that.context), statisticsSweep(that.statisticsSweep),
      smearing(that.smearing), innerBandStructure(that.innerBandStructure),
      outerBandStructure(that.outerBandStructure),
      constantRTA(that.constantRTA), highMemory(that.highMemory),
      internalDiagonal(that.internalDiagonal), theMatrix(that.theMatrix),
      numStates(that.numStates), numPoints(that.numPoints),
      numCalcs(that.numCalcs), dimensionality_(that.dimensionality_),
      excludeIndeces(that.excludeIndeces) {}

// assignment operator
ScatteringMatrix &ScatteringMatrix::operator=(const ScatteringMatrix &that) {
  if (this != &that) {
    context = that.context;
    statisticsSweep = that.statisticsSweep;
    smearing = that.smearing;
    innerBandStructure = that.innerBandStructure;
    outerBandStructure = that.outerBandStructure;
    constantRTA = that.constantRTA;
    highMemory = that.highMemory;
    internalDiagonal = that.internalDiagonal;
    theMatrix = that.theMatrix;
    numStates = that.numStates;
    numPoints = that.numPoints;
    numCalcs = that.numCalcs;
    excludeIndeces = that.excludeIndeces;
    dimensionality_ = that.dimensionality_;
  }
  return *this;
}

void ScatteringMatrix::setup() {
  // note: here we want to build the matrix or its diagonal
  // builder is a pure virtual function, which is implemented in subclasses
  // c++ discourages calls to pure virtual functions in the constructor

  if (constantRTA)
    return; // nothing to construct

  std::vector<VectorBTE> emptyVector;

  if (highMemory) {
    if (numCalcs > 1) {
      // note: one could write code around this
      // but the methods are very memory intensive for production runs
      Error e("High memory BTE methods can only work with one "
              "temperature and/or chemical potential in a single run");
    }
    if (context.getUseSymmetries()) {
      theMatrix = ParallelMatrix<double>(3*numStates, 3*numStates);
    } else {
      theMatrix = ParallelMatrix<double>(numStates, numStates);
    }

    // calc matrix and linew.
    builder(&internalDiagonal, emptyVector, emptyVector);
  } else {
    // calc linewidths only
    builder(&internalDiagonal, emptyVector, emptyVector);
  }
}

VectorBTE ScatteringMatrix::diagonal() {
  if (constantRTA) {
    double crt = context.getConstantRelaxationTime();
    VectorBTE diag(statisticsSweep, outerBandStructure, 1);
    diag.setConst(1. / crt);
    return diag;
  } else {
    return internalDiagonal;
  }
}

VectorBTE ScatteringMatrix::offDiagonalDot(VectorBTE &inPopulation) {
  // outPopulation = outPopulation - internalDiagonal * inPopulation;
  VectorBTE outPopulation = dot(inPopulation);
#pragma omp parallel for collapse(3)
  for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
    for (int iDim : {0, 1, 2}) {
      for (long ibte = 0; ibte < numStates; ibte++) {
        outPopulation(iCalc, iDim, ibte) -=
            internalDiagonal(iCalc, 0, ibte) * inPopulation(iCalc, iDim, ibte);
      }
    }
  }
  return outPopulation;
}

std::vector<VectorBTE>
ScatteringMatrix::offDiagonalDot(std::vector<VectorBTE> &inPopulations) {
  // outPopulation = outPopulation - internalDiagonal * inPopulation;
  std::vector<VectorBTE> outPopulations = dot(inPopulations);
  for (unsigned int iVec = 0; iVec < inPopulations.size(); iVec++) {
#pragma omp parallel for collapse(3)
    for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
      for (int iDim : {0, 1, 2}) {
        for (long ibte = 0; ibte < numStates; ibte++) {
          outPopulations[iVec](iCalc, iDim, ibte) -=
              internalDiagonal(iCalc, 0, ibte) *
              inPopulations[iVec](iCalc, iDim, ibte);
        }
      }
    }
  }
  return outPopulations;
}

VectorBTE ScatteringMatrix::dot(VectorBTE &inPopulation) {
  if (highMemory) {
    VectorBTE outPopulation(statisticsSweep, outerBandStructure,
                            inPopulation.dimensionality);
    // note: we are assuming that ScatteringMatrix has numCalcs = 1

    if (context.getUseSymmetries()) {
      for (auto tup : theMatrix.getAllLocalStates()) {
        long iMat1 = std::get<0>(tup);
        long iMat2 = std::get<1>(tup);
        auto t1 = getSMatrixIndex(iMat1);
        auto t2 = getSMatrixIndex(iMat2);
        long ibte1 = std::get<0>(t1).get();
        long ibte2 = std::get<0>(t2).get();
        long i = std::get<1>(t1).get();
        long j = std::get<1>(t2).get();
        outPopulation(0, i, ibte1) +=
            theMatrix(iMat1, iMat2) * inPopulation(0, j, ibte2);
      }
    } else {
      for (auto tup : theMatrix.getAllLocalStates()) {
        auto ibte1 = std::get<0>(tup);
        auto ibte2 = std::get<1>(tup);
        for (int i : {0, 1, 2}) {
          outPopulation(0, i, ibte1) +=
              theMatrix(ibte1, ibte2) * inPopulation(0, i, ibte2);
        }
      }
    }

    mpi->allReduceSum(&outPopulation.data);
    return outPopulation;
  } else {
    VectorBTE outPopulation(statisticsSweep, outerBandStructure,
                            inPopulation.dimensionality);
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
      VectorBTE outPopulation(statisticsSweep, outerBandStructure,
                              inPopulation.dimensionality);
      outPopulation = dot(inPopulation);
      outPopulations.push_back(outPopulation);
    }
    return outPopulations;
  } else {
    std::vector<VectorBTE> outPopulations;
    for (auto inPopulation : inPopulations) {
      VectorBTE outPopulation(statisticsSweep, outerBandStructure,
                              inPopulation.dimensionality);
      outPopulations.push_back(outPopulation);
    }
    builder(nullptr, inPopulations, outPopulations);
    return outPopulations;
  }
}

// set,unset the scaling of omega = A/sqrt(bose1*bose1+1)/sqrt(bose2*bose2+1)
void ScatteringMatrix::a2Omega() {
  if (!highMemory) {
    Error e("a2Omega only works if the matrix is stored in memory");
  }

  if (theMatrix.rows() == 0) {
    Error e("The scattering matrix hasn't been built yet");
  }

  if (isMatrixOmega) { // it's already with the scaling of omega
    return;
  }

  long iCalc = 0; // as there can only be one temperature

  auto particle = outerBandStructure.getParticle();
  auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
  double temp = calcStatistics.temperature;
  double chemPot = calcStatistics.chemicalPotential;

  for (auto tup : theMatrix.getAllLocalStates()) {
    long ibte1, ibte2, iMat1, iMat2, is1, is2;
    if (context.getUseSymmetries()) {
      iMat1 = std::get<0>(tup);
      iMat2 = std::get<1>(tup);
      auto tup1 = getSMatrixIndex(iMat1);
      auto tup2 = getSMatrixIndex(iMat2);
      BteIndex ibte1Idx = std::get<0>(tup1);
      BteIndex ibte2Idx = std::get<0>(tup2);
      ibte1 = ibte1Idx.get();
      ibte2 = ibte2Idx.get();
      is1 = outerBandStructure.bteToState(ibte1Idx).get();
      is2 = outerBandStructure.bteToState(ibte2Idx).get();
    } else {
      ibte1 = std::get<0>(tup);
      ibte2 = std::get<1>(tup);
      iMat1 = ibte1;
      iMat2 = ibte2;
      is1 = ibte1;
      is2 = ibte2;
    }
    if (std::find(excludeIndeces.begin(), excludeIndeces.end(), ibte1) !=
        excludeIndeces.end())
      continue;
    if (std::find(excludeIndeces.begin(), excludeIndeces.end(), ibte2) !=
        excludeIndeces.end())
      continue;

    double en1 = outerBandStructure.getEnergy(is1);
    double en2 = outerBandStructure.getEnergy(is2);

    // n(n+1) for bosons, n(1-n) for fermions
    double term1 = particle.getPopPopPm1(en1, temp, chemPot);
    double term2 = particle.getPopPopPm1(en2, temp, chemPot);

    if (ibte1 == ibte2) {
      internalDiagonal(0, 0, ibte1) /= term1;
    }

    theMatrix(iMat1, iMat2) /= sqrt(term1 * term2);
  }
  isMatrixOmega = true;
}

// add a flag to remember if we have A or Omega
//void ScatteringMatrix::omega2A() {
//  if (!highMemory) {
//    Error e("a2Omega only works if the matrix is stored in memory");
//  }
//
//  if (theMatrix.rows() == 0) {
//    Error e("The scattering matrix hasn't been built yet");
//  }
//
//  if (!isMatrixOmega) { // it's already with the scaling of A
//    return;
//  }
//
//  long iCalc = 0; // as there can only be one temperature
//
//  auto particle = outerBandStructure.getParticle();
//  auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
//  double temp = calcStatistics.temperature;
//  double chemPot = calcStatistics.chemicalPotential;
//
//  for (auto tup : theMatrix.getAllLocalStates()) {
//    auto ind1 = std::get<0>(tup);
//    auto ind2 = std::get<1>(tup);
//    if (std::find(excludeIndeces.begin(), excludeIndeces.end(), ind1) !=
//        excludeIndeces.end())
//      continue;
//    if (std::find(excludeIndeces.begin(), excludeIndeces.end(), ind2) !=
//        excludeIndeces.end())
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
VectorBTE ScatteringMatrix::getSingleModeTimes() {
  if (constantRTA) {
    double crt = context.getConstantRelaxationTime();
    VectorBTE times(statisticsSweep, outerBandStructure, 1);
    times.setConst(crt);
    times.excludeIndeces = excludeIndeces;
    return times;
  } else {
    if (isMatrixOmega) {
      VectorBTE times = internalDiagonal.reciprocal();
      times.excludeIndeces = excludeIndeces;
      return times;
    } else { // A_nu,nu = N(1+-N) / tau
      VectorBTE times = internalDiagonal;
      auto particle = outerBandStructure.getParticle();
      for (long ibte = 0; ibte < numStates; ibte++) {
        auto ibteIdx = BteIndex(ibte);
        long is = outerBandStructure.bteToState(ibteIdx).get();
        double en = outerBandStructure.getEnergy(is);
        for (long iCalc = 0; iCalc < internalDiagonal.numCalcs; iCalc++) {
          auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
          double temp = calcStatistics.temperature;
          double chemPot = calcStatistics.chemicalPotential;
          // n(n+1) for bosons, n(1-n) for fermions
          double popTerm = particle.getPopPopPm1(en, temp, chemPot);
          times(iCalc, 0, ibte) = popTerm / internalDiagonal(iCalc, 0, ibte);
        }
      }
      times.excludeIndeces = excludeIndeces;
      return times;
    }
  }
}

// to compute the RTA, get the single mode relaxation times
VectorBTE ScatteringMatrix::getLinewidths() {
  if (constantRTA) {
    double crt = context.getConstantRelaxationTime();
    VectorBTE linewidths(statisticsSweep, outerBandStructure, 1);
    linewidths.setConst(1. / crt);
    linewidths.excludeIndeces = excludeIndeces;
    return linewidths;
  } else {
    VectorBTE linewidths = internalDiagonal;
    linewidths.excludeIndeces = excludeIndeces;
    if (isMatrixOmega) {
      // Important note: don't use this for fermions, or you can get rubbish!
      // the factor popTerm could be = 0!
      auto particle = outerBandStructure.getParticle();
      if (particle.isElectron()) {
        Error e("Attempting to use a numerically unstable quantity");
      }
      for (long ibte = 0; ibte < numStates; ibte++) {
        auto ibteIdx = BteIndex(ibte);
        long is = outerBandStructure.bteToState(ibteIdx).get();
        double en = outerBandStructure.getEnergy(is);
        for (long iCalc = 0; iCalc < internalDiagonal.numCalcs; iCalc++) {
          auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
          double temp = calcStatistics.temperature;
          double chemPot = calcStatistics.chemicalPotential;
          // n(n+1) for bosons, n(1-n) for fermions
          double popTerm = particle.getPopPopPm1(en, temp, chemPot);
          linewidths(iCalc, 0, ibte) /= popTerm;
        }
      }
      return linewidths;

    } else { // A_nu,nu = N(1+-N) / tau
      return linewidths;
    }
  }
}

void ScatteringMatrix::outputToJSON(std::string outFileName) {

  if (!mpi->mpiHead())
    return;

  VectorBTE times = getSingleModeTimes();
  VectorBTE tmpLinewidths = getLinewidths();

  std::string energyUnit;
  std::string particleType;
  double energyConversion;
  auto particle = outerBandStructure.getParticle();
  energyConversion = energyRyToEv;
  energyUnit = "eV";
  if (particle.isPhonon()) {
    particleType = "phonon";
  } else {
    particleType = "electron";
  }

  // need to store as a vector format with dimensions
  // icalc, ik. ib, idim (where istate is unfolded into
  // ik, ib) for the velocities and lifetimes, no dim for energies
  std::vector<std::vector<std::vector<double>>> outTimes;
  std::vector<std::vector<std::vector<double>>> outLinewidths;
  std::vector<std::vector<std::vector<std::vector<double>>>> velocities;
  std::vector<std::vector<std::vector<double>>> energies;
  std::vector<double> temps;
  std::vector<double> chemPots;

  for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
    auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStatistics.temperature;
    double chemPot = calcStatistics.chemicalPotential;
    temps.push_back(temp * temperatureAuToSi);
    chemPots.push_back(chemPot * energyConversion);

    std::vector<std::vector<double>> wavevectorsT;
    std::vector<std::vector<double>> wavevectorsL;
    std::vector<std::vector<std::vector<double>>> wavevectorsV;
    std::vector<std::vector<double>> wavevectorsE;
    // loop over wavevectors
    for (long ik : outerBandStructure.irrPointsIterator()) {
      auto ikIndex = WavevectorIndex(ik);

      std::vector<double> bandsT;
      std::vector<double> bandsL;
      std::vector<std::vector<double>> bandsV;
      std::vector<double> bandsE;
      // loop over bands here
      // get numBands at this point, in case it's an active bandstructure
      for (int ib = 0; ib < outerBandStructure.getNumBands(ikIndex); ib++) {
        auto ibIndex = BandIndex(ib);
        long is = outerBandStructure.getIndex(ikIndex, ibIndex);
        double ene = outerBandStructure.getEnergy(is);
        auto vel = outerBandStructure.getGroupVelocity(is);
        bandsE.push_back(ene * energyConversion);
        auto isIdx = StateIndex(is);
        long ibte = outerBandStructure.stateToBte(isIdx).get();
        double tau = times(iCalc, 0, ibte); // only zero dim is meaningful
        bandsT.push_back(tau * timeRyToFs);
        double linewidth =
            tmpLinewidths(iCalc, 0, ibte); // only zero dim is meaningful
        bandsL.push_back(linewidth * energyRyToEv);

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

  auto points = outerBandStructure.getPoints();
  std::vector<std::vector<double>> meshCoords;
  for (long ik : outerBandStructure.irrPointsIterator()) {
    // save the wavevectors
    auto ikIndex = WavevectorIndex(ik);
    auto coord =
        points.cartesianToCrystal(outerBandStructure.getWavevector(ikIndex));
    meshCoords.push_back({coord[0], coord[1], coord[2]});
  }

  // output to json
  nlohmann::json output;
  output["temperatures"] = temps;
  output["temperatureUnit"] = "K";
  output["chemicalPotentials"] = chemPots;
  output["linewidths"] = outLinewidths;
  output["linewidthsUnit"] = "eV";
  output["relaxationTimes"] = outTimes;
  output["relaxationTimeUnit"] = "fs";
  output["velocities"] = velocities;
  output["velocityUnit"] = "m/s";
  output["energies"] = energies;
  output["energyUnit"] = energyUnit;
  output["wavevectorCoordinates"] = meshCoords;
  output["coordsType"] = "lattice";
  output["particleType"] = particleType;
  std::ofstream o(outFileName);
  o << std::setw(3) << output << std::endl;
  o.close();
}

std::tuple<Eigen::VectorXd, ParallelMatrix<double>> ScatteringMatrix::diagonalize() {
  //    std::vector<double> eigenvalues;
  //    ParallelMatrix<double> eigenvectors;
  auto tup = theMatrix.diagonalize();
  auto eigenvalues = std::get<0>(tup);
  auto eigenvectors = std::get<1>(tup);

  // place eigenvalues in an VectorBTE object
  Eigen::VectorXd eigvals(theMatrix.rows());
  for (long is = 0; is < eigvals.size(); is++) {
    eigvals(is) = eigenvalues[is];
  }

  return {eigvals, eigenvectors};
}

std::vector<std::tuple<std::vector<long>, long>>
ScatteringMatrix::getIteratorWavevectorPairs(const int &switchCase,
                                             const bool &rowMajor) {
  if (rowMajor) { // case for el-ph scattering

    if (switchCase == 2) { // case for linewidth construction
      // here I parallelize over ik1
      // which is the outer loop on q-points
      size_t a = innerBandStructure.getNumPoints();
      std::vector<long> q2Iterator = mpi->divideWorkIter(a);

      // I don't parallelize the outer band structure (iq1 in phonons)
      // which is the inner loop on q-points
      std::vector<long> q1Iterator(outerBandStructure.getNumPoints());
      // populate vector with integers from 0 to numPoints-1
      std::iota(std::begin(q1Iterator), std::end(q1Iterator), 0);

      std::vector<std::tuple<std::vector<long>, long>> pairIterator(
          q2Iterator.size());
      int i = 0;
      for (long iq2 : q2Iterator) {
        auto t = std::make_tuple(q1Iterator, iq2);
        pairIterator[i] = t;
        i++;
      }
      return pairIterator;

    } else if (switchCase == 1) {
      // this case is used by el_scattering, to compute lifetimes, or A.f
      std::vector<long> k1Iterator = outerBandStructure.parallelIrrPointsIterator();

      // Note: phScatteringMatrix needs iq2 to be the outer loop
      // in order to be efficient!
      std::vector<long> k2Iterator(innerBandStructure.getNumPoints());
      // populate vector with integers from 0 to numPoints-1
      std::iota(std::begin(k2Iterator), std::end(k2Iterator), 0);

      std::vector<std::tuple<std::vector<long>, long>> pairIterator;
      for (long ik1 : k1Iterator) {
        auto t = std::make_tuple(k2Iterator, ik1);
        pairIterator.push_back(t);
      }
      return pairIterator;

    } else { // case el-ph scattering, matrix in memory
      // here we operate assuming innerBandStructure=outerBandStructure

      // list in form [[0,0],[1,0],[2,0],...]
      std::set<std::pair<int, int>> localPairs;
      // first unpack Bloch Index and get all wavevector pairs
      for (auto tup0 : theMatrix.getAllLocalStates()) {
        auto iMat1 = std::get<0>(tup0);
        auto iMat2 = std::get<1>(tup0);
        auto tup1 = getSMatrixIndex(iMat1);
        auto tup2 = getSMatrixIndex(iMat2);
        BteIndex ibte1 = std::get<0>(tup1);
        BteIndex ibte2 = std::get<0>(tup2);
        // map the index on the irreducible points of BTE to bandstructure index
        StateIndex is1 = outerBandStructure.bteToState(ibte1);
        StateIndex is2 = outerBandStructure.bteToState(ibte2);
        auto tupp1 = outerBandStructure.getIndex(is1);
        auto tupp2 = outerBandStructure.getIndex(is2);
        WavevectorIndex ik1Index = std::get<0>(tupp1);
        WavevectorIndex ik2Index = std::get<0>(tupp2);
        long ik1Irr = ik1Index.get();
        long ik2Irr = ik2Index.get();
        for ( long ik2 : outerBandStructure.getReduciblesFromIrreducible(ik2Irr) ) {
          std::pair<int, int> xx = std::make_pair(ik1Irr, ik2);
          localPairs.insert(xx);
        }
      }

      // find set of q1
      // this because we have duplicates with different band indices.
      std::set<long> q1Indexes;
      for (auto tup : localPairs) {
        auto iq1 = std::get<0>(tup);
        q1Indexes.insert(iq1);
      }

      std::vector<std::tuple<std::vector<long>, long>> pairIterator;
      // find all q2 for fixed q1
      for (long iq1 : q1Indexes) {
        std::vector<long> q2Indexes;

        // note: I'm assuming that the pairs in x are unique
        for (auto tup : localPairs) {
          auto iq1_ = std::get<0>(tup);
          auto iq2 = std::get<1>(tup);
          if (iq1 == iq1_) {
            q2Indexes.push_back(iq2);
          }
        }

        auto t = std::make_tuple(q2Indexes, iq1);
        pairIterator.push_back(t);
      }
      return pairIterator;
    }

  } else { // case for ph_scattering

    if (switchCase == 2) { // case for linewidth
      // must parallelize over the inner band structure (iq2 in phonons)
      // which is the outer loop on q-points
      size_t a = innerBandStructure.getNumPoints();
      std::vector<long> q2Iterator = mpi->divideWorkIter(a);

      // I don't parallelize the outer band structure (iq1 in phonons)
      // which is the inner loop on q-points
      std::vector<long> q1Iterator(outerBandStructure.getNumPoints());
      // populate vector with integers from 0 to numPoints-1
      std::iota(std::begin(q1Iterator), std::end(q1Iterator), 0);

      std::vector<std::tuple<std::vector<long>, long>> pairIterator(
          q2Iterator.size());
      int i = 0;
      for (long iq2 : q2Iterator) {
        auto t = std::make_tuple(q1Iterator, iq2);
        pairIterator[i] = t;
        i++;
      }
      return pairIterator;

    } else if (switchCase == 1) { // case for dot
      // must parallelize over the inner band structure (iq2 in phonons)
      // which is the outer loop on q-points
      size_t a = innerBandStructure.getNumPoints();
      std::vector<long> q2Iterator = mpi->divideWorkIter(a);
      std::vector<long> q1Iterator = outerBandStructure.irrPointsIterator();

      std::vector<std::tuple<std::vector<long>, long>> pairIterator(
          q2Iterator.size());
      int i = 0;
      for (long iq2 : q2Iterator) {
        auto t = std::make_tuple(q1Iterator, iq2);
        pairIterator[i] = t;
        i++;
      }
      return pairIterator;

    } else { // case for constructing A matrix

      // list in form [[0,0],[1,0],[2,0],...]
      std::set<std::pair<int, int>> x;
      auto points = outerBandStructure.getPoints();
      for (auto tup0 : theMatrix.getAllLocalStates()) {
        auto iMat1 = std::get<0>(tup0);
        auto iMat2 = std::get<1>(tup0);
        auto tup1 = getSMatrixIndex(iMat1);
        auto tup2 = getSMatrixIndex(iMat2);
        BteIndex ibte1 = std::get<0>(tup1);
        BteIndex ibte2 = std::get<0>(tup2);
        // map the index on the irreducible points of BTE to bandstructure index
        StateIndex is1 = outerBandStructure.bteToState(ibte1);
        StateIndex is2 = outerBandStructure.bteToState(ibte2);
        auto tupp1 = outerBandStructure.getIndex(is1);
        auto tupp2 = outerBandStructure.getIndex(is2);
        WavevectorIndex ik1Index = std::get<0>(tupp1);
        WavevectorIndex ik2Index = std::get<0>(tupp2);
        long ik1Irr = ik1Index.get();
        long ik2Irr = ik2Index.get();
        for ( long ik2 : outerBandStructure.getReduciblesFromIrreducible(ik2Irr) ) {
          std::pair<int, int> xx = std::make_pair(ik1Irr, ik2);
          x.insert(xx);
        }
      }
      std::vector<std::pair<int, int>> wavevectorPair;
      for (auto t : x)
        wavevectorPair.push_back(t);

      // find set of q2
      std::set<long> q2Indexes;
      for (auto tup : x) {
        auto iq2 = std::get<1>(tup);
        q2Indexes.insert(iq2);
      }

      std::vector<std::tuple<std::vector<long>, long>> pairIterator;
      // find all q1 for fixes q2
      for (long iq2 : q2Indexes) {
        std::vector<long> q1Indexes;

        // note: I'm assuming that the pairs in x are unique
        for (auto tup : x) {
          auto iq1 = std::get<0>(tup);
          auto iq2_ = std::get<1>(tup);
          if (iq2 == iq2_) {
            q1Indexes.push_back(iq1);
          }
        }

        auto t = std::make_tuple(q1Indexes, iq2);
        pairIterator.push_back(t);
      }
      return pairIterator;
    }
  }
}

long ScatteringMatrix::getSMatrixIndex(BteIndex &bteIndex,
                                       CartIndex &cartIndex) {
  if (context.getUseSymmetries()) {
    long is = bteIndex.get();
    long alpha = cartIndex.get();
    return compress2Indeces(is, alpha, numStates, 3);
  } else {
    return bteIndex.get();
  }
}

std::tuple<BteIndex, CartIndex>
ScatteringMatrix::getSMatrixIndex(const long &iMat) {
  if (context.getUseSymmetries()) {
    auto t = decompress2Indeces(iMat, numStates, 3);
    return {BteIndex(std::get<0>(t)), CartIndex(std::get<1>(t))};
  } else {
    return {BteIndex(iMat), CartIndex(0)};
  }
}
