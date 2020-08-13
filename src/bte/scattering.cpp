#include "scattering.h"
#include <algorithm>
#include <numeric>  // std::iota
#include <set>
#include "constants.h"
#include "mpiHelper.h"

ScatteringMatrix::ScatteringMatrix(Context &context_,
                                   StatisticsSweep &statisticsSweep_,
                                   BaseBandStructure &innerBandStructure_,
                                   BaseBandStructure &outerBandStructure_)
    : context(context_),
      statisticsSweep(statisticsSweep_),
      innerBandStructure(innerBandStructure_),
      outerBandStructure(outerBandStructure_),
      internalDiagonal(statisticsSweep, outerBandStructure, 1) {
  numStates = outerBandStructure.getNumStates();
  numPoints = outerBandStructure.getNumPoints();

  double constantRelaxationTime = context.getConstantRelaxationTime();
  if (constantRelaxationTime > 0.) {
    constantRTA = true;
    return;
  }

  dimensionality_ = context.getDimensionality();

  highMemory = context.getScatteringMatrixInMemory();

  smearing = DeltaFunction::smearingFactory(context, innerBandStructure);

  // the difference arises from the fact that vectors in the BTE have
  // numCalcs set to nTemp * 3, whereas the matrix only depends on nTemp.
  // so, when we compute the action of the scattering matrix, we must take
  // into account for this misalignment
  if (highMemory) {
    numCalcs = statisticsSweep.getNumCalcs();
  } else {
    numCalcs = statisticsSweep.getNumCalcs() * dimensionality_;
  }

  // we want to know the state index of acoustic modes at gamma,
  // so that we can set their populations to zero
  if (outerBandStructure.getParticle().isPhonon()) {
    for (long is = 0; is < numStates; is++) {
      double en = outerBandStructure.getEnergy(is);
      if (en < 0.1 / ryToCmm1) {  // cutoff at 0.1 cm^-1
        excludeIndeces.push_back(is);
      }
    }
  }
}

ScatteringMatrix::~ScatteringMatrix() { delete smearing; }

// copy constructor
ScatteringMatrix::ScatteringMatrix(const ScatteringMatrix &that)
    : context(that.context),
      statisticsSweep(that.statisticsSweep),
      smearing(that.smearing),
      innerBandStructure(that.innerBandStructure),
      outerBandStructure(that.outerBandStructure),
      constantRTA(that.constantRTA),
      highMemory(that.highMemory),
      internalDiagonal(that.internalDiagonal),
      theMatrix(that.theMatrix),
      numStates(that.numStates),
      numPoints(that.numPoints),
      numCalcs(that.numCalcs),
      dimensionality_(that.dimensionality_),
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

  if (constantRTA) return;  // nothing to construct

  if (highMemory) {
    if (numCalcs > 1) {
      // note: one could write code around this
      // but the methods are very memory intensive for production runs
      Error e(
          "High memory BTE methods can only work with one "
          "temperature and/or chemical potential in a single run");
    }
    theMatrix = ParallelMatrix<double>(numStates, numStates);

    // calc matrix and linew.
    builder(theMatrix, &internalDiagonal, nullptr, nullptr);
  } else {
    // calc linewidths only
    builder(theMatrix, &internalDiagonal, nullptr, nullptr);
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

// Matrix<double> ScatteringMatrix::dot(const Matrix<double> &otherMatrix) {
//    return theMatrix * otherMatrix;
//}

VectorBTE ScatteringMatrix::offDiagonalDot(VectorBTE &inPopulation) {
  VectorBTE outPopulation = dot(inPopulation);
  // outPopulation = outPopulation - internalDiagonal * inPopulation;
  for (long i = 0; i < outPopulation.numCalcs; i++) {
    auto tup = inPopulation.loc2Glob(i);
    auto imu = std::get<0>(tup);
    auto it = std::get<1>(tup);
    auto idim = std::get<2>(tup);
    auto j = internalDiagonal.glob2Loc(imu, it, DimIndex(0));
    for (long is = 0; is < numStates; is++) {
      outPopulation.data(i, is) -=
          internalDiagonal.data(j, is) * inPopulation.data(i, is);
    }
  }
  return outPopulation;
}

VectorBTE ScatteringMatrix::dot(VectorBTE &inPopulation) {
  if (highMemory) {
    VectorBTE outPopulation(statisticsSweep, outerBandStructure,
                            inPopulation.dimensionality);
    outPopulation.data.setZero();
    // note: we are assuming that ScatteringMatrix has numCalcs = 1
    for (int idim = 0; idim < inPopulation.dimensionality; idim++) {
      for (auto tup : theMatrix.getAllLocalStates()) {
        auto i1 = std::get<0>(tup);
        auto j1 = std::get<1>(tup);
        outPopulation.data(idim, i1) +=
            theMatrix(i1, j1) * inPopulation.data(idim, j1);
      }
    }
    mpi->allReduceSum(&outPopulation.data);
    return outPopulation;
  } else {
    VectorBTE outPopulation(statisticsSweep, outerBandStructure,
                            inPopulation.dimensionality);
    builder(theMatrix, nullptr, &inPopulation, &outPopulation);
    return outPopulation;
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

  if (isMatrixOmega) {  // it's already with the scaling of omega
    return;
  }

  long iCalc = 0;  // as there can only be one temperature

  auto particle = outerBandStructure.getParticle();
  auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
  double temp = calcStatistics.temperature;
  double chemPot = calcStatistics.chemicalPotential;

  for (auto tup : theMatrix.getAllLocalStates()) {
    auto ind1 = std::get<0>(tup);
    auto ind2 = std::get<1>(tup);
    if (std::find(excludeIndeces.begin(), excludeIndeces.end(), ind1) !=
        excludeIndeces.end())
      continue;
    if (std::find(excludeIndeces.begin(), excludeIndeces.end(), ind2) !=
        excludeIndeces.end())
      continue;

    double en1 = outerBandStructure.getEnergy(ind1);
    double en2 = outerBandStructure.getEnergy(ind2);

    // n(n+1) for bosons, n(1-n) for fermions
    double term1 = particle.getPopPopPm1(en1, temp, chemPot);
    double term2 = particle.getPopPopPm1(en2, temp, chemPot);

    if (ind1 == ind2) {
      internalDiagonal.data(iCalc, ind1) /= term1;
    }

    theMatrix(ind1, ind2) /= sqrt(term1 * term2);
  }
  isMatrixOmega = true;
}

// add a flag to remember if we have A or Omega
void ScatteringMatrix::omega2A() {
  if (!highMemory) {
    Error e("a2Omega only works if the matrix is stored in memory");
  }

  if (theMatrix.rows() == 0) {
    Error e("The scattering matrix hasn't been built yet");
  }

  if (!isMatrixOmega) {  // it's already with the scaling of A
    return;
  }

  long iCalc = 0;  // as there can only be one temperature

  auto particle = outerBandStructure.getParticle();
  auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
  double temp = calcStatistics.temperature;
  double chemPot = calcStatistics.chemicalPotential;

  for (auto tup : theMatrix.getAllLocalStates()) {
    auto ind1 = std::get<0>(tup);
    auto ind2 = std::get<1>(tup);
    if (std::find(excludeIndeces.begin(), excludeIndeces.end(), ind1) !=
        excludeIndeces.end())
      continue;
    if (std::find(excludeIndeces.begin(), excludeIndeces.end(), ind2) !=
        excludeIndeces.end())
      continue;

    double en1 = outerBandStructure.getEnergy(ind1);
    double en2 = outerBandStructure.getEnergy(ind2);

    // n(n+1) for bosons, n(1-n) for fermions
    double term1 = particle.getPopPopPm1(en1, temp, chemPot);
    double term2 = particle.getPopPopPm1(en2, temp, chemPot);

    if (ind1 == ind2) {
      internalDiagonal.data(iCalc, ind1) *= term1;
    }

    theMatrix(ind1, ind2) *= sqrt(term1 * term2);
  }
  isMatrixOmega = false;
}

// to compute the RTA, get the single mode relaxation times
VectorBTE ScatteringMatrix::getSingleModeTimes() {
  if (constantRTA) {
    double crt = context.getConstantRelaxationTime();
    VectorBTE times(statisticsSweep, outerBandStructure, 1);
    times.setConst(crt);
    times.excludeIndeces = excludeIndeces;
    return times;
  } else {
    VectorBTE times = internalDiagonal;
    if (isMatrixOmega) {
      for (long iCalc = 0; iCalc < internalDiagonal.numCalcs; iCalc++) {
        for (long is = 0; is < internalDiagonal.numStates; is++) {
          times.data(iCalc, is) = 1. / times.data(iCalc, is);
        }
      }
      times.excludeIndeces = excludeIndeces;
      return times;
    } else {  // A_nu,nu = N(1+-N) / tau
      auto particle = outerBandStructure.getParticle();
      for (long iCalc = 0; iCalc < internalDiagonal.numCalcs; iCalc++) {
        auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
        double temp = calcStatistics.temperature;
        double chemPot = calcStatistics.chemicalPotential;

        for (long is = 0; is < internalDiagonal.numStates; is++) {
          double en = outerBandStructure.getEnergy(is);

          // n(n+1) for bosons, n(1-n) for fermions
          double popTerm = particle.getPopPopPm1(en, temp, chemPot);

          times.data(iCalc, is) = popTerm / times.data(iCalc, is);
        }
      }
      times.excludeIndeces = excludeIndeces;
      return times;
    }
  }
}

std::tuple<VectorBTE, ParallelMatrix<double>> ScatteringMatrix::diagonalize() {
  //    std::vector<double> eigenvalues;
  //    ParallelMatrix<double> eigenvectors;
  auto tup = theMatrix.diagonalize();
  auto eigenvalues = std::get<0>(tup);
  auto eigenvectors = std::get<1>(tup);

  // place eigenvalues in an VectorBTE object
  VectorBTE eigvals(statisticsSweep, outerBandStructure, 1);
  for (long is = 0; is < numStates; is++) {
    eigvals.data(0, is) = eigenvalues[is];
  }
  eigvals.excludeIndeces = excludeIndeces;

  // correct normalization of eigenvectors
  double volume = outerBandStructure.getPoints().getCrystal().getVolumeUnitCell(
      context.getDimensionality());
  eigenvectors *= sqrt(innerBandStructure.getNumPoints(true) * volume);

  return {eigvals, eigenvectors};
}

std::vector<std::tuple<std::vector<long>, long>>
ScatteringMatrix::getIteratorWavevectorPairs(const int &switchCase) {
  if (switchCase != 0) {
    std::vector<std::tuple<std::vector<long>, long>> pairIterator;

    size_t a = innerBandStructure.getNumPoints();
    std::vector<int> wavevectorIterator = mpi->divideWorkIter(a);
    // Note: phScatteringMatrix needs iq2 to be the outer loop
    // in order to be efficient!

    std::vector<long> outerIterator(outerBandStructure.getNumPoints());
    // populate vector with integers from 0 to numPoints-1
    std::iota(std::begin(outerIterator), std::end(outerIterator), 0);

    for (long iq2 : wavevectorIterator) {
      auto t = std::make_tuple(outerIterator, iq2);
      pairIterator.push_back(t);
    }
    return pairIterator;

  } else {

    // list in form [[0,0],[1,0],[2,0],...]
    std::set<std::pair<int,int>> x;
    for (auto tup0 : theMatrix.getAllLocalStates()) { 
       auto is1 = std::get<0>(tup0);
       auto is2 = std::get<1>(tup0);
       auto tup1 = outerBandStructure.getIndex(is1);
       auto ik1 = std::get<0>(tup1); 
       auto tup2 = outerBandStructure.getIndex(is2);
       auto ik2 = std::get<0>(tup2); 
       std::pair<int,int> xx = std::make_pair(ik1.get(), ik2.get());
       x.insert(xx);
    }
    std::vector<std::pair<int, int>> wavevectorPair;
    for ( auto t : x ) wavevectorPair.push_back(t);

    // find set of q2
    std::set<long> q2Indexes;
    for ( auto tup : x ) {
      auto iq1 = std::get<0>(tup);
      auto iq2 = std::get<1>(tup);
      q2Indexes.insert(iq2);
    }

    std::vector<std::tuple<std::vector<long>, long>> pairIterator;
    // find all q1 for fixes q2
    for ( long iq2 : q2Indexes ) {
      std::vector<long> q1Indexes;

      // note: I'm assuming that the pairs in x are unique
      for ( auto tup : x ) {
        auto iq1 = std::get<0>(tup);
        auto iq2_ = std::get<1>(tup);
        if ( iq2 == iq2_ ) {
          q1Indexes.push_back(iq1);
        }
      }

      auto t = std::make_tuple(q1Indexes, iq2);
      pairIterator.push_back(t);
    }
    return pairIterator;
  }
}
