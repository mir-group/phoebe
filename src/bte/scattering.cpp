#include "scattering.h"
#include <algorithm>
#include <numeric>  // std::iota
#include <set>
#include <nlohmann/json.hpp>
#include <fstream>
#include <iomanip>
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

  numCalcs = statisticsSweep.getNumCalcs();

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

  std::vector<VectorBTE> emptyVector;

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
  for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
    for (int iDim = 0; iDim < dimensionality_; iDim++) {
      for (long is = 0; is < numStates; is++) {
        outPopulation(iCalc, iDim, is) -=
            internalDiagonal(iCalc, 0, is) * inPopulation(iCalc, iDim, is);
      }
    }
  }
  return outPopulation;
}

std::vector<VectorBTE> ScatteringMatrix::offDiagonalDot(
    std::vector<VectorBTE> &inPopulations) {
  // outPopulation = outPopulation - internalDiagonal * inPopulation;
  std::vector<VectorBTE> outPopulations = dot(inPopulations);
  for (unsigned int iVec = 0; iVec < inPopulations.size(); iVec++) {
    for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
      for (int iDim = 0; iDim < dimensionality_; iDim++) {
        for (long is = 0; is < numStates; is++) {
          outPopulations[iVec](iCalc, iDim, is) -=
              internalDiagonal(iCalc, 0, is) *
              inPopulations[iVec](iCalc, iDim, is);
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
    outPopulation.data.setZero();
    // note: we are assuming that ScatteringMatrix has numCalcs = 1
    for (int idim = 0; idim < inPopulation.dimensionality; idim++) {
      for (auto tup : theMatrix.getAllLocalStates()) {
        auto i1 = std::get<0>(tup);
        auto j1 = std::get<1>(tup);
        outPopulation(0, idim, i1) +=
            theMatrix(i1, j1) * inPopulation(0, idim, j1);
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

std::vector<VectorBTE> ScatteringMatrix::dot(
    std::vector<VectorBTE> &inPopulations) {
  if (highMemory) {
    std::vector<VectorBTE> outPopulations;
    for (auto inPopulation : inPopulations) {
      VectorBTE outPopulation(statisticsSweep, outerBandStructure,
                              inPopulation.dimensionality);
      outPopulation.data.setZero();
      // note: we are assuming that ScatteringMatrix has numCalcs = 1
      for (int idim = 0; idim < inPopulation.dimensionality; idim++) {
        for (auto tup : theMatrix.getAllLocalStates()) {
          auto i1 = std::get<0>(tup);
          auto j1 = std::get<1>(tup);
          outPopulation(0, idim, i1) +=
              theMatrix(i1, j1) * inPopulation(0, idim, j1);
        }
      }
      mpi->allReduceSum(&outPopulation.data);
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
      internalDiagonal(0, 0, ind1) /= term1;
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
      internalDiagonal(iCalc, 0, ind1) *= term1;
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
    if (isMatrixOmega || innerBandStructure.getParticle().isElectron()) {
      for (long iCalc = 0; iCalc < internalDiagonal.numCalcs; iCalc++) {
        for (long is = 0; is < internalDiagonal.numStates; is++) {
          times(iCalc, 0, is) = 1. / times(iCalc, 0, is);
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

          times(iCalc, 0, is) = popTerm / times(iCalc, 0, is);
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
    linewidths.setConst(1./crt);
    linewidths.excludeIndeces = excludeIndeces;
    return linewidths;
  } else {
    VectorBTE linewidths = internalDiagonal;
    if (isMatrixOmega || innerBandStructure.getParticle().isElectron()) {
      linewidths.excludeIndeces = excludeIndeces;
      return linewidths;
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

          linewidths(iCalc, 0, is) /= popTerm;
        }
      }
      linewidths.excludeIndeces = excludeIndeces;
      return linewidths;
    }
  }
}

void ScatteringMatrix::outputToJSON(std::string outFileName) {

  if(!mpi->mpiHead()) return;

  VectorBTE times = getSingleModeTimes();
  VectorBTE tmpLinewidths = getLinewidths();

  std::string energyUnit;
  std::string particleType;
  double energyConversion;
  auto particle = outerBandStructure.getParticle();
  energyConversion = energyRyToEv;
  energyUnit = "eV";
  if(particle.isPhonon()) {
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
    for (int ik = 0; ik < outerBandStructure.getNumPoints(); ik++) {
      auto ikIndex = WavevectorIndex(ik);

      std::vector<double> bandsT;
      std::vector<double> bandsL;
      std::vector<std::vector<double>> bandsV;
      std::vector<double> bandsE;
      // loop over bands here
      // get numBands at this point, in case it's an active bandstructure
      for (int ib = 0; ib < outerBandStructure.getNumBands(ikIndex); ib++) {
        auto ibIndex = BandIndex(ib);
        int is = outerBandStructure.getIndex(ikIndex,ibIndex);
        double ene = outerBandStructure.getEnergy(is);
        auto vel = outerBandStructure.getGroupVelocity(is);
        bandsE.push_back(ene * energyConversion);
        double tau = times(iCalc, 0, is); // only zero dim is meaningful
        bandsT.push_back(tau * timeRyToFs);
        double linewidth = tmpLinewidths(iCalc, 0, is); // only zero dim is meaningful
        bandsL.push_back(linewidth * energyRyToEv);

        std::vector<double> iDimsV;
        // loop over dimensions
        for (int iDim = 0; iDim < dimensionality_; iDim++) {
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
  for (long ik = 0; ik < outerBandStructure.getNumPoints(); ik++) {
    // save the wavevectors
    auto coord = points.getPointCoords(ik);
    meshCoords.push_back({coord[0],coord[1],coord[2]});
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

std::tuple<VectorBTE, ParallelMatrix<double>> ScatteringMatrix::diagonalize() {
  //    std::vector<double> eigenvalues;
  //    ParallelMatrix<double> eigenvectors;
  auto tup = theMatrix.diagonalize();
  auto eigenvalues = std::get<0>(tup);
  auto eigenvectors = std::get<1>(tup);

  // place eigenvalues in an VectorBTE object
  VectorBTE eigvals(statisticsSweep, outerBandStructure, 1);
  for (long is = 0; is < numStates; is++) {
    eigvals(0, 0, is) = eigenvalues[is];
  }
  eigvals.excludeIndeces = excludeIndeces;

  // correct normalization of eigenvectors
  double volume = outerBandStructure.getPoints().getCrystal().getVolumeUnitCell(
      context.getDimensionality());
  eigenvectors *= sqrt(innerBandStructure.getNumPoints(true) * volume);

  return {eigvals, eigenvectors};
}

std::vector<std::tuple<std::vector<long>, long>>
ScatteringMatrix::getIteratorWavevectorPairs(const int &switchCase,
                                             const bool &rowMajor) {
  if ( rowMajor ) {
    if (switchCase != 0) {
      std::vector<std::tuple<std::vector<long>, long>> pairIterator;

      size_t a = outerBandStructure.getNumPoints();
      std::vector<long> outerIterator = mpi->divideWorkIter(a);
      // Note: phScatteringMatrix needs iq2 to be the outer loop
      // in order to be efficient!
      std::vector<long> innerIterator(innerBandStructure.getNumPoints());
      // populate vector with integers from 0 to numPoints-1
      std::iota(std::begin(innerIterator), std::end(innerIterator), 0);

      for (long iq1 : innerIterator) {
        auto t = std::make_tuple(outerIterator, iq1);
        pairIterator.push_back(t);
      }
      return pairIterator;

    } else {

      // list in form [[0,0],[1,0],[2,0],...]
      std::set<std::pair<int, int>> x;
      for (auto tup0 : theMatrix.getAllLocalStates()) {
        auto is1 = std::get<0>(tup0);
        auto is2 = std::get<1>(tup0);
        auto tup1 = outerBandStructure.getIndex(is1);
        auto ik1 = std::get<0>(tup1);
        auto tup2 = outerBandStructure.getIndex(is2);
        auto ik2 = std::get<0>(tup2);
        std::pair<int, int> xx = std::make_pair(ik1.get(), ik2.get());
        x.insert(xx);
      }
      std::vector<std::pair<int, int>> wavevectorPair;
      for (auto t : x) {
        wavevectorPair.push_back(t);
      }

      // find set of q1
      std::set<long> q1Indexes;
      for (auto tup : x) {
        auto iq1 = std::get<1>(tup);
        q1Indexes.insert(iq1);
      }

      std::vector<std::tuple<std::vector<long>, long>> pairIterator;
      // find all q2 for fixed q1
      for (long iq1 : q1Indexes) {
        std::vector<long> q2Indexes;

        // note: I'm assuming that the pairs in x are unique
        for (auto tup : x) {
          auto iq2 = std::get<0>(tup);
          auto iq1_ = std::get<1>(tup);
          if (iq1 == iq1_) {
            q2Indexes.push_back(iq2);
          }
        }

        auto t = std::make_tuple(q2Indexes, iq1);
        pairIterator.push_back(t);
      }
      return pairIterator;
    }

  } else {

    if (switchCase != 0) {
      size_t a = innerBandStructure.getNumPoints();
      std::vector<long> wavevectorIterator = mpi->divideWorkIter(a);
      // Note: phScatteringMatrix needs iq2 to be the outer loop
      // in order to be efficient!

      std::vector<long> outerIterator(outerBandStructure.getNumPoints());
      // populate vector with integers from 0 to numPoints-1
      std::iota(std::begin(outerIterator), std::end(outerIterator), 0);

      std::vector<std::tuple<std::vector<long>, long>> pairIterator(a);
      int i=0;
      for (long iq2 : wavevectorIterator) {
        auto t = std::make_tuple(outerIterator, iq2);
        pairIterator[i] = t;
        i++;
      }

      return pairIterator;

    } else {

      // list in form [[0,0],[1,0],[2,0],...]
      std::set<std::pair<int, int>> x;
      for (auto tup0 : theMatrix.getAllLocalStates()) {
        auto is1 = std::get<0>(tup0);
        auto is2 = std::get<1>(tup0);
        auto tup1 = outerBandStructure.getIndex(is1);
        auto ik1 = std::get<0>(tup1);
        auto tup2 = outerBandStructure.getIndex(is2);
        auto ik2 = std::get<0>(tup2);
        std::pair<int, int> xx = std::make_pair(ik1.get(), ik2.get());
        x.insert(xx);
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
