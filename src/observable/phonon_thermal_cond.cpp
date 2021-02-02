#include "phonon_thermal_cond.h"

#include "constants.h"
#include "mpiHelper.h"
#include <ctime>
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>

PhononThermalConductivity::PhononThermalConductivity(
    Context &context_, StatisticsSweep &statisticsSweep_, Crystal &crystal_,
    BaseBandStructure &bandStructure_)
    : Observable(context_, statisticsSweep_, crystal_),
      bandStructure(bandStructure_) {
  tensordxd =
      Eigen::Tensor<double, 3>(numCalculations, dimensionality, dimensionality);
  tensordxd.setZero();
}

// copy constructor
PhononThermalConductivity::PhononThermalConductivity(
    const PhononThermalConductivity &that)
    : Observable(that), bandStructure(that.bandStructure) {}

// copy assignment
PhononThermalConductivity &
PhononThermalConductivity::operator=(const PhononThermalConductivity &that) {
  Observable::operator=(that);
  if (this != &that) {
    bandStructure = that.bandStructure;
  }
  return *this;
}

PhononThermalConductivity
PhononThermalConductivity::operator-(const PhononThermalConductivity &that) {
  PhononThermalConductivity newObservable(context, statisticsSweep, crystal,
                                          bandStructure);
  baseOperatorMinus(newObservable, that);
  return newObservable;
}

void PhononThermalConductivity::calcFromCanonicalPopulation(VectorBTE &f) {
  VectorBTE n = f;
  n.canonical2Population(); // n = bose (bose+1) f
  calcFromPopulation(n);
}

void PhononThermalConductivity::calcFromPopulation(VectorBTE &n) {
  double norm = 1. / context.getQMesh().prod() /
                crystal.getVolumeUnitCell(dimensionality);

  auto excludeIndices = n.excludeIndices;

  tensordxd.setZero();

#pragma omp parallel default(none) shared(excludeIndices, n, norm)
  {
    // we do manually the reduction, to avoid custom type declaration
    // which is not always allowed by the compiler e.g. by clang

    // first omp parallel for on a private variable
    Eigen::Tensor<double, 3> tensorPrivate(numCalculations, dimensionality,
                                           dimensionality);
    tensorPrivate.setZero();

#pragma omp for nowait
    for (int is : bandStructure.parallelIrrStateIterator()) {
      StateIndex isIdx(is);
      double en = bandStructure.getEnergy(isIdx);
      Eigen::Vector3d velIrr = bandStructure.getGroupVelocity(isIdx);

      int iBte = bandStructure.stateToBte(isIdx).get();

      // skip the acoustic phonons
      if (std::find(excludeIndices.begin(), excludeIndices.end(), iBte) !=
          excludeIndices.end())
        continue;

      auto rots = bandStructure.getRotationsStar(isIdx);
      for (const Eigen::Matrix3d &rot : rots) {
        auto vel = rot * velIrr;

        for (int iCalc = 0; iCalc < statisticsSweep.getNumCalculations();
             iCalc++) {

          Eigen::Vector3d nRot;
          for (int i : {0, 1, 2}) {
            nRot(i) = n(iCalc, i, iBte);
          }
          nRot = rot * nRot;

          for (int j : {0, 1, 2}) {
            for (int i : {0, 1, 2}) {
              tensorPrivate(iCalc, i, j) += nRot(i) * vel(j) * en * norm;
            }
          }
        }
      }
    }

// now we do the reduction thread by thread
#pragma omp critical
    {
      for (int j : {0, 1, 2}) {
        for (int i : {0, 1, 2}) {
          for (int iCalc = 0; iCalc < statisticsSweep.getNumCalculations();
               iCalc++) {
            tensordxd(iCalc, i, j) += tensorPrivate(iCalc, i, j);
          }
        }
      }
    }
  }
  // lastly, the states were distributed with MPI
  mpi->allReduceSum(&tensordxd);
}

void PhononThermalConductivity::calcVariational(VectorBTE &af, VectorBTE &f,
                                                VectorBTE &scalingCG) {
  VectorBTE fUnscaled = f / scalingCG;
  calcFromCanonicalPopulation(fUnscaled);

  double norm = 1. / context.getQMesh().prod() /
                crystal.getVolumeUnitCell(dimensionality);
  auto excludeIndices = f.excludeIndices;

  int numCalculations = statisticsSweep.getNumCalculations();

  tensordxd *= tensordxd.constant(2.);

  Eigen::Tensor<double, 3> tmpTensor = tensordxd.constant(0.);
#pragma omp parallel default(none) shared(                                     \
    excludeIndices, bandStructure, tmpTensor, norm, af, f, numCalculations)
  {
    Eigen::Tensor<double, 3> tmpTensorPrivate(numCalculations, 3, 3);
    tmpTensorPrivate.setConstant(0.);
#pragma omp for nowait
    for (int is : bandStructure.parallelIrrStateIterator()) {

      // skip the acoustic phonons
      if (std::find(excludeIndices.begin(), excludeIndices.end(), is) !=
          excludeIndices.end()) {
        continue;
      }

      auto isIndex = StateIndex(is);
      BteIndex iBteIndex = bandStructure.stateToBte(isIndex);
      int isBte = iBteIndex.get();
      auto rots = bandStructure.getRotationsStar(isIndex);

      for (const Eigen::Matrix3d &rot : rots) {

        for (int iCalc = 0; iCalc < numCalculations; iCalc++) {

          auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
          double temp = calcStat.temperature;
          double norm2 = norm * temp * temp;

          Eigen::Vector3d fRot, afRot;
          for (int i : {0, 1, 2}) {
            fRot(i) = f(iCalc, i, isBte);
            afRot(i) = af(iCalc, i, isBte);
          }
          fRot = rot * fRot;
          afRot = rot * afRot;

          for (int i : {0, 1, 2}) {
            for (int j : {0, 1, 2}) {
              tmpTensorPrivate(iCalc, i, j) += fRot(i) * afRot(j) * norm2;
            }
          }
        }
      }
    }
#pragma omp critical
    for (int j = 0; j < dimensionality; j++) {
      for (int i = 0; i < dimensionality; i++) {
        for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
          tmpTensor(iCalc, i, j) += tmpTensorPrivate(iCalc, i, j);
        }
      }
    }
  }
  mpi->allReduceSum(&tmpTensor);
  tensordxd -= tmpTensor;
}

void PhononThermalConductivity::calcFromRelaxons(
    Context &context, StatisticsSweep &statisticsSweep,
    ParallelMatrix<double> &eigenvectors, PhScatteringMatrix &scatteringMatrix,
    const Eigen::VectorXd &eigenvalues) {

  int dimensionality = context.getDimensionality();
  auto particle = bandStructure.getParticle();

  int iCalc = 0;
  double temp = statisticsSweep.getCalcStatistics(0).temperature;
  double chemPot = statisticsSweep.getCalcStatistics(0).chemicalPotential;

  VectorBTE population(statisticsSweep, bandStructure, dimensionality);

  if (context.getUseSymmetries()) {

    Eigen::VectorXd relPopulation(eigenvalues.size());
    relPopulation.setZero();

#pragma omp parallel default(none)                                             \
    shared(relPopulation, temp, chemPot, scatteringMatrix, eigenvectors,       \
           eigenvalues, particle)
    {
      Eigen::VectorXd relPopPrivate(eigenvalues.size());
      relPopPrivate.setZero();
#pragma omp for nowait
      for (auto tup0 : eigenvectors.getAllLocalStates()) {
        int iMat1 = std::get<0>(tup0);
        int alpha = std::get<1>(tup0);
        auto tup1 = scatteringMatrix.getSMatrixIndex(iMat1);
        BteIndex iBteIndex = std::get<0>(tup1);
        CartIndex dimIndex = std::get<1>(tup1);
        int iDim = dimIndex.get();
        StateIndex isIndex = bandStructure.bteToState(iBteIndex);
        //
        auto vel = bandStructure.getGroupVelocity(isIndex);
        double en = bandStructure.getEnergy(isIndex);
        double term = sqrt(particle.getPopPopPm1(en, temp, chemPot));
        double dndt = particle.getDndt(en, temp, chemPot);
        if (eigenvalues(alpha) > 0.) {
          relPopPrivate(alpha) += dndt / term * vel(iDim) *
                                  eigenvectors(iMat1, alpha) /
                                  eigenvalues(alpha);
        }
      }
#pragma omp critical
      for (int alpha = 0; alpha < eigenvalues.size(); alpha++) {
        relPopulation(alpha) += relPopPrivate(alpha);
      }
    }
    mpi->allReduceSum(&relPopulation);

    // back rotate to phonon coordinates

#pragma omp parallel default(none)                                             \
    shared(bandStructure, eigenvectors, relPopulation, population,             \
           scatteringMatrix, iCalc)
    {
      Eigen::MatrixXd popPrivate(3, bandStructure.irrStateIterator().size());
      popPrivate.setZero();

#pragma omp for nowait
      for (auto tup0 : eigenvectors.getAllLocalStates()) {
        int iMat1 = std::get<0>(tup0);
        int alpha = std::get<1>(tup0);
        auto tup1 = scatteringMatrix.getSMatrixIndex(iMat1);
        BteIndex iBteIndex = std::get<0>(tup1);
        CartIndex dimIndex = std::get<1>(tup1);
        int iBte = iBteIndex.get();
        int iDim = dimIndex.get();
        popPrivate(iDim, iBte) +=
            eigenvectors(iMat1, alpha) * relPopulation(alpha);
      }

#pragma omp critical
      for (int is = 0; is < popPrivate.cols(); is++) {
        for (int iDim : {0, 1, 2}) {
          population(iCalc, iDim, is) += popPrivate(iDim, is);
        }
      }
    }
    mpi->allReduceSum(&population.data);

  } else { // case without symmetries ------------------------------------------

    Eigen::MatrixXd relPopulation(eigenvalues.size(), 3);
    relPopulation.setZero();
#pragma omp parallel default(none)                                             \
    shared(eigenvectors, eigenvalues, temp, chemPot, particle, relPopulation)
    {
      Eigen::MatrixXd relPopPrivate(eigenvalues.size(), 3);
      relPopPrivate.setZero();
#pragma omp for nowait
      for (auto tup0 : eigenvectors.getAllLocalStates()) {
        int is = std::get<0>(tup0);
        StateIndex isIdx(is);
        int alpha = std::get<1>(tup0);
        //
        double en = bandStructure.getEnergy(isIdx);
        if (eigenvalues(alpha) > 0. && en >= 0.) {
          auto vel = bandStructure.getGroupVelocity(isIdx);
          double term = sqrt(particle.getPopPopPm1(en, temp, chemPot));
          double dndt = particle.getDndt(en, temp, chemPot);
          for (int i : {0, 1, 2}) {
            relPopPrivate(alpha, i) += dndt / term * vel(i) *
                                       eigenvectors(is, alpha) /
                                       eigenvalues(alpha);
          }
        }
      }
#pragma omp critical
      for (int alpha = 0; alpha < eigenvalues.size(); alpha++) {
        for (int i : {0, 1, 2}) {
          relPopulation(alpha, i) += relPopPrivate(alpha, i);
        }
      }
    }
    mpi->allReduceSum(&relPopulation);
    // back rotate to phonon coordinates
#pragma omp parallel default(none)                                             \
    shared(eigenvectors, relPopulation, population, iCalc)
    {
      Eigen::MatrixXd popPrivate(3, bandStructure.getNumStates());
      popPrivate.setZero();
#pragma omp for nowait
      for (auto tup0 : eigenvectors.getAllLocalStates()) {
        int is = std::get<0>(tup0);
        int alpha = std::get<1>(tup0);
        for (int i : {0, 1, 2}) {
          popPrivate(i, is) +=
              eigenvectors(is, alpha) * relPopulation(alpha, i);
        }
      }
#pragma omp critical
      for (int is = 0; is < bandStructure.getNumStates(); is++) {
        for (int i : {0, 1, 2}) {
          population(iCalc, i, is) += popPrivate(i, is);
        }
      }
    }
    mpi->allReduceSum(&population.data);
  }
  // put back the rescaling factor

#pragma omp parallel for default(none)                                         \
    shared(population, particle, temp, chemPot, iCalc, bandStructure)
  for (int is : bandStructure.irrStateIterator()) {
    auto isIdx = StateIndex(is);
    double en = bandStructure.getEnergy(isIdx);
    int iBte = bandStructure.stateToBte(isIdx).get();
    if (en > 0.) {
      double term = particle.getPopPopPm1(en, temp, chemPot);
      for (int iDim : {0, 1, 2}) {
        population(iCalc, iDim, iBte) *= sqrt(term);
      }
    }
  }

  calcFromPopulation(population);
}

void PhononThermalConductivity::print() {
  if (!mpi->mpiHead())
    return;

  std::string units;
  if (dimensionality == 1) {
    units = "W m / K";
  } else if (dimensionality == 2) {
    units = "W / K";
  } else {
    units = "W / m / K";
  }

  std::cout << "\n";
  std::cout << "Thermal Conductivity (" << units << ")\n";

  for (int iCalc = 0; iCalc < statisticsSweep.getNumCalculations(); iCalc++) {
    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStat.temperature;

    std::cout << std::fixed;
    std::cout.precision(2);
    std::cout << "Temperature: " << temp * temperatureAuToSi << " (K)\n";
    std::cout.precision(5);
    for (int i = 0; i < dimensionality; i++) {
      std::cout << "  " << std::scientific;
      for (int j = 0; j < dimensionality; j++) {
        std::cout << " " << std::setw(13) << std::right;
        std::cout << tensordxd(iCalc, i, j) * thConductivityAuToSi;
      }
      std::cout << "\n";
    }
    std::cout << std::endl;
  }
}

void PhononThermalConductivity::outputToJSON(const std::string &outFileName) {
  if (!mpi->mpiHead())
    return;

  std::string units;
  if (dimensionality == 1) {
    units = "W m / K";
  } else if (dimensionality == 2) {
    units = "W / K";
  } else {
    units = "W / m / K";
  }

  std::vector<double> temps;
  std::vector<std::vector<std::vector<double>>> conductivities;
  for (int iCalc = 0; iCalc < statisticsSweep.getNumCalculations(); iCalc++) {

    // store temperatures
    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStat.temperature;
    temps.push_back(temp * temperatureAuToSi);

    // store conductivity
    std::vector<std::vector<double>> rows;
    for (int i = 0; i < dimensionality; i++) {
      std::vector<double> cols;
      for (int j = 0; j < dimensionality; j++) {
        cols.push_back(tensordxd(iCalc, i, j) * thConductivityAuToSi);
      }
      rows.push_back(cols);
    }
    conductivities.push_back(rows);
  }

  // output to json
  nlohmann::json output;
  output["temperatures"] = temps;
  output["thermalConductivity"] = conductivities;
  output["temperatureUnit"] = "K";
  output["thermalConductivityUnit"] = units;
  output["particleType"] = "phonon";
  std::ofstream o(outFileName);
  o << std::setw(3) << output << std::endl;
  o.close();
}

void PhononThermalConductivity::print(const int &iter) {
  if (!mpi->mpiHead())
    return;

  // get the time
  time_t currentTime;
  currentTime = time(nullptr);
  // and format the time nicely
  char s[200];
  struct tm *p = localtime(&currentTime);
  strftime(s, 200, "%F, %T", p);

  std::cout << "Iteration: " << iter << " | " << s << "\n";
  for (int iCalc = 0; iCalc < statisticsSweep.getNumCalculations(); iCalc++) {
    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStat.temperature;
    std::cout << std::fixed;
    std::cout.precision(2);
    std::cout << "T = " << temp * temperatureAuToSi << ", k = ";
    std::cout.precision(5);
    for (int i = 0; i < dimensionality; i++) {
      std::cout << std::scientific;
      std::cout << tensordxd(iCalc, i, i) * thConductivityAuToSi << " ";
    }
    std::cout << "\n";
  }
  std::cout << std::endl;
}

int PhononThermalConductivity::whichType() { return is2Tensor; }
