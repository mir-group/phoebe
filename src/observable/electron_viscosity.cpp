#include "electron_viscosity.h"
#include "constants.h"
#include "mpiHelper.h"
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>

ElectronViscosity::ElectronViscosity(StatisticsSweep &statisticsSweep_,
                                     Crystal &crystal_,
                                     BaseBandStructure &bandStructure_)
    : Observable(statisticsSweep_, crystal_), bandStructure(bandStructure_) {

  tensordxdxdxd = Eigen::Tensor<double, 5>(
      numCalcs, dimensionality, dimensionality, dimensionality, dimensionality);
  tensordxdxdxd.setZero();
};

// copy constructor
ElectronViscosity::ElectronViscosity(const ElectronViscosity &that)
    : Observable(that), bandStructure(that.bandStructure) {}

// copy assigmnent
ElectronViscosity &ElectronViscosity::operator=(const ElectronViscosity &that) {
  Observable::operator=(that);
  if (this != &that) {
    bandStructure = that.bandStructure;
  }
  return *this;
}

void ElectronViscosity::calcRTA(VectorBTE &tau) {
  double norm = 1. / bandStructure.getNumPoints(true) /
                crystal.getVolumeUnitCell(dimensionality);

  auto particle = bandStructure.getParticle();
  tensordxdxdxd.setZero();

  auto excludeIndeces = tau.excludeIndeces;

#pragma omp parallel
  {
    Eigen::Tensor<double, 5> tmpTensor = tensordxdxdxd.constant(0.);

#pragma omp for nowait
    for (long is : bandStructure.parallelStateIterator()) {
      auto en = bandStructure.getEnergy(is);

      auto vel = bandStructure.getGroupVelocity(is);
      auto q = bandStructure.getWavevector(is);

      // should I skip some states?
      //      if (std::find(excludeIndeces.begin(), excludeIndeces.end(), is)
      //              != excludeIndeces.end())
      //          continue;

      for (long iCalc = 0; iCalc < numCalcs; iCalc++) {

        auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
        double temperature = calcStat.temperature;
        double chemPot = calcStat.chemicalPotential;
        double bosep1 = particle.getPopPopPm1(en, temperature, chemPot);

        for (long i = 0; i < dimensionality; i++) {
          for (long j = 0; j < dimensionality; j++) {
            for (long k = 0; k < dimensionality; k++) {
              for (long l = 0; l < dimensionality; l++) {
                tmpTensor(iCalc, i, j, k, l) += q(i) * vel(j) * q(k) * vel(l) *
                                                bosep1 * tau(iCalc, 0, is) /
                                                temperature * norm;
              }
            }
          }
        }
      }
    }
#pragma omp critical
    for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
      for (long i = 0; i < dimensionality; i++) {
        for (long j = 0; j < dimensionality; j++) {
          for (long k = 0; k < dimensionality; k++) {
            for (long l = 0; l < dimensionality; l++) {
              tensordxdxdxd(iCalc, i, j, k, l) += tmpTensor(iCalc, i, j, k, l);
            }
          }
        }
      }
    }
  }
  mpi->allReduceSum(&tensordxdxdxd);
}

void ElectronViscosity::calcFromRelaxons(Vector0 &vector0, VectorBTE &relTimes,
                                         ElScatteringMatrix &sMatrix,
                                         ParallelMatrix<double> &eigenvectors) {
  Error e("Viscosity from relaxons untested");
  if (numCalcs > 1) {
    Error e("Viscosity for relaxons only for 1 temperature");
  }

  // we decide to skip relaxon states
  // 1) there is a relaxon with zero (or epsilon) eigenvalue -> infinite tau
  // 2) there might be other states with infinite lifetimes, we skip them
  long firstState = 1;
  firstState += relTimes.excludeIndeces.size();

  double volume = crystal.getVolumeUnitCell(dimensionality);
  double numPoints = double(bandStructure.getNumPoints(true));
  long numStates = bandStructure.getNumStates();
  auto particle = bandStructure.getParticle();

  // to simplify, here I do everything considering there is a single
  // temperature (due to memory constraints)

  Eigen::VectorXd A(dimensionality);
  A.setZero();

  long iCalc = 0;
  auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
  double temp = calcStat.temperature;
  double chemPot = calcStat.chemicalPotential;

  for (long is = firstState; is < numStates; is++) {
    auto en = bandStructure.getEnergy(is);
    double bosep1 = particle.getPopPopPm1(en, temp, chemPot); // = n(n+1)
    auto q = bandStructure.getWavevector(is);
    for (long idim = 0; idim < dimensionality; idim++) {
      A(idim) += bosep1 * q(idim) * q(idim);
    }
  }
  A /= temp * numPoints * volume;

  VectorBTE driftEigenvector(statisticsSweep, bandStructure, 3);
  for (long is = firstState; is < numStates; is++) {
    auto en = bandStructure.getEnergy(is);
    double bosep1 = particle.getPopPopPm1(en, temp, chemPot); // = n(n+1)
    auto q = bandStructure.getWavevector(is);
    for (auto idim : {0, 1, 2}) {
      driftEigenvector(0, idim, is) = q(idim) * sqrt(bosep1 / temp / A(idim));
    }
  }

  Eigen::MatrixXd D(3, 3);
  D.setZero();
  //    D = driftEigenvector * sMatrix.dot(driftEigenvector.transpose());
  VectorBTE tmp = sMatrix.dot(driftEigenvector);
  for (long is = firstState; is < numStates; is++) {
    for (auto i : {0, 1, 2}) {
      for (auto j : {0, 1, 2}) {
        D(i, j) += driftEigenvector(0, i, is) * tmp(0, j, is);
      }
    }
  }
  D /= volume * numPoints;

  Eigen::Tensor<double, 3> tmpDriftEigvecs(3, 3, numStates);
  tmpDriftEigvecs.setZero();
  Eigen::MatrixXd W(3, 3);
  W.setZero();
  for (long is = firstState; is < bandStructure.getNumStates(); is++) {
    auto v = bandStructure.getGroupVelocity(is);
    for (auto i : {0, 1, 2}) {
      for (auto j : {0, 1, 2}) {
        tmpDriftEigvecs(i, j, is) = driftEigenvector(0, j, is) * v(i);
        W(i, j) += vector0(0, 0, is) * v(i) * driftEigenvector(0, j, is);
      }
    }
  }
  W /= volume * numPoints;

  Eigen::Tensor<double, 3> w(3, 3, numStates);
  w.setZero();
  for (auto i : {0, 1, 2}) {
    for (auto j : {0, 1, 2}) {
      Eigen::VectorXd x(numStates);
      for (long is = 0; is < numStates; is++) {
        x(is) = tmpDriftEigvecs(i, j, is);
      }

      // auto x2 = x.transpose() * eigenvectors;

      std::vector<double> x2(numStates, 0.);
      for (auto tup : eigenvectors.getAllLocalStates()) {
        auto is1 = std::get<0>(tup);
        auto is2 = std::get<1>(tup);
        x2[is2] += x(is1) * eigenvectors(is1, is2);
      }
      mpi->allReduceSum(&x2);

      for (long is = 0; is < numStates; is++) {
        w(i, j, is) = x2[is] / volume / numPoints;
      }
    }
  }

  // Eq. 9, Simoncelli PRX (2019)
  tensordxdxdxd.setZero();
#pragma omp parallel
  {
    Eigen::Tensor<double, 5> tmpTensor = tensordxdxdxd.constant(0.);
#pragma omp for nowait
    for (long is = firstState; is < numStates; is++) {
      for (long i = 0; i < dimensionality; i++) {
        for (long j = 0; j < dimensionality; j++) {
          for (long k = 0; k < dimensionality; k++) {
            for (long l = 0; l < dimensionality; l++) {
              tmpTensor(iCalc, i, j, k, l) +=
                  0.5 *
                  (w(i, j, is) * w(k, l, is) + w(i, l, is) * w(k, j, is)) *
                  A(i) * A(k) * relTimes(0, 0, is);
            }
          }
        }
      }
    }
#pragma omp critical
    for (long i = 0; i < dimensionality; i++) {
      for (long j = 0; j < dimensionality; j++) {
        for (long k = 0; k < dimensionality; k++) {
          for (long l = 0; l < dimensionality; l++) {
            tensordxdxdxd(iCalc, i, j, k, l) += tmpTensor(iCalc, i, j, k, l);
          }
        }
      }
    }
  }
}

void ElectronViscosity::print() {
  if (!mpi->mpiHead())
    return;

  std::string units;
  if (dimensionality == 1) {
    units = "Pa s / m^2";
  } else if (dimensionality == 2) {
    units = "Pa s / m";
  } else {
    units = "Pa s";
  }

  std::cout << "\n";
  std::cout << "Electron Viscosity (" << units << ")\n";
  std::cout << "i, j, k, eta[i,j,k,1,0], eta[i,j,k,1], eta[i,j,k,2]\n";

  double conversion = pow(hBarSi, 2) // momentum is hbar q
                      / pow(distanceRyToSi, dimensionality) // volume conversion
                      * rydbergSi /
                      hBarSi       // conversion time (q^2 v^2 tau = [time])
                      / rydbergSi; // temperature conversion

  for (long iCalc = 0; iCalc < numCalcs; iCalc++) {

    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStat.temperature;

    std::cout << std::fixed;
    std::cout.precision(2);
    std::cout << "Temperature: " << temp * temperatureAuToSi << " (K)\n";
    std::cout.precision(5);
    std::cout << std::scientific;
    for (long i = 0; i < dimensionality; i++) {
      for (long j = 0; j < dimensionality; j++) {
        for (long k = 0; k < dimensionality; k++) {
          std::cout << i << " " << j << " " << k;
          for (long l = 0; l < dimensionality; l++) {
            std::cout << " " << std::setw(12) << std::right
                      << tensordxdxdxd(iCalc, i, j, k, l) * conversion;
          }
          std::cout << "\n";
        }
      }
    }
    std::cout << std::endl;
  }
}

void ElectronViscosity::outputToJSON(std::string outFileName) {

  if (mpi->mpiHead()) {

    std::string units;
    if (dimensionality == 1) {
      units = "Pa s / m^2";
    } else if (dimensionality == 2) {
      units = "Pa s / m";
    } else {
      units = "Pa s";
    }

    double conversion =
        pow(hBarSi, 2)                        // momentum is hbar q
        / pow(distanceRyToSi, dimensionality) // volume conversion
        * rydbergSi / hBarSi // conversion time (q^2 v^2 tau = [time])
        / rydbergSi;         // temperature conversion

    std::vector<double> temps;
    // this vector mess is of shape (iCalcs, irows, icols, k, l)
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>
        viscosity;

    for (long iCalc = 0; iCalc < numCalcs; iCalc++) {

      // store temperatures
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      double temp = calcStat.temperature;
      temps.push_back(temp * temperatureAuToSi);

      // store viscosity
      std::vector<std::vector<std::vector<std::vector<double>>>> rows;
      for (long i = 0; i < dimensionality; i++) {
        std::vector<std::vector<std::vector<double>>> cols;
        for (long j = 0; j < dimensionality; j++) {
          std::vector<std::vector<double>> ijk;
          for (long k = 0; k < dimensionality; k++) {
            std::vector<double> ijkl;
            for (long l = 0; l < dimensionality; l++) {
              ijkl.push_back(tensordxdxdxd(iCalc, i, j, k, l) * conversion);
            }
            ijk.push_back(ijkl);
          }
          cols.push_back(ijk);
        }
        rows.push_back(cols);
      }
      viscosity.push_back(rows);
    }

    // output to json
    nlohmann::json output;
    output["temperatures"] = temps;
    output["electronViscosity"] = viscosity;
    output["temperatureUnit"] = "K";
    output["electronViscosityUnit"] = units;
    output["particleType"] = "electron";
    std::ofstream o(outFileName);
    o << std::setw(3) << output << std::endl;
    o.close();
  }
}

int ElectronViscosity::whichType() { return is4Tensor; }