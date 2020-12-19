#include "electron_viscosity.h"
#include "constants.h"
#include "mpiHelper.h"
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>

ElectronViscosity::ElectronViscosity(Context &context_,
                                     StatisticsSweep &statisticsSweep_,
                                     Crystal &crystal_,
                                     BaseBandStructure &bandStructure_)
    : Observable(context_, statisticsSweep_, crystal_),
      bandStructure(bandStructure_) {

  tensordxdxdxd = Eigen::Tensor<double, 5>(
      numCalcs, dimensionality, dimensionality, dimensionality, dimensionality);
  tensordxdxdxd.setZero();
}

// copy constructor
ElectronViscosity::ElectronViscosity(const ElectronViscosity &that)
    : Observable(that), bandStructure(that.bandStructure) {}

// copy assignment
ElectronViscosity &ElectronViscosity::operator=(const ElectronViscosity &that) {
  Observable::operator=(that);
  if (this != &that) {
    bandStructure = that.bandStructure;
  }
  return *this;
}

void ElectronViscosity::calcRTA(VectorBTE &tau) {
  double spinFactor = 2.;
  if (context.getHasSpinOrbit()) {
    spinFactor = 1.;
  }
  double norm = spinFactor / context.getKMesh().prod() /
                crystal.getVolumeUnitCell(dimensionality);
  auto particle = bandStructure.getParticle();
  tensordxdxdxd.setZero();
  auto excludeIndices = tau.excludeIndices;

#pragma omp parallel
  {
    Eigen::Tensor<double, 5> tmpTensor = tensordxdxdxd.constant(0.);

#pragma omp for nowait
    for (long is : bandStructure.parallelIrrStateIterator()) {
      auto isIdx = StateIndex(is);

      double en = bandStructure.getEnergy(isIdx);
      Eigen::Vector3d velIrr = bandStructure.getGroupVelocity(isIdx);
      Eigen::Vector3d qIrr = bandStructure.getWavevector(isIdx);

      long ibte = bandStructure.stateToBte(isIdx).get();

      // skip the acoustic phonons
      if (std::find(excludeIndices.begin(), excludeIndices.end(), ibte) !=
          excludeIndices.end())
        continue;

      auto rots = bandStructure.getRotationsStar(isIdx);
      for (Eigen::Matrix3d rot : rots) {
        Eigen::Vector3d vel = rot * velIrr;
        Eigen::Vector3d q = rot * qIrr;

        for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
          auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
          double temperature = calcStat.temperature;
          double chemPot = calcStat.chemicalPotential;
          double bosep1 = particle.getPopPopPm1(en, temperature, chemPot);

          for (long l = 0; l < dimensionality; l++) {
            for (long k = 0; k < dimensionality; k++) {
              for (long j = 0; j < dimensionality; j++) {
                for (long i = 0; i < dimensionality; i++) {
                  tmpTensor(iCalc, i, j, k, l) +=
                      q(i) * vel(j) * q(k) * vel(l) * bosep1 *
                      tau(iCalc, 0, ibte) / temperature * norm;
                }
              }
            }
          }
        }
      }
    }
#pragma omp critical
    for (long l = 0; l < dimensionality; l++) {
      for (long k = 0; k < dimensionality; k++) {
        for (long j = 0; j < dimensionality; j++) {
          for (long i = 0; i < dimensionality; i++) {
            for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
              tensordxdxdxd(iCalc, i, j, k, l) += tmpTensor(iCalc, i, j, k, l);
            }
          }
        }
      }
    }
  }
  mpi->allReduceSum(&tensordxdxdxd);
}

void ElectronViscosity::calcFromRelaxons(Eigen::VectorXd &eigenvalues,
                                         ParallelMatrix<double> &eigenvectors) {
  if (numCalcs > 1) {
    Error e("Viscosity for relaxons only for 1 temperature");
  }

  // we decide to skip relaxon states
  // 1) there is a relaxon with zero (or epsilon) eigenvalue -> infinite tau
  // 2) there might be other states with infinite lifetimes, we skip them

  double volume = crystal.getVolumeUnitCell(dimensionality);
  auto particle = bandStructure.getParticle();

  // to simplify, here I do everything considering there is a single
  // temperature (due to memory constraints)

  long iCalc = 0;
  auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
  double temp = calcStat.temperature;
  double chemPot = calcStat.chemicalPotential;

  Eigen::Tensor<double,3> fRelaxons(3, 3, bandStructure.getNumStates());
  fRelaxons.setZero();
  for (auto tup0 : eigenvectors.getAllLocalStates()) {
    long is = std::get<0>(tup0);
    long alpha = std::get<1>(tup0);
    if ( eigenvalues(alpha) <= 0. ) {
      continue;
    }
    StateIndex isIdx(is);
    Eigen::Vector3d vec = bandStructure.getWavevector(isIdx);
    Eigen::Vector3d vel = bandStructure.getGroupVelocity(isIdx);
    double en = bandStructure.getEnergy(isIdx);
    double pop = particle.getPopPopPm1(en, temp, chemPot);
    for (long k = 0; k < dimensionality; k++) {
      for (long l = 0; l < dimensionality; l++) {
        fRelaxons(k, l, alpha) += vec(k) * vel(l) * pop / temp /
                                  eigenvalues(alpha) * eigenvectors(is, alpha);
      }
    }
  }

  Eigen::Tensor<double,3> f(3, 3, bandStructure.getNumStates());
  f.setZero();
  for (auto tup0 : eigenvectors.getAllLocalStates()) {
    long is = std::get<0>(tup0);
    long alpha = std::get<1>(tup0);
    for (int i : {0, 1, 2}) {
      for (int j : {0, 1, 2}) {
        f(i, j, is) += eigenvectors(is, alpha) * fRelaxons(i, j, alpha);
      }
    }
  }

  double norm = 1. / volume / context.getKMesh().prod();
  tensordxdxdxd.setZero();
  for (long is : bandStructure.parallelStateIterator()) {
    StateIndex isIdx(is);
    Eigen::Vector3d vec = bandStructure.getWavevector(isIdx);
    Eigen::Vector3d vel = bandStructure.getGroupVelocity(isIdx);
    double en = bandStructure.getEnergy(isIdx);
    double pop = particle.getPopPopPm1(en, temp, chemPot);

    for (long i = 0; i < dimensionality; i++) {
      for (long j = 0; j < dimensionality; j++) {
        for (long k = 0; k < dimensionality; k++) {
          for (long l = 0; l < dimensionality; l++) {
            tensordxdxdxd(iCalc, i, j, k, l) +=
                0.5 * pop * norm *
                (vec(i) * vel(j) * f(k, l, is) +
                 vec(i) * vel(l) * f(k, j, is));
          }
        }
      }
    }
  }
  mpi->allReduceSum(&tensordxdxdxd);
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
