#include "phonon_viscosity.h"
#include "constants.h"
#include "mpiHelper.h"
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>

PhononViscosity::PhononViscosity(Context &context_,
                                 StatisticsSweep &statisticsSweep_,
                                 Crystal &crystal_,
                                 BaseBandStructure &bandStructure_)
    : Observable(context_, statisticsSweep_, crystal_),
      bandStructure(bandStructure_) {

  tensordxdxdxd = Eigen::Tensor<double, 5>(numCalculations, dimensionality, dimensionality, dimensionality, dimensionality);
  tensordxdxdxd.setZero();
}

// copy constructor
PhononViscosity::PhononViscosity(const PhononViscosity &that)
    : Observable(that), bandStructure(that.bandStructure) {}

// copy assignment
PhononViscosity &PhononViscosity::operator=(const PhononViscosity &that) {
  Observable::operator=(that);
  if (this != &that) {
    bandStructure = that.bandStructure;
  }
  return *this;
}

void PhononViscosity::calcRTA(VectorBTE &tau) {
  double norm = 1. / context.getQMesh().prod() /
                crystal.getVolumeUnitCell(dimensionality);

  auto particle = bandStructure.getParticle();
  tensordxdxdxd.setZero();

  auto excludeIndices = tau.excludeIndices;

#pragma omp parallel default(none) shared(tensordxdxdxd,bandStructure,excludeIndices,numCalculations,statisticsSweep,particle,norm,tau)
  {
    Eigen::Tensor<double, 5> tmpTensor = tensordxdxdxd.constant(0.);

#pragma omp for nowait
    for (int is : bandStructure.parallelIrrStateIterator()) {

      auto isIdx = StateIndex(is);
      int iBte = bandStructure.stateToBte(isIdx).get();

      // skip the acoustic phonons
      if (std::find(excludeIndices.begin(), excludeIndices.end(), iBte) !=
          excludeIndices.end())
        continue;

      auto en = bandStructure.getEnergy(isIdx);
      auto velIrr = bandStructure.getGroupVelocity(isIdx);
      auto qIrr = bandStructure.getWavevector(isIdx);

      auto rotations = bandStructure.getRotationsStar(isIdx);
      for (const Eigen::Matrix3d& rotation : rotations) {

        Eigen::Vector3d q = rotation * qIrr;
        Eigen::Vector3d vel = rotation * velIrr;

        for (int iCalc = 0; iCalc < numCalculations; iCalc++) {

          auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
          double temperature = calcStat.temperature;
          double chemPot = calcStat.chemicalPotential;
          double boseP1 = particle.getPopPopPm1(en, temperature, chemPot);

          for (int i = 0; i < dimensionality; i++) {
            for (int j = 0; j < dimensionality; j++) {
              for (int k = 0; k < dimensionality; k++) {
                for (int l = 0; l < dimensionality; l++) {
                  tmpTensor(iCalc, i, j, k, l) +=
                      q(i) * vel(j) * q(k) * vel(l) * boseP1 *
                      tau(iCalc, 0, iBte) / temperature * norm;
                }
              }
            }
          }
        }
      }
    }
#pragma omp critical
    for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
      for (int i = 0; i < dimensionality; i++) {
        for (int j = 0; j < dimensionality; j++) {
          for (int k = 0; k < dimensionality; k++) {
            for (int l = 0; l < dimensionality; l++) {
              tensordxdxdxd(iCalc, i, j, k, l) += tmpTensor(iCalc, i, j, k, l);
            }
          }
        }
      }
    }
  }
  mpi->allReduceSum(&tensordxdxdxd);
}

void PhononViscosity::calcFromRelaxons(Vector0 &vector0,
                                       Eigen::VectorXd &eigenvalues,
                                       PhScatteringMatrix &sMatrix,
                                       ParallelMatrix<double> &eigenvectors) {

  if (numCalculations > 1) {
    Error("Viscosity for relaxons only for 1 temperature");
  }

  double volume = crystal.getVolumeUnitCell(dimensionality);
  int numStates = bandStructure.getNumStates();
  auto particle = bandStructure.getParticle();

  // to simplify, here I do everything considering there is a single
  // temperature (due to memory constraints)

  Eigen::VectorXd A(dimensionality);
  A.setZero();

  int iCalc = 0;
  auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
  double temp = calcStat.temperature;
  double chemPot = calcStat.chemicalPotential;

  for (int is : bandStructure.parallelStateIterator()) {
    auto isIdx = StateIndex(is);
    auto en = bandStructure.getEnergy(isIdx);
    double boseP1 = particle.getPopPopPm1(en, temp, chemPot); // = n(n+1)
    auto q = bandStructure.getWavevector(isIdx);
    for (int iDim = 0; iDim < dimensionality; iDim++) {
      A(iDim) += boseP1 * q(iDim) * q(iDim);
    }
  }
  A /= temp * context.getQMesh().prod() * volume;
  mpi->allReduceSum(&A);

  VectorBTE driftEigenvector(statisticsSweep, bandStructure, 3);
  for (int is : bandStructure.parallelStateIterator()) {
    auto isIdx = StateIndex(is);
    auto en = bandStructure.getEnergy(isIdx);
    double boseP1 = particle.getPopPopPm1(en, temp, chemPot); // = n(n+1)
    auto q = bandStructure.getWavevector(isIdx);
    for (auto iDim : {0, 1, 2}) {
      driftEigenvector(0, iDim, is) = q(iDim) * sqrt(boseP1 / temp / A(iDim));
    }
  }
  mpi->allReduceSum(&driftEigenvector.data);

  Eigen::MatrixXd D(3, 3);
  D.setZero();
  //    D = driftEigenvector * sMatrix.dot(driftEigenvector.transpose());
  VectorBTE tmp = sMatrix.dot(driftEigenvector);
  for (int is : bandStructure.parallelStateIterator()) {
    for (auto i : {0, 1, 2}) {
      for (auto j : {0, 1, 2}) {
        D(i, j) += driftEigenvector(0, i, is) * tmp(0, j, is);
      }
    }
  }
  D /= volume * context.getQMesh().prod();
  mpi->allReduceSum(&D);

  Eigen::Tensor<double, 3> tmpDriftEigvecs(3, 3, numStates);
  tmpDriftEigvecs.setZero();
  Eigen::MatrixXd W(3, 3);
  W.setZero();
  for (int is : bandStructure.parallelStateIterator()) {
    auto isIdx = StateIndex(is);
    auto v = bandStructure.getGroupVelocity(isIdx);
    for (auto i : {0, 1, 2}) {
      for (auto j : {0, 1, 2}) {
        tmpDriftEigvecs(i, j, is) = driftEigenvector(0, j, is) * v(i);
        W(i, j) += vector0(0, 0, is) * v(i) * driftEigenvector(0, j, is);
      }
    }
  }
  W /= volume * context.getQMesh().prod();
  mpi->allReduceSum(&W);
  mpi->allReduceSum(&tmpDriftEigvecs);

  Eigen::Tensor<double, 3> w(3, 3, numStates);
  w.setZero();
  for (auto i : {0, 1, 2}) {
    for (auto j : {0, 1, 2}) {
      Eigen::VectorXd x(numStates);
      for (int is = 0; is < numStates; is++) {
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

      for (int is = 0; is < numStates; is++) {
        w(i, j, is) = x2[is] / sqrt(volume) / sqrt(context.getQMesh().prod());
      }
      // note: in Eq. 9 of PRX, w is normalized by V*N_q
      // here however I normalize the eigenvectors differently:
      // \sum_state theta_s^2 = 1, instead of 1/VN_q \sum_state theta_s^2 = 1
    }
  }

  // Eq. 9, Simoncelli PRX (2019)
  tensordxdxdxd.setZero();
#pragma omp parallel default(none) shared(tensordxdxdxd,bandStructure,dimensionality,eigenvalues,w,A,iCalc)
  {
    Eigen::Tensor<double, 5> tmpTensor = tensordxdxdxd.constant(0.);
#pragma omp for nowait
    for (int is : bandStructure.parallelStateIterator()) {
      if (eigenvalues(is) <= 0.) { // avoid division by zero
        continue;
      }
      for (int i = 0; i < dimensionality; i++) {
        for (int j = 0; j < dimensionality; j++) {
          for (int k = 0; k < dimensionality; k++) {
            for (int l = 0; l < dimensionality; l++) {
              tmpTensor(iCalc, i, j, k, l) +=
                  0.5 *
                  (w(i, j, is) * w(k, l, is) + w(i, l, is) * w(k, j, is)) *
                  A(i) * A(k) / eigenvalues(is);
            }
          }
        }
      }
    }
#pragma omp critical
    for (int i = 0; i < dimensionality; i++) {
      for (int j = 0; j < dimensionality; j++) {
        for (int k = 0; k < dimensionality; k++) {
          for (int l = 0; l < dimensionality; l++) {
            tensordxdxdxd(iCalc, i, j, k, l) += tmpTensor(iCalc, i, j, k, l);
          }
        }
      }
    }
  }
  mpi->allReduceSum(&tensordxdxdxd);
}

void PhononViscosity::print() {
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
  std::cout << "Thermal Viscosity (" << units << ")\n";
  std::cout << "i, j, k, eta[i,j,k,0], eta[i,j,k,1], eta[i,j,k,2]\n";

  double conversion = pow(hBarSi, 2) // momentum is hBar q
                      / pow(distanceRyToSi, dimensionality) // volume conversion
                      / twoPi // because angular frequencies
                      * rydbergSi /
                      hBarSi       // conversion time (q^2 v^2 tau = [time])
                      / rydbergSi; // temperature conversion

  for (int iCalc = 0; iCalc < numCalculations; iCalc++) {

    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStat.temperature;

    std::cout << std::fixed;
    std::cout.precision(2);
    std::cout << "Temperature: " << temp * temperatureAuToSi << " (K)\n";
    std::cout.precision(5);
    std::cout << std::scientific;
    for (int i = 0; i < dimensionality; i++) {
      for (int j = 0; j < dimensionality; j++) {
        for (int k = 0; k < dimensionality; k++) {
          std::cout << i << " " << j << " " << k;
          for (int l = 0; l < dimensionality; l++) {
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

void PhononViscosity::outputToJSON(const std::string& outFileName) {

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
        pow(hBarSi, 2)                        // momentum is hBar q
        / pow(distanceRyToSi, dimensionality) // volume conversion
        / twoPi // because angular frequencies
        * rydbergSi / hBarSi // conversion time (q^2 v^2 tau = [time])
        / rydbergSi;         // temperature conversion

    std::vector<double> temps;
    // this vector mess is of shape (iCalculations, iRows, iColumns, k, l)
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>
        viscosity;

    for (int iCalc = 0; iCalc < numCalculations; iCalc++) {

      // store temperatures
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      double temp = calcStat.temperature;
      temps.push_back(temp * temperatureAuToSi);

      // store viscosity
      std::vector<std::vector<std::vector<std::vector<double>>>> rows;
      for (int i = 0; i < dimensionality; i++) {
        std::vector<std::vector<std::vector<double>>> cols;
        for (int j = 0; j < dimensionality; j++) {
          std::vector<std::vector<double>> ijk;
          for (int k = 0; k < dimensionality; k++) {
            std::vector<double> ijkl;
            for (int l = 0; l < dimensionality; l++) {
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
    output["phononViscosity"] = viscosity;
    output["temperatureUnit"] = "K";
    output["phononViscosityUnit"] = units;
    output["particleType"] = "phonon";
    std::ofstream o(outFileName);
    o << std::setw(3) << output << std::endl;
    o.close();
  }
}

int PhononViscosity::whichType() { return is4Tensor; }
