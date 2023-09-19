#include "electron_viscosity.h"
#include "constants.h"
#include "mpiHelper.h"
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>

ElectronViscosity::ElectronViscosity(Context &context_, StatisticsSweep &statisticsSweep_,
                                 Crystal &crystal_, BaseBandStructure &bandStructure_)
    : Observable(context_, statisticsSweep_, crystal_), bandStructure(bandStructure_) {

  tensordxdxdxd = Eigen::Tensor<double, 5>(numCalculations, dimensionality, dimensionality, dimensionality, dimensionality);
  tensordxdxdxd.setZero();
}

void ElectronViscosity::calcRTA(VectorBTE &tau) {

  // add a relevant spin factor
  double spinFactor = 2.;
  if (context.getHasSpinOrbit()) {
    spinFactor = 1.;
  }

  double Nk = context.getKMesh().prod();
  double norm = spinFactor / Nk / crystal.getVolumeUnitCell(dimensionality);
  auto particle = bandStructure.getParticle();
  tensordxdxdxd.setZero();
  //auto excludeIndices = tau.excludeIndices; // not used for electrons

  std::vector<int> iss = bandStructure.parallelIrrStateIterator();
  int niss = iss.size();

  Kokkos::View<double*****, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> tensordxdxdxd_k(tensordxdxdxd.data(), numCalculations, dimensionality, dimensionality, dimensionality, dimensionality);
  Kokkos::Experimental::ScatterView<double*****, Kokkos::LayoutLeft, Kokkos::HostSpace> scatter_tensordxdxdxd(tensordxdxdxd_k);

  Kokkos::parallel_for("electron_viscosity", Kokkos::RangePolicy<Kokkos::HostSpace::execution_space>(0, niss), [&] (int iis){

    auto tmpTensor = scatter_tensordxdxdxd.access();
    int is = iss[iis];
    auto isIdx = StateIndex(is);
    int iBte = bandStructure.stateToBte(isIdx).get();

    auto en = bandStructure.getEnergy(isIdx);
    auto velIrr = bandStructure.getGroupVelocity(isIdx);
    auto kIrr = bandStructure.getWavevector(isIdx);

    auto rotations = bandStructure.getRotationsStar(isIdx);
    for (const Eigen::Matrix3d& rotation : rotations) {

      Eigen::Vector3d k = rotation * kIrr;
      k = bandStructure.getPoints().bzToWs(k,Points::cartesianCoordinates);
      Eigen::Vector3d vel = rotation * velIrr;

      for (int iCalc = 0; iCalc < numCalculations; iCalc++) {

        auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
        double kBT = calcStat.temperature;
        double chemPot = calcStat.chemicalPotential;
        double fermiP1 = particle.getPopPopPm1(en, kBT, chemPot);

        for (int i = 0; i < dimensionality; i++) {
          for (int j = 0; j < dimensionality; j++) {
            for (int k = 0; k < dimensionality; k++) {
              for (int l = 0; l < dimensionality; l++) {
                tmpTensor(iCalc, i, j, k, l) +=
                  k(i) * vel(j) * k(k) * vel(l) * fermiP1 * tau(iCalc, 0, iBte) / kBT * norm;
              }
            }
          }
        }
      }
    }
  });
  Kokkos::Experimental::contribute(tensordxdxdxd_k, scatter_tensordxdxdxd);
  mpi->allReduceSum(&tensordxdxdxd);
}


void ElectronViscosity::calcFromRelaxons(Eigen::VectorXd &eigenvalues, ParallelMatrix<double> &eigenvectors) {

  if (numCalculations > 1) {
    Error("Developer error: Relaxons electron viscosity cannot be calculated for more than one T.");
  }

  // NOTE: view phonon viscosity for notes about which equations are calculated here.

  // we decide to skip relaxon states
  // 1) there is a relaxon with zero (or epsilon) eigenvalue -> infinite tau
  // 2) there might be other states with infinite lifetimes, we skip them
  // 3) states which are alpha > numRelaxons, which were not calculated to
  //    save computational expense.

  double volume = crystal.getVolumeUnitCell(dimensionality);
  auto particle = bandStructure.getParticle();
  int numRelaxons = eigenvalues.size();
  double Nk = context.getKMesh().prod();
  double numStates = bandStructure.getNumStates();

  int iCalc = 0; // set to zero because of relaxons
  auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
  double kBT = calcStat.temperature;
  double chemPot = calcStat.chemicalPotential;

  Eigen::Tensor<double, 3> fRelaxons(3, 3, numStates);
  fRelaxons.setZero();
  for (auto tup0 : eigenvectors.getAllLocalStates()) {

    int is = std::get<0>(tup0);
    int alpha = std::get<1>(tup0);
    if (eigenvalues(alpha) <= 0. || alpha >= numRelaxons) { continue; }
    StateIndex isIdx(is);
    Eigen::Vector3d vec = bandStructure.getWavevector(isIdx);
    Eigen::Vector3d vel = bandStructure.getGroupVelocity(isIdx);
    double en = bandStructure.getEnergy(isIdx);
    double pop = particle.getPopPopPm1(en, kBT, chemPot);
    // true sets a sqrt term
    for (int k = 0; k < dimensionality; k++) {
      for (int l = 0; l < dimensionality; l++) {
        fRelaxons(k, l, alpha) += vec(k) * vel(l) * sqrt(pop) / kBT /
                                  eigenvalues(alpha) * eigenvectors(is, alpha);
      }
    }
  }

  // transform from relaxon to electron populations
  Eigen::Tensor<double, 3> f(3, 3, bandStructure.getNumStates());
  f.setZero();
  for (auto tup0 : eigenvectors.getAllLocalStates()) {
    int is = std::get<0>(tup0);
    int alpha = std::get<1>(tup0);
    if (eigenvalues(alpha) <= 0. || alpha >= numRelaxons) { continue; }
    for (int i : {0, 1, 2}) {
      for (int j : {0, 1, 2}) {
        f(i, j, is) += eigenvectors(is, alpha) * fRelaxons(i, j, alpha);
      }
    }
  }

  double norm = 1. / volume / Nk;
  tensordxdxdxd.setZero();
  for (int is : bandStructure.parallelStateIterator()) {

    StateIndex isIdx(is);
    Eigen::Vector3d vec = bandStructure.getWavevector(isIdx);
    Eigen::Vector3d vel = bandStructure.getGroupVelocity(isIdx);
    double en = bandStructure.getEnergy(isIdx);
    double pop = particle.getPopPopPm1(en, kBT, chemPot);

    for (int i = 0; i < dimensionality; i++) {
      for (int j = 0; j < dimensionality; j++) {
        for (int k = 0; k < dimensionality; k++) {
          for (int l = 0; l < dimensionality; l++) {
            // note: the sqrt(pop) is to rescale the population from the symmetrized exact BTE
            tensordxdxdxd(iCalc, i, j, k, l) +=
                0.5 * pop * norm * sqrt(pop) *
                (vec(i) * vel(j) * f(k, l, is) + vec(i) * vel(l) * f(k, j, is));
          }
        }
      }
    }
  }
  mpi->allReduceSum(&tensordxdxdxd);
}

void ElectronViscosity::print() {

  if (!mpi->mpiHead()) return;

  std::string units;
  if (dimensionality == 1) {      units = "Pa s / m^2"; } // 3d
  else if (dimensionality == 2) { units = "Pa s / m";   } // 2d
  else {                          units = "Pa s";       } // 1d

  std::cout << "\n";
  std::cout << "Electron Viscosity (" << units << ")\n";
  std::cout << "i, j, k, eta[i,j,k,0], eta[i,j,k,1], eta[i,j,k,2]\n";

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
                      << tensordxdxdxd(iCalc, i, j, k, l) * viscosityAuToSi;
          }
          std::cout << "\n";
        }
      }
    }
    std::cout << std::endl;
  }
}

// TODO replace with general one
void ElectronViscosity::outputToJSON(const std::string &outFileName) {

  if (!mpi->mpiHead()) return;

  std::string units;
  if (dimensionality == 1) {      units = "Pa s / m^2"; } // 3d
  else if (dimensionality == 2) { units = "Pa s / m";   } // 2d
  else {                          units = "Pa s";       } // 1d

  std::vector<double> temps;
  // this vector mess is of shape (iCalculations, iRows, iColumns, k, l)
  std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> viscosity;

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
            ijkl.push_back(tensordxdxdxd(iCalc, i, j, k, l) * viscosityAuToSi);
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

int ElectronViscosity::whichType() { return is4Tensor; }
