#include "onsager.h"
#include "constants.h"
#include "io.h"
#include "mpiHelper.h"
#include "particle.h"
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>

OnsagerCoefficients::OnsagerCoefficients(StatisticsSweep &statisticsSweep_,
                                         Crystal &crystal_,
                                         BaseBandStructure &bandStructure_,
                                         Context &context_)
    : statisticsSweep(statisticsSweep_), crystal(crystal_),
      bandStructure(bandStructure_), context(context_) {
  dimensionality = crystal.getDimensionality();

  if (context.getHasSpinOrbit()) {
    spinFactor = 1.;
  } else {
    // TODO: for spin polarized calculations, we will have to set this
    // to 1.
    spinFactor = 2.;
  }

  numCalculations = statisticsSweep.getNumCalculations();

  sigma.resize(numCalculations, dimensionality, dimensionality);
  seebeck.resize(numCalculations, dimensionality, dimensionality);
  kappa.resize(numCalculations, dimensionality, dimensionality);
  mobility.resize(numCalculations, dimensionality, dimensionality);
  sigma.setZero();
  seebeck.setZero();
  kappa.setZero();
  mobility.setZero();
  LEE.resize(numCalculations, dimensionality, dimensionality);
  LTE.resize(numCalculations, dimensionality, dimensionality);
  LET.resize(numCalculations, dimensionality, dimensionality);
  LTT.resize(numCalculations, dimensionality, dimensionality);
  LEE.setZero();
  LTE.setZero();
  LET.setZero();
  LTT.setZero();
}

// copy constructor
OnsagerCoefficients::OnsagerCoefficients(const OnsagerCoefficients &that)
    : statisticsSweep(that.statisticsSweep), crystal(that.crystal),
      bandStructure(that.bandStructure), context(that.context),
      dimensionality(that.dimensionality), spinFactor(that.spinFactor),
      numCalculations(that.numCalculations), sigma(that.sigma), seebeck(that.seebeck),
      kappa(that.kappa), mobility(that.mobility), LEE(that.LEE), LET(that.LET),
      LTE(that.LTE), LTT(that.LTT) {}

// copy assignment
OnsagerCoefficients &
OnsagerCoefficients::operator=(const OnsagerCoefficients &that) {
  if (this != &that) {
    statisticsSweep = that.statisticsSweep;
    crystal = that.crystal;
    bandStructure = that.bandStructure;
    context = that.context;
    dimensionality = that.dimensionality;
    spinFactor = that.spinFactor;
    numCalculations = that.numCalculations;
    sigma = that.sigma;
    seebeck = that.seebeck;
    kappa = that.kappa;
    mobility = that.mobility;
    LEE = that.LEE;
    LET = that.LET;
    LTE = that.LTE;
    LTT = that.LTT;
  }
  return *this;
}

void OnsagerCoefficients::calcFromEPA(
    VectorEPA &scatteringRates,
    Eigen::Tensor<double, 3> &energyProjVelocity, Eigen::VectorXd &energies) {

  Particle particle(Particle::electron);
  double factor = spinFactor / pow(twoPi, dimensionality);
  double energyStep = energies(1) - energies(0);

  LEE.setZero();
  LET.setZero();
  LTE.setZero();
  LTT.setZero();
  for (int iCalc = 0; iCalc < numCalculations; ++iCalc) {
    double chemPot = statisticsSweep.getCalcStatistics(iCalc).chemicalPotential;
    double temp = statisticsSweep.getCalcStatistics(iCalc).temperature;
    for (int iBeta = 0; iBeta < dimensionality; ++iBeta) {
      for (int iAlpha = 0; iAlpha < dimensionality; ++iAlpha) {
        for (int iEnergy = 0; iEnergy < energies.size(); ++iEnergy) {

          double pop = particle.getPopPopPm1(energies(iEnergy), temp, chemPot);
          double en = energies(iEnergy) - chemPot;
          if (scatteringRates.data(iCalc, iEnergy) <= 1.0e-10 ||
              pop <= 1.0e-20) {
            continue;
          }

          double term = energyProjVelocity(iAlpha, iBeta, iEnergy) /
                        scatteringRates.data(iCalc, iEnergy) * factor * pop *
                        energyStep;

          LEE(iCalc, iAlpha, iBeta) += term / temp;
          LET(iCalc, iAlpha, iBeta) -= term * en / pow(temp, 2);
          LTE(iCalc, iAlpha, iBeta) -= term * en / temp;
          LTT(iCalc, iAlpha, iBeta) -= term * pow(en, 2) / pow(temp, 2);
        }
      }
    }
  }
  calcTransportCoefficients();
}

void OnsagerCoefficients::calcFromCanonicalPopulation(VectorBTE &fE,
                                                      VectorBTE &fT) {
  VectorBTE nE = fE;
  VectorBTE nT = fT;
  nE.canonical2Population(); // n = bose (bose+1) f
  nT.canonical2Population(); // n = bose (bose+1) f
  calcFromPopulation(nE, nT);
}

void OnsagerCoefficients::calcFromSymmetricPopulation(VectorBTE &nE, VectorBTE &nT) {
  VectorBTE nE2 = nE;
  VectorBTE nT2 = nT;

  Particle electron = bandStructure.getParticle();

  for (int is : bandStructure.parallelIrrStateIterator()) {
    StateIndex isIdx(is);
    double energy = bandStructure.getEnergy(isIdx);
    int iBte = bandStructure.stateToBte(isIdx).get();

    for (int iCalc = 0; iCalc < statisticsSweep.getNumCalculations(); iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      double temp = calcStat.temperature;
      double chemPot = calcStat.chemicalPotential;
      double x = sqrt(electron.getPopPopPm1(energy, temp, chemPot));
      for (int i : {0, 1, 2}) {
        nE2(iCalc, i, iBte) *= x;
        nT2(iCalc, i, iBte) *= x;
      }
    }
  }
  calcFromPopulation(nE2, nT2);
}

void OnsagerCoefficients::calcFromPopulation(VectorBTE &nE, VectorBTE &nT) {
  double norm = spinFactor / context.getKMesh().prod() /
                crystal.getVolumeUnitCell(dimensionality);
  LEE.setZero();
  LET.setZero();
  LTE.setZero();
  LTT.setZero();

  auto points = bandStructure.getPoints();

  std::vector<int> states = bandStructure.parallelIrrStateIterator();
  int numStates = states.size();

  Kokkos::View<double***, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      LEE_k(LEE.data(), numCalculations, 3, 3);
  Kokkos::View<double***, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      LET_k(LET.data(), numCalculations, 3, 3);
  Kokkos::View<double***, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      LTE_k(LTE.data(), numCalculations, 3, 3);
  Kokkos::View<double***, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      LTT_k(LTT.data(), numCalculations, 3, 3);
  Kokkos::Experimental::ScatterView<double***, Kokkos::LayoutLeft,
                                    Kokkos::HostSpace> scatter_LEE_k(LEE_k);
  Kokkos::Experimental::ScatterView<double***, Kokkos::LayoutLeft,
                                    Kokkos::HostSpace> scatter_LET_k(LET_k);
  Kokkos::Experimental::ScatterView<double***, Kokkos::LayoutLeft,
                                    Kokkos::HostSpace> scatter_LTE_k(LTE_k);
  Kokkos::Experimental::ScatterView<double***, Kokkos::LayoutLeft,
                                    Kokkos::HostSpace> scatter_LTT_k(LTT_k);
  Kokkos::parallel_for("onsager_coefficients",
                       Kokkos::RangePolicy<Kokkos::HostSpace::execution_space>(0, numStates),
                       [&] (int iis) {

    auto LEE_private = scatter_LEE_k.access();
    auto LET_private = scatter_LET_k.access();
    auto LTE_private = scatter_LTE_k.access();
    auto LTT_private = scatter_LTT_k.access();

    int is = states[iis];

    StateIndex isIdx(is);
    double energy = bandStructure.getEnergy(isIdx);
    Eigen::Vector3d velIrr = bandStructure.getGroupVelocity(isIdx);
    int iBte = bandStructure.stateToBte(isIdx).get();
    auto rotations = bandStructure.getRotationsStar(isIdx);

    for (int iCalc = 0; iCalc < statisticsSweep.getNumCalculations(); iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      double en = energy - calcStat.chemicalPotential;

      for (const Eigen::Matrix3d& r : rotations) {
        Eigen::Vector3d thisNE = Eigen::Vector3d::Zero();
        Eigen::Vector3d thisNT = Eigen::Vector3d::Zero();
        for (int i : {0, 1, 2}) {
          thisNE(i) += nE(iCalc, i, iBte);
          thisNT(i) += nT(iCalc, i, iBte);
        }
        thisNE = r * thisNE;
        thisNT = r * thisNT;
        Eigen::Vector3d vel = r * velIrr;

        for (int i : {0, 1, 2}) {
          for (int j : {0, 1, 2}) {
            LEE_private(iCalc, i, j) += thisNE(i) * vel(j) * norm;
            LET_private(iCalc, i, j) += thisNT(i) * vel(j) * norm;
            LTE_private(iCalc, i, j) += thisNE(i) * vel(j) * en * norm;
            LTT_private(iCalc, i, j) += thisNT(i) * vel(j) * en * norm;
          }
        }
      }
    }
  });
  // safely sum the contribution to tensors from all different threads
  Kokkos::Experimental::contribute(LEE_k, scatter_LEE_k);
  Kokkos::Experimental::contribute(LET_k, scatter_LET_k);
  Kokkos::Experimental::contribute(LTE_k, scatter_LTE_k);
  Kokkos::Experimental::contribute(LTT_k, scatter_LTT_k);
  // same code as above, but without Kokkos
  // in this code, we aren't using OMP threads
  /*
  for (int is : bandStructure.parallelIrrStateIterator()) {
    StateIndex isIdx(is);
    double energy = bandStructure.getEnergy(isIdx);
    Eigen::Vector3d velIrr = bandStructure.getGroupVelocity(isIdx);
    int iBte = bandStructure.stateToBte(isIdx).get();
    auto rotations = bandStructure.getRotationsStar(isIdx);

    for (int iCalc = 0; iCalc < statisticsSweep.getNumCalculations(); iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      double en = energy - calcStat.chemicalPotential;

      for (const Eigen::Matrix3d& r : rotations) {
        Eigen::Vector3d thisNE = Eigen::Vector3d::Zero();
        Eigen::Vector3d thisNT = Eigen::Vector3d::Zero();
        for (int i : {0, 1, 2}) {
          thisNE(i) += nE(iCalc, i, iBte);
          thisNT(i) += nT(iCalc, i, iBte);
        }
        thisNE = r * thisNE;
        thisNT = r * thisNT;
        Eigen::Vector3d vel = r * velIrr;

        for (int i : {0, 1, 2}) {
          for (int j : {0, 1, 2}) {
            LEE(iCalc, i, j) += thisNE(i) * vel(j) * norm;
            LET(iCalc, i, j) += thisNT(i) * vel(j) * norm;
            LTE(iCalc, i, j) += thisNE(i) * vel(j) * en * norm;
            LTT(iCalc, i, j) += thisNT(i) * vel(j) * en * norm;
          }
        }
      }
    }
  }
  */
  // lastly, the states were distributed with MPI
  mpi->allReduceSum(&LEE);
  mpi->allReduceSum(&LTE);
  mpi->allReduceSum(&LET);
  mpi->allReduceSum(&LTT);
  calcTransportCoefficients();
}

void OnsagerCoefficients::calcTransportCoefficients() {
  sigma = LEE;
  mobility = sigma;

  for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
    Eigen::MatrixXd thisLEE(dimensionality, dimensionality);
    Eigen::MatrixXd thisLET(dimensionality, dimensionality);
    Eigen::MatrixXd thisLTE(dimensionality, dimensionality);
    Eigen::MatrixXd thisLTT(dimensionality, dimensionality);
    for (int i = 0; i < dimensionality; i++) {
      for (int j = 0; j < dimensionality; j++) {
        thisLEE(i, j) = LEE(iCalc, i, j);
        thisLET(i, j) = LET(iCalc, i, j);
        thisLTE(i, j) = LTE(iCalc, i, j);
        thisLTT(i, j) = LTT(iCalc, i, j);
      }
    }

    // seebeck = - matmul(L_EE_inv, L_ET)
    Eigen::MatrixXd thisSeebeck = -(thisLEE.inverse()) * thisLET;
    // note: in the unit conversion, I have to consider that S, proportional
    // to 1/e, gets a negative sign coming from the negative electron charge
    // Note that below, the correction on kappa is always positive (L_TE L_ET
    // carries a factor e^2)
    thisSeebeck *= -1;

    // k = L_tt - L_TE L_EE^-1 L_ET
    Eigen::MatrixXd thisKappa = thisLTE * (thisLEE.inverse());
    thisKappa = -(thisLTT - thisKappa * thisLET);
    double doping = abs(statisticsSweep.getCalcStatistics(iCalc).doping);
    doping *= pow(distanceBohrToCm, dimensionality); // from cm^-3 to bohr^-3
    for (int i = 0; i < dimensionality; i++) {
      for (int j = 0; j < dimensionality; j++) {
        seebeck(iCalc, i, j) = thisSeebeck(i, j);
        kappa(iCalc, i, j) = thisKappa(i, j);
        if (doping > 0.) {
          mobility(iCalc, i, j) /= doping;
        }
      }
    }
  }
}

void OnsagerCoefficients::calcFromRelaxons(
    Eigen::VectorXd &eigenvalues, ParallelMatrix<double> &eigenvectors,
    ElScatteringMatrix &scatteringMatrix) {
  int iCalc = 0;
  double chemPot = statisticsSweep.getCalcStatistics(iCalc).chemicalPotential;
  double temp = statisticsSweep.getCalcStatistics(iCalc).temperature;
  auto particle = bandStructure.getParticle();

  VectorBTE nE(statisticsSweep, bandStructure, 3);
  VectorBTE nT(statisticsSweep, bandStructure, 3);

  if (!context.getUseSymmetries()) {

    VectorBTE fE(statisticsSweep, bandStructure, 3);
    VectorBTE fT(statisticsSweep, bandStructure, 3);
    for (auto tup0 : eigenvectors.getAllLocalStates()) {
      int is = std::get<0>(tup0);
      int alfa = std::get<1>(tup0);
      if (eigenvalues(alfa) <= 0.) {
        continue;
      }
      auto isIndex = StateIndex(is);
      double en = bandStructure.getEnergy(isIndex);
      auto vel = bandStructure.getGroupVelocity(isIndex);
      for (int i : {0, 1, 2}) {
        fE(iCalc, i, alfa) += -particle.getDnde(en, temp, chemPot, true)
            * vel(i) * eigenvectors(is, alfa) / eigenvalues(alfa);
        fT(iCalc, i, alfa) += -particle.getDndt(en, temp, chemPot, true)
            * vel(i) * eigenvectors(is, alfa) / eigenvalues(alfa);
      }
    }
    mpi->allReduceSum(&fE.data);
    mpi->allReduceSum(&fT.data);

    for (auto tup0 : eigenvectors.getAllLocalStates()) {
      int is = std::get<0>(tup0);
      int alfa = std::get<1>(tup0);
      for (int i : {0, 1, 2}) {
        nE(iCalc, i, is) += fE(iCalc, i, alfa) * eigenvectors(is, alfa);
        nT(iCalc, i, is) += fT(iCalc, i, alfa) * eigenvectors(is, alfa);
      }
    }
    mpi->allReduceSum(&nE.data);
    mpi->allReduceSum(&nT.data);

  } else { // with symmetries
    Error("Not sure relaxons with symmetries would work");
    Eigen::MatrixXd fE(3, eigenvectors.cols());
    Eigen::MatrixXd fT(3, eigenvectors.cols());
    fE.setZero();
    fT.setZero();
    for (auto tup0 : eigenvectors.getAllLocalStates()) {
      int iMat1 = std::get<0>(tup0);
      int alpha = std::get<1>(tup0);
      auto tup1 = scatteringMatrix.getSMatrixIndex(iMat1);
      BteIndex iBteIndex = std::get<0>(tup1);
      CartIndex dimIndex = std::get<1>(tup1);
      int iDim = dimIndex.get();
      StateIndex isIndex = bandStructure.bteToState(iBteIndex);

      auto vel = bandStructure.getGroupVelocity(isIndex);
      double en = bandStructure.getEnergy(isIndex);
      double dndt = particle.getDndt(en, temp, chemPot);
      double dnde = particle.getDnde(en, temp, chemPot);

      if (eigenvalues(alpha) <= 0.) {
        continue;
      }
      fE(iDim, alpha) +=
          -sqrt(dnde) * vel(iDim) * eigenvectors(iMat1, alpha) / eigenvalues(alpha);
      fT(iDim, alpha) +=
          -sqrt(dndt) * vel(iDim) * eigenvectors(iMat1, alpha) / eigenvalues(alpha);
    }
    mpi->allReduceSum(&fT);
    mpi->allReduceSum(&fE);

    // back rotate to Bloch electron coordinates

    for (auto tup0 : eigenvectors.getAllLocalStates()) {
      int iMat1 = std::get<0>(tup0);
      int alpha = std::get<1>(tup0);
      auto tup1 = scatteringMatrix.getSMatrixIndex(iMat1);
      BteIndex iBteIndex = std::get<0>(tup1);
      CartIndex dimIndex = std::get<1>(tup1);
      int iBte = iBteIndex.get();
      int iDim = dimIndex.get();
      nE(iCalc, iDim, iBte) += fE(iDim, alpha) * eigenvectors(iMat1, alpha);
      nT(iCalc, iDim, iBte) += fT(iDim, alpha) * eigenvectors(iMat1, alpha);
    }
    mpi->allReduceSum(&nE.data);
    mpi->allReduceSum(&nT.data);
  }
  calcFromSymmetricPopulation(nE, nT);
}

void OnsagerCoefficients::print() {
  if (!mpi->mpiHead())
    return;

  std::string unitsSigma, unitsKappa;
  double convSigma, convKappa;
  if (dimensionality == 1) {
    unitsSigma = "S m";
    unitsKappa = "W m / K";
    convSigma = elConductivityAuToSi * rydbergSi * rydbergSi;
    convKappa = thConductivityAuToSi * rydbergSi * rydbergSi;
  } else if (dimensionality == 2) {
    unitsSigma = "S";
    unitsKappa = "W / K";
    convSigma = elConductivityAuToSi * rydbergSi;
    convKappa = thConductivityAuToSi * rydbergSi;
  } else {
    unitsSigma = "S / m";
    unitsKappa = "W / m / K";
    convSigma = elConductivityAuToSi;
    convKappa = thConductivityAuToSi;
  }

  double convMobility = mobilityAuToSi * 100 * 100; // from m^2/Vs to cm^2/Vs
  std::string unitsMobility = "cm^2 / V / s";

  double convSeebeck = thermopowerAuToSi * 1e6;
  std::string unitsSeebeck = "muV / K";

  std::cout << "\n";
  for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStat.temperature;
    double doping = calcStat.doping;
    double chemPot = calcStat.chemicalPotential;

    std::cout << std::fixed;
    std::cout.precision(2);
    std::cout << "iCalc = " << iCalc << ", T = " << temp * temperatureAuToSi;
    std::cout << "\n";

    std::cout.precision(4);
    std::cout << " (K)"
              << ", mu = " << chemPot * energyRyToEv << " (eV)";

    std::cout << std::scientific;
    std::cout << ", n = " << doping << " (cm^-3)" << std::endl;

    std::cout << "Electrical Conductivity (" << unitsSigma << ")\n";
    std::cout.precision(5);
    for (int i = 0; i < dimensionality; i++) {
      std::cout << "  " << std::scientific;
      for (int j = 0; j < dimensionality; j++) {
        std::cout << " " << std::setw(13) << std::right;
        std::cout << sigma(iCalc, i, j) * convSigma;
      }
      std::cout << "\n";
    }
    std::cout << "\n";

    // Note: in metals, one has conductivity without doping
    // and the mobility = sigma / doping-density is ill-defined
    if (abs(doping) > 0.) {
      std::cout << "Carrier mobility (" << unitsMobility << ")\n";
      std::cout.precision(5);
      for (int i = 0; i < dimensionality; i++) {
        std::cout << "  " << std::scientific;
        for (int j = 0; j < dimensionality; j++) {
          std::cout << " " << std::setw(13) << std::right;
          std::cout << mobility(iCalc, i, j) * convMobility;
        }
        std::cout << "\n";
      }
      std::cout << "\n";
    }

    std::cout << "Thermal Conductivity (" << unitsKappa << ")\n";
    std::cout.precision(5);
    for (int i = 0; i < dimensionality; i++) {
      std::cout << "  " << std::scientific;
      for (int j = 0; j < dimensionality; j++) {
        std::cout << " " << std::setw(13) << std::right;
        std::cout << kappa(iCalc, i, j) * convKappa;
      }
      std::cout << "\n";
    }
    std::cout << "\n";

    std::cout << "Seebeck Coefficient (" << unitsSeebeck << ")\n";
    std::cout.precision(5);
    for (int i = 0; i < dimensionality; i++) {
      std::cout << "  " << std::scientific;
      for (int j = 0; j < dimensionality; j++) {
        std::cout << " " << std::setw(13) << std::right;
        std::cout << seebeck(iCalc, i, j) * convSeebeck;
      }
      std::cout << "\n";
    }
    std::cout << std::endl;
  }
}

void OnsagerCoefficients::print(const int &iter) {
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
  for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStat.temperature;
    std::cout << std::fixed;
    std::cout.precision(2);
    std::cout << "T = " << temp * temperatureAuToSi << ", sigma = ";
    std::cout.precision(5);
    for (int i = 0; i < dimensionality; i++) {
      std::cout << std::scientific;
      std::cout << sigma(iCalc, i, i) * elConductivityAuToSi << " ";
    }
    std::cout << "\n";
    std::cout << std::fixed;
    std::cout.precision(2);
    std::cout << "T = " << temp * temperatureAuToSi << ", k = ";
    std::cout.precision(5);
    for (int i = 0; i < dimensionality; i++) {
      std::cout << std::scientific;
      std::cout << kappa(iCalc, i, i) * thConductivityAuToSi << " ";
    }
    std::cout << "\n";
  }
  std::cout << std::endl;
}

void OnsagerCoefficients::outputToJSON(const std::string &outFileName) {
  if (!mpi->mpiHead())
    return;

  std::string unitsSigma, unitsKappa;
  double convSigma, convKappa;
  if (dimensionality == 1) {
    unitsSigma = "S m";
    unitsKappa = "W m / K";
    convSigma = elConductivityAuToSi * rydbergSi * rydbergSi;
    convKappa = thConductivityAuToSi * rydbergSi * rydbergSi;
  } else if (dimensionality == 2) {
    unitsSigma = "S";
    unitsKappa = "W / K";
    convSigma = elConductivityAuToSi * rydbergSi;
    convKappa = thConductivityAuToSi * rydbergSi;
  } else {
    unitsSigma = "S / m";
    unitsKappa = "W / m / K";
    convSigma = elConductivityAuToSi;
    convKappa = thConductivityAuToSi;
  }

  double convMobility = mobilityAuToSi * pow(100., 2); // from m^2/Vs to cm^2/Vs
  std::string unitsMobility = "cm^2 / V / s";

  double convSeebeck = thermopowerAuToSi * 1.0e6;
  std::string unitsSeebeck = "muV / K";

  std::vector<double> temps, dopings, chemPots;
  std::vector<std::vector<std::vector<double>>> sigmaOut;
  std::vector<std::vector<std::vector<double>>> mobilityOut;
  std::vector<std::vector<std::vector<double>>> kappaOut;
  std::vector<std::vector<std::vector<double>>> seebeckOut;
  for (int iCalc = 0; iCalc < numCalculations; iCalc++) {

    // store temperatures
    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStat.temperature;
    temps.push_back(temp * temperatureAuToSi);
    double doping = calcStat.doping;
    dopings.push_back(doping); // output in (cm^-3)
    double chemPot = calcStat.chemicalPotential;
    chemPots.push_back(chemPot * energyRyToEv); // output in eV

    // store the electrical conductivity for output
    std::vector<std::vector<double>> rows;
    for (int i = 0; i < dimensionality; i++) {
      std::vector<double> cols;
      for (int j = 0; j < dimensionality; j++) {
        cols.push_back(sigma(iCalc, i, j) * convSigma);
      }
      rows.push_back(cols);
    }
    sigmaOut.push_back(rows);

    // store the carrier mobility for output
    // Note: in metals, one has conductivity without doping
    // and the mobility = sigma / doping-density is ill-defined
    if (abs(doping) > 0.) {
      rows.clear();
      for (int i = 0; i < dimensionality; i++) {
        std::vector<double> cols;
        for (int j = 0; j < dimensionality; j++) {
          cols.push_back(mobility(iCalc, i, j) * convMobility);
        }
        rows.push_back(cols);
      }
      mobilityOut.push_back(rows);
    }

    // store thermal conductivity for output
    rows.clear();
    for (int i = 0; i < dimensionality; i++) {
      std::vector<double> cols;
      for (int j = 0; j < dimensionality; j++) {
        cols.push_back(kappa(iCalc, i, j) * convKappa);
      }
      rows.push_back(cols);
    }
    kappaOut.push_back(rows);

    // store seebeck coefficient for output
    rows.clear();
    for (int i = 0; i < dimensionality; i++) {
      std::vector<double> cols;
      for (int j = 0; j < dimensionality; j++) {
        cols.push_back(seebeck(iCalc, i, j) * convSeebeck);
      }
      rows.push_back(cols);
    }
    seebeckOut.push_back(rows);
  }

  // output to json
  nlohmann::json output;
  output["temperatures"] = temps;
  output["temperatureUnit"] = "K";
  output["dopingConcentrations"] = dopings;
  output["dopingConcentrationUnit"] = "cm$^{-" + std::to_string(dimensionality) + "}$";
  output["chemicalPotentials"] = chemPots;
  output["chemicalPotentialUnit"] = "eV";
  output["electricalConductivity"] = sigmaOut;
  output["electricalConductivityUnit"] = unitsSigma;
  output["mobility"] = mobilityOut;
  output["mobilityUnit"] = unitsMobility;
  output["electronicThermalConductivity"] = kappaOut;
  output["electronicThermalConductivityUnit"] = unitsKappa;
  output["seebeckCoefficient"] = seebeckOut;
  output["seebeckCoefficientUnit"] = unitsSeebeck;
  output["particleType"] = "electron";
  std::ofstream o(outFileName);
  o << std::setw(3) << output << std::endl;
  o.close();
}

Eigen::Tensor<double, 3> OnsagerCoefficients::getElectricalConductivity() {
  return sigma;
}

Eigen::Tensor<double, 3> OnsagerCoefficients::getThermalConductivity() {
  return kappa;
}

void OnsagerCoefficients::calcVariational(VectorBTE &afE, VectorBTE &afT,
                                          VectorBTE &fE, VectorBTE &fT,
                                          VectorBTE &bE, VectorBTE &bT,
                                          VectorBTE &scalingCG) {
  double norm = spinFactor / context.getKMesh().prod() /
      crystal.getVolumeUnitCell(dimensionality);
  (void) scalingCG;
  int numCalculations = statisticsSweep.getNumCalculations();

  sigma.setConstant(0.);
  kappa.setConstant(0.);

  Eigen::Tensor<double, 3> y1E = sigma.constant(0.);
  Eigen::Tensor<double, 3> y2E = kappa.constant(0.);
  Eigen::Tensor<double, 3> y1T = sigma.constant(0.);
  Eigen::Tensor<double, 3> y2T = kappa.constant(0.);

  std::vector<int> iss = bandStructure.parallelIrrStateIterator();
  int niss = iss.size();

#pragma omp parallel
  {
    Eigen::Tensor<double, 3> x1E(numCalculations, 3, 3);
    Eigen::Tensor<double, 3> x2E(numCalculations, 3, 3);
    Eigen::Tensor<double, 3> x1T(numCalculations, 3, 3);
    Eigen::Tensor<double, 3> x2T(numCalculations, 3, 3);
    x1E.setConstant(0.);
    x2E.setConstant(0.);
    x1T.setConstant(0.);
    x2T.setConstant(0.);

#pragma omp for nowait
    for (int iis=0; iis<niss; ++iis) {
      int is = iss[iis];
      // skip the acoustic phonons
      if (std::find(fE.excludeIndices.begin(), fE.excludeIndices.end(),
                    is) != fE.excludeIndices.end()) {
        continue;
      }

      StateIndex isIndex(is);
      BteIndex iBteIndex = bandStructure.stateToBte(isIndex);
      int isBte = iBteIndex.get();
      auto rots = bandStructure.getRotationsStar(isIndex);

      for (const Eigen::Matrix3d &rot : rots) {

        for (int iCalc = 0; iCalc < numCalculations; iCalc++) {

          auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
          double temp = calcStat.temperature;

          Eigen::Vector3d fRotE, afRotE, bRotE;
          Eigen::Vector3d fRotT, afRotT, bRotT;
          for (int i : {0, 1, 2}) {
            fRotE(i) = fE(iCalc, i, isBte);
            afRotE(i) = afE(iCalc, i, isBte);
            bRotE(i) = bE(iCalc, i, isBte);
            fRotT(i) = fT(iCalc, i, isBte);
            afRotT(i) = afT(iCalc, i, isBte);
            bRotT(i) = bT(iCalc, i, isBte);
          }
          fRotE = rot * fRotE;
          afRotE = rot * afRotE;
          bRotE = rot * bRotE;

          fRotT = rot * fRotT;
          afRotT = rot * afRotT;
          bRotT = rot * bRotT;

          for (int i : {0, 1, 2}) {
            for (int j : {0, 1, 2}) {
              x1E(iCalc, i, j) += fRotE(i) * afRotE(j) * norm * temp;
              x2E(iCalc, i, j) += fRotE(i) * bRotE(j) * norm * temp;
              x1T(iCalc, i, j) += fRotT(i) * afRotT(j) * norm * temp * temp;
              x2T(iCalc, i, j) += fRotT(i) * bRotT(j) * norm * temp * temp;
            }
          }
        }
      }
    }

#pragma omp critical
    for (int j = 0; j < dimensionality; j++) {
      for (int i = 0; i < dimensionality; i++) {
        for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
          y1E(iCalc, i, j) += x1E(iCalc, i, j);
          y2E(iCalc, i, j) += x2E(iCalc, i, j);
          y1T(iCalc, i, j) += x1T(iCalc, i, j);
          y2T(iCalc, i, j) += x2T(iCalc, i, j);
        }
      }
    }
  }
  mpi->allReduceSum(&y1E);
  mpi->allReduceSum(&y2E);
  mpi->allReduceSum(&y1T);
  mpi->allReduceSum(&y2T);
  sigma = 2. * y2E - y1E;
  kappa = 2. * y2T - y1T;
}

