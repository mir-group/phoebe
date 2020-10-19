#include "onsager.h"

#include <iomanip>
#include <fstream>
#include <nlohmann/json.hpp>
#include "constants.h"
#include "mpiHelper.h"

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

  numCalcs = statisticsSweep.getNumCalcs();

  sigma = Eigen::Tensor<double, 3>(numCalcs, dimensionality, dimensionality);
  seebeck = Eigen::Tensor<double, 3>(numCalcs, dimensionality, dimensionality);
  kappa = Eigen::Tensor<double, 3>(numCalcs, dimensionality, dimensionality);
  mobility = Eigen::Tensor<double, 3>(numCalcs, dimensionality, dimensionality);
  sigma.setZero();
  seebeck.setZero();
  kappa.setZero();
  mobility.setZero();
  LEE = Eigen::Tensor<double, 3>(numCalcs, dimensionality, dimensionality);
  LTE = Eigen::Tensor<double, 3>(numCalcs, dimensionality, dimensionality);
  LET = Eigen::Tensor<double, 3>(numCalcs, dimensionality, dimensionality);
  LTT = Eigen::Tensor<double, 3>(numCalcs, dimensionality, dimensionality);
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
      numCalcs(that.numCalcs), sigma(that.sigma), seebeck(that.seebeck),
      kappa(that.kappa), mobility(that.mobility), LEE(that.LEE), LET(that.LET),
      LTE(that.LTE), LTT(that.LTT) {}

// copy assigmnent
OnsagerCoefficients &
OnsagerCoefficients::operator=(const OnsagerCoefficients &that) {
  if (this != &that) {
    statisticsSweep = that.statisticsSweep;
    crystal = that.crystal;
    bandStructure = that.bandStructure;
    context = that.context;
    dimensionality = that.dimensionality;
    spinFactor = that.spinFactor;
    numCalcs = that.numCalcs;
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

void OnsagerCoefficients::calcFromCanonicalPopulation(VectorBTE &fE,
                                                      VectorBTE &fT) {
  VectorBTE nE = fE;
  VectorBTE nT = fT;
  nE.canonical2Population(); // n = bose (bose+1) f
  nT.canonical2Population(); // n = bose (bose+1) f
  calcFromPopulation(nE, nT);
}

void OnsagerCoefficients::calcFromEPA() {}

void OnsagerCoefficients::calcFromPopulation(VectorBTE &nE, VectorBTE &nT) {
  double norm = spinFactor / bandStructure.getNumPoints(true) /
                crystal.getVolumeUnitCell(dimensionality);

  LEE.setZero();
  LET.setZero();
  LTE.setZero();
  LTT.setZero();

  for (long is : bandStructure.parallelStateIterator()) {
    double energy = bandStructure.getEnergy(is);
    Eigen::Vector3d vel = bandStructure.getGroupVelocity(is);

    auto t = bandStructure.getIndex(is);
    auto ik = std::get<0>(t).get();

    auto rotations = bandStructure.getPoints().getRotationsStar(ik);

    for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
      double chemicalPotential =
          statisticsSweep.getCalcStatistics(iCalc).chemicalPotential;
      double en = energy - chemicalPotential;

      for (Eigen::Matrix3d rotation : rotations) {
        Eigen::Vector3d thisNE, thisNT, thisVel;
        thisNE.setZero();
        thisNT.setZero();
        thisVel.setZero();
        for (long i = 0; i < dimensionality; i++) {
          for (long j = 0; j < dimensionality; j++) {
            thisNE(i) += rotation(i, j) * nE(iCalc, j, is);
            thisNT(i) += rotation(i, j) * nT(iCalc, j, is);
            thisVel(i) += rotation(i, j) * vel(j);
          }
        }

        for (int i = 0; i < dimensionality; i++) {
          for (int j = 0; j < dimensionality; j++) {
            LEE(iCalc, i, j) += thisNE(i) * thisVel(j) * norm;
            LET(iCalc, i, j) += thisNT(i) * thisVel(j) * norm;
            LTE(iCalc, i, j) += thisNE(i) * thisVel(j) * en * norm;
            LTT(iCalc, i, j) += thisNT(i) * thisVel(j) * en * norm;
          }
        }
      }
    }
  }
  mpi->allReduceSum(&LEE);
  mpi->allReduceSum(&LTE);
  mpi->allReduceSum(&LET);
  mpi->allReduceSum(&LTT);
  calcTransportCoefficients();
}

void OnsagerCoefficients::calcTransportCoefficients() {
  sigma = LEE;
  mobility = sigma;

  for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
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
    // k = L_tt - L_TE L_EE^-1 L_ET
    Eigen::MatrixXd thisKappa = thisLTE * (thisLEE.inverse());
    thisKappa = thisLTT - thisKappa * thisLET;

    double doping = abs(statisticsSweep.getCalcStatistics(iCalc).doping);
    doping *= pow(distanceBohrToCm, dimensionality); // from cm^-3 to bohr^-3
    for (int i = 0; i < dimensionality; i++) {
      for (int j = 0; j < dimensionality; j++) {
        seebeck(iCalc, i, j) = thisSeebeck(i, j);
        kappa(iCalc, i, j) = thisKappa(i, j);
        if ( doping > 0. ) {
          mobility(iCalc, i, j) /= doping;
        }
      }
    }
  }
}

void OnsagerCoefficients::print() {
  if (!mpi->mpiHead()) return;

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
  std::string unitsMobility = "cm^2 / V s";

  double convSeebeck = thermopowerAuToSi * 10.0e6;
  std::string unitsSeebeck = "muV / K";

  std::cout << "\n";
  for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
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
    if ( abs(doping) > 0. ) {
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
  currentTime = time(NULL);
  // and format the time nicely
  char s[200];
  struct tm *p = localtime(&currentTime);
  strftime(s, 200, "%F, %T", p);

  std::cout << "Iteration: " << iter << " | " << s << "\n";
  for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStat.temperature;
    std::cout << std::fixed;
    std::cout.precision(2);
    std::cout << "T = " << temp * temperatureAuToSi << ", k = ";
    std::cout.precision(5);
    for (int i = 0; i < dimensionality; i++) {
      std::cout << std::scientific;
      std::cout << sigma(iCalc, i, i) * elConductivityAuToSi << " ";
    }
    std::cout << "\n";
  }
  std::cout << std::endl;
}

void OnsagerCoefficients::outputToJSON(std::string outFileName) {
  if(!mpi->mpiHead()) return;

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
  std::string unitsMobility = "cm^2 / V s";

  double convSeebeck = thermopowerAuToSi * 10.0e6;
  std::string unitsSeebeck = "muV / K";

  std::vector<double> temps, dopings, chemPots;
  std::vector<std::vector<std::vector<double>>> sigmaOut;
  std::vector<std::vector<std::vector<double>>> mobilityOut;
  std::vector<std::vector<std::vector<double>>> kappaOut;
  std::vector<std::vector<std::vector<double>>> seebeckOut;
  for (long iCalc = 0; iCalc < numCalcs; iCalc++) {

    // store temperatures
    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStat.temperature;
    temps.push_back(temp * temperatureAuToSi);
    double doping = calcStat.doping;
    dopings.push_back(doping); // output in (cm^-3)
    double chemPot = calcStat.chemicalPotential;
    chemPots.push_back(chemPot * energyRyToEv); // output in eV

    //std::cout << "Electrical Conductivity (" << unitsSigma << ")\n";
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
    if ( abs(doping) > 0. ) {
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

