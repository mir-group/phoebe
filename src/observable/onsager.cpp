#include "onsager.h"

#include <iomanip>

#include "constants.h"
#include "mpiHelper.h"
#include "io.h"
#include "particle.h"

OnsagerCoefficients::OnsagerCoefficients(StatisticsSweep &statisticsSweep_,
                                         Crystal &crystal_,
                                         BaseBandStructure &bandStructure_,
                                         Context &context_)
    : statisticsSweep(statisticsSweep_),
      crystal(crystal_),
      bandStructure(bandStructure_),
      context(context_) {
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
  LEE = Eigen::Tensor<double, 3>(numCalcs, dimensionality, dimensionality);
  LTE = Eigen::Tensor<double, 3>(numCalcs, dimensionality, dimensionality);
  LET = Eigen::Tensor<double, 3>(numCalcs, dimensionality, dimensionality);
  LTT = Eigen::Tensor<double, 3>(numCalcs, dimensionality, dimensionality);
  sigma.setZero();
  seebeck.setZero();
  kappa.setZero();
  LEE.setZero();
  LTE.setZero();
  LET.setZero();
  LTT.setZero();
}

// copy constructor
OnsagerCoefficients::OnsagerCoefficients(const OnsagerCoefficients &that)
    : statisticsSweep(that.statisticsSweep),
      crystal(that.crystal),
      bandStructure(that.bandStructure),
      context(that.context) {}

// copy assigmnent
OnsagerCoefficients &OnsagerCoefficients::operator=(
    const OnsagerCoefficients &that) {
  if (this != &that) {
    statisticsSweep = that.statisticsSweep;
    crystal = that.crystal;
    bandStructure = that.bandStructure;
    context = that.context;
  }
  return *this;
}

void OnsagerCoefficients::calcFromCanonicalPopulation(VectorBTE &fE,
                                                      VectorBTE &fT) {
  VectorBTE nE = fE;
  VectorBTE nT = fT;
  nE.canonical2Population();  // n = bose (bose+1) f
  nT.canonical2Population();  // n = bose (bose+1) f
  calcFromPopulation(nE, nT);
}

void OnsagerCoefficients::calcFromEPA(BaseVectorBTE &scatteringRates, Eigen::Tensor<double,3> &energyProjVelocity, Eigen::VectorXd &energies, double &energyStep, Particle &particle) {
    
    double factor = spinFactor/pow(twoPi,3)/(crystal.getVolumeUnitCell(dimensionality));
    
    double LEEtoSi = pow(electronSi,2)*pow(velocityRyToSi,2)*hBarSi/temperatureAuToSi/kBoltzmannSi/pow(bohrRadiusSi,3)/rydbergSi;
    
    double LETtoSi = electronSi*pow(velocityRyToSi,2)*hBarSi/pow(temperatureAuToSi,2)/kBoltzmannSi/pow(bohrRadiusSi,3);
    
    double LTEtoSi = electronSi*pow(velocityRyToSi,2)*hBarSi/temperatureAuToSi/kBoltzmannSi/pow(bohrRadiusSi,3);
    
    double LTTtoSi = pow(velocityRyToSi,2)*hBarSi*rydbergSi/pow(temperatureAuToSi,2)/kBoltzmannSi/pow(bohrRadiusSi,3);
    
    LEE.setZero();
    LET.setZero();
    LTE.setZero();
    LTT.setZero();
    
    LoopPrint loopPrint("to calculate electronic transport properties", "pairs of temperatures and chemical potentials" , numCalcs);
    
    for (long iCalc = 0; iCalc != numCalcs; ++iCalc) {
        
        int numChemPots = statisticsSweep.getNumChemicalPotentials();
        int numTemps = statisticsSweep.getNumTemperatures();
        auto [imu, it] = decompress2Indeces(iCalc, numChemPots, numTemps);
            
        double chemicalPotential =
            statisticsSweep.getCalcStatistics(iCalc).chemicalPotential;
        double temperature =
        statisticsSweep.getCalcStatistics(iCalc).temperature;
        
        for ( long iEnergy = 0; iEnergy != energies.size(); ++iEnergy ) {
        double en = energies(iEnergy) - chemicalPotential;
            
            for (long iBeta = 0; iBeta != dimensionality; ++iBeta) {
                for (long iAlpha = 0; iAlpha != dimensionality; ++iAlpha) {
                    
                    double transportDistFunc = energyProjVelocity(iAlpha,iBeta,iEnergy)*(1/scatteringRates.data(iCalc,iEnergy));
                    
                    LEE(iCalc, iAlpha, iBeta) += -LEEtoSi*factor*transportDistFunc*particle.getDnde(energies(iEnergy),it,imu)*energyStep; //in Siemens/meter (Coulomb^2*sec/kg/meter^3)
                    
                    LET(iCalc, iAlpha, iBeta) += -LETtoSi*factor*transportDistFunc*particle.getDndt(energies(iEnergy),it,imu)*energyStep; //in Coulomb/sec/meter/Kelvin
                    
                    LTE(iCalc,iAlpha,iBeta) += -LTEtoSi*factor*transportDistFunc*particle.getDnde(energies(iEnergy),it,imu)*en*energyStep; //in Coulomb/meter/sec
                    
                    LTT(iCalc, iAlpha, iBeta) +=  LTTtoSi*factor*transportDistFunc*particle.getDndt(energies(iEnergy),it,imu)*en*energyStep; //in Joule/sec/Kelvin/meter
                }
            }
        }
    }
    loopPrint.close();
}

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
    for (long iCalc = 0; iCalc < numCalcs; iCalc++) {

      int numChemPots = statisticsSweep.getNumChemicalPotentials();
      int numTemps = statisticsSweep.getNumTemperatures();
      auto [imu, it] = decompress2Indeces(iCalc, numChemPots, numTemps);

      double chemicalPotential =
          statisticsSweep.getCalcStatistics(iCalc).chemicalPotential;
      double en = energy - chemicalPotential;
      for (long i = 0; i < dimensionality; i++) {
        long icPop = nE.glob2Loc(ChemPotIndex(imu), TempIndex(it), DimIndex(i));
        for (long j = 0; j < dimensionality; j++) {
          LEE(iCalc, i, j) += nE.data(icPop, is) * vel(j) * norm;
          LET(iCalc, i, j) += nT.data(icPop, is) * vel(j) * norm;
          LTE(iCalc, i, j) += nE.data(icPop, is) * vel(j) * en * norm;
          LTT(iCalc, i, j) += nT.data(icPop, is) * vel(j) * en * norm;
        }
      }
    }
  }
  mpi->allReduceSum(&LEE);
  mpi->allReduceSum(&LTE);
  mpi->allReduceSum(&LET);
  mpi->allReduceSum(&LTT);
}

void OnsagerCoefficients::calcTransportCoefficients() {
  sigma = LEE;
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
    Eigen::MatrixXd thisSeebeck = - (thisLEE.inverse()) * thisLET;
    // k = L_tt - L_TE L_EE^-1 L_ET
//    Eigen::MatrixXd thisKappa =
//        thisLTT - thisLTE.dot(thisLEE.inverse().dot(thisLET));
    Eigen::MatrixXd thisKappa = thisLTE * (thisLEE.inverse());
    thisKappa = thisLTT - thisKappa * thisLET;

    for (int i = 0; i < dimensionality; i++) {
      for (int j = 0; j < dimensionality; j++) {
        seebeck(iCalc, i, j) = thisSeebeck(i, j);
        kappa(iCalc, i, j) = thisKappa(i, j);
      }
    }
  }
}

void OnsagerCoefficients::print() {
  if (!mpi->mpiHead()) return;  // debugging now

//  std::string units;
//  if (dimensionality == 1) {
//    units = "W m / K";
//  } else if (dimensionality == 2) {
//    units = "W / K";
//  } else {
//    units = "W / m / K";
//  }
//
//  std::cout << "\n";
//  std::cout << "Thermal Conductivity (" << units << ")\n";
//
//  for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
//    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
//    double temp = calcStat.temperature;
//
//    std::cout << std::fixed;
//    std::cout.precision(2);
//    std::cout << "Temperature: " << temp * temperatureAuToSi << " (K)\n";
//    std::cout.precision(5);
//    for (long i = 0; i < dimensionality; i++) {
//      std::cout << "  " << std::scientific;
//      for (long j = 0; j < dimensionality; j++) {
//        std::cout << " " << std::setw(13) << std::right;
//        std::cout << tensordxd(iCalc, i, j) * thConductivityAuToSi;
//      }
//      std::cout << "\n";
//    }
//    std::cout << "\n";
//  }
}
