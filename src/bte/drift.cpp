#include "drift.h"

#include <math.h>

BulkTDrift::BulkTDrift(StatisticsSweep &statisticsSweep_,
                       BaseBandStructure &bandStructure_,
                       const long &dimensionality_)
    : VectorBTE(statisticsSweep_, bandStructure_, dimensionality_) {
  Particle particle = bandStructure.getParticle();
  for (long is = 0; is < numStates; is++) {
    double energy = bandStructure.getEnergy(is);
    Eigen::Vector3d vel = bandStructure.getGroupVelocity(is);
    for (long iCalc = 0; iCalc < statisticsSweep.getNumCalcs(); iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      auto chemicalPotential = calcStat.chemicalPotential;
      auto temperature = calcStat.temperature;
      for (int i : {0,1,2}) {
        operator()(iCalc, i, is) = particle.getDndt(energy, temperature, chemicalPotential) * vel(i);
      }
    }
  }
}

BulkEDrift::BulkEDrift(StatisticsSweep &statisticsSweep_,
                       BaseBandStructure &bandStructure_,
                       const long &dimensionality_)
    : VectorBTE(statisticsSweep_, bandStructure_, dimensionality_) {
  Particle particle = bandStructure.getParticle();
  for (long is = 0; is < numStates; is++) {
    double energy = bandStructure.getEnergy(is);
    Eigen::Vector3d vel = bandStructure.getGroupVelocity(is);
    for (long iCalc = 0; iCalc < statisticsSweep.getNumCalcs(); iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      auto chemicalPotential = calcStat.chemicalPotential;
      auto temperature = calcStat.temperature;
      for (int i : {0,1,2}) {
        operator()(iCalc, i, is) = -particle.getDnde(energy, temperature, chemicalPotential) * vel(i);
      }
    }
  }
}

Vector0::Vector0(StatisticsSweep &statisticsSweep_,
                 BaseBandStructure &bandStructure_, SpecificHeat &specificHeat)
    : VectorBTE(statisticsSweep_, bandStructure_, 1) {
  Particle particle = bandStructure.getParticle();
  for (long is = 0; is < numStates; is++) {
    double energy = bandStructure.getEnergy(is);
    for (long iCalc = 0; iCalc < statisticsSweep.getNumCalcs(); iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      double temp = calcStat.temperature;
      double chemPot = calcStat.chemicalPotential;
      double dnde = particle.getDnde(energy, temp, chemPot);
      // note dnde = n(n+1)/T  (for bosons)
      auto c = specificHeat.get(iCalc);
      double x = -dnde / temp / c;
      operator()(iCalc, 0, is) += std::sqrt(x) * energy;
      // we use std::sqrt because we overwrote sqrt() in the base class
    }
  }
}
