#include "drift.h"

#include <math.h>

BulkTDrift::BulkTDrift(StatisticsSweep &statisticsSweep_,
                       BaseBandStructure &bandStructure_,
                       const long &dimensionality_)
    : VectorBTE(statisticsSweep_, bandStructure_, dimensionality_) {
  Particle particle = bandStructure.getParticle();
  for (long is = 0; is < numStates; is++) {
    double energy = bandStructure.getEnergy(is);
    Eigen::Vector3d velocity = bandStructure.getGroupVelocity(is);

    for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
      auto tup = loc2Glob(iCalc);
      auto imu = std::get<0>(tup);
      auto it = std::get<1>(tup);
      auto idim = std::get<2>(tup);
      int idimIndex = idim.get();
      double vel = velocity(idimIndex);
      auto calcStat = statisticsSweep.getCalcStatistics(it, imu);
      auto chemicalPotential = calcStat.chemicalPotential;
      auto temperature = calcStat.temperature;

      double x = particle.getDndt(energy, temperature, chemicalPotential) * vel;
      data(iCalc, is) = x;
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
    Eigen::Vector3d velocity = bandStructure.getGroupVelocity(is);

    for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
      auto tup = loc2Glob(iCalc);
      auto imu = std::get<0>(tup);
      auto it = std::get<1>(tup);
      auto idim = std::get<2>(tup);
      int idimIndex = idim.get();
      double vel = velocity(idimIndex);
      auto calcStat = statisticsSweep.getCalcStatistics(it, imu);
      auto chemicalPotential = calcStat.chemicalPotential;
      auto temperature = calcStat.temperature;

      double x = particle.getDnde(energy, temperature, chemicalPotential) * vel;
      data(iCalc, is) = -x;
    }
  }
}

Vector0::Vector0(StatisticsSweep &statisticsSweep_,
                 BaseBandStructure &bandStructure_, SpecificHeat &specificHeat)
    : VectorBTE(statisticsSweep_, bandStructure_, 1) {
  data.setZero();

  Particle particle = bandStructure.getParticle();

  for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
    auto tup = loc2Glob(iCalc);
    auto imu = std::get<0>(tup);
    auto it = std::get<1>(tup);
    auto calcStat = statisticsSweep.getCalcStatistics(it, imu);
    double temp = calcStat.temperature;
    double chemPot = calcStat.chemicalPotential;

    for (long is = 0; is < numStates; is++) {
      double energy = bandStructure.getEnergy(is);
      double dnde = particle.getDnde(energy, temp, chemPot);
      // note dnde = n(n+1)/T  (for bosons)
      auto c = specificHeat.get(imu, it);
      double x = -dnde / temp / c;
      data(iCalc, is) += std::sqrt(x) * energy;
      // we use std::sqrt because we overwrote sqrt() in the base class
    }
  }
}
