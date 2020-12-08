#include "drift.h"

BulkTDrift::BulkTDrift(StatisticsSweep &statisticsSweep_,
                       BaseBandStructure &bandStructure_,
                       const long &dimensionality_)
    : VectorBTE(statisticsSweep_, bandStructure_, dimensionality_) {
  Particle particle = bandStructure.getParticle();
#pragma omp parallel for
  for (long is : bandStructure.parallelIrrStateIterator()) {
    StateIndex isIdx(is);
    double energy = bandStructure.getEnergy(isIdx);
    Eigen::Vector3d vel = bandStructure.getGroupVelocity(isIdx);
    long ibte = bandStructure.stateToBte(isIdx).get();
    for (long iCalc = 0; iCalc < statisticsSweep.getNumCalcs(); iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      auto chemPot = calcStat.chemicalPotential;
      auto temp = calcStat.temperature;
      for (int i : {0, 1, 2}) {
        operator()(iCalc, i, ibte) =
            particle.getDndt(energy, temp, chemPot) * vel(i);
      }
    }
  }
  mpi->allReduceSum(&data);
}

BulkEDrift::BulkEDrift(StatisticsSweep &statisticsSweep_,
                       BaseBandStructure &bandStructure_,
                       const long &dimensionality_)
    : VectorBTE(statisticsSweep_, bandStructure_, dimensionality_) {
  Particle particle = bandStructure.getParticle();
#pragma omp parallel for
  for (long is : bandStructure.parallelIrrStateIterator()) {
    StateIndex isIdx(is);
    double energy = bandStructure.getEnergy(isIdx);
    Eigen::Vector3d vel = bandStructure.getGroupVelocity(isIdx);
    long ibte = bandStructure.stateToBte(isIdx).get();
    for (long iCalc = 0; iCalc < statisticsSweep.getNumCalcs(); iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      auto chemPot = calcStat.chemicalPotential;
      auto temp = calcStat.temperature;
      double x = particle.getDnde(energy, temp, chemPot);
      for (int i : {0, 1, 2}) {
        operator()(iCalc, i, ibte) = - x * vel(i);
      }
    }
  }
  mpi->allReduceSum(&data);
}

Vector0::Vector0(StatisticsSweep &statisticsSweep_,
                 BaseBandStructure &bandStructure_, SpecificHeat &specificHeat)
    : VectorBTE(statisticsSweep_, bandStructure_, 1) {
  Particle particle = bandStructure.getParticle();
#pragma omp parallel for
  for (long is : bandStructure.parallelIrrStateIterator()) {
    StateIndex isIdx(is);
    double energy = bandStructure.getEnergy(isIdx);
    long ibte = bandStructure.stateToBte(isIdx).get();
    for (long iCalc = 0; iCalc < statisticsSweep.getNumCalcs(); iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      double temp = calcStat.temperature;
      double chemPot = calcStat.chemicalPotential;
      double dnde = particle.getDnde(energy, temp, chemPot);
      // note dnde = n(n+1)/T  (for bosons)
      auto c = specificHeat.get(iCalc);
      double x = -dnde / temp / c;
      operator()(iCalc, 0, ibte) = std::sqrt(x) * energy;
      // we use std::sqrt because we overwrote sqrt() in the base class
    }
  }
  mpi->allReduceSum(&data);
}
