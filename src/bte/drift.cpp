#include "drift.h"

BulkTDrift::BulkTDrift(StatisticsSweep &statisticsSweep_,
                       BaseBandStructure &bandStructure_,
                       const int &dimensionality_)
    : VectorBTE(statisticsSweep_, bandStructure_, dimensionality_) {
  Particle particle = bandStructure.getParticle();
#pragma omp parallel for default(none) shared(bandStructure,statisticsSweep,particle)
  for (int is : bandStructure.parallelIrrStateIterator()) {
    StateIndex isIdx(is);
    double energy = bandStructure.getEnergy(isIdx);
    Eigen::Vector3d vel = bandStructure.getGroupVelocity(isIdx);
    int iBte = bandStructure.stateToBte(isIdx).get();
    for (int iCalc = 0; iCalc < statisticsSweep.getNumCalculations(); iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      auto chemPot = calcStat.chemicalPotential;
      auto temp = calcStat.temperature;
      for (int i : {0, 1, 2}) {
        operator()(iCalc, i, iBte) =
            particle.getDndt(energy, temp, chemPot) * vel(i);
      }
    }
  }
  mpi->allReduceSum(&data);
}

BulkEDrift::BulkEDrift(StatisticsSweep &statisticsSweep_,
                       BaseBandStructure &bandStructure_,
                       const int &dimensionality_)
    : VectorBTE(statisticsSweep_, bandStructure_, dimensionality_) {
  Particle particle = bandStructure.getParticle();
#pragma omp parallel for default(none) shared(bandStructure,statisticsSweep,particle)
  for (int is : bandStructure.parallelIrrStateIterator()) {
    StateIndex isIdx(is);
    double energy = bandStructure.getEnergy(isIdx);
    Eigen::Vector3d vel = bandStructure.getGroupVelocity(isIdx);
    int iBte = bandStructure.stateToBte(isIdx).get();
    for (int iCalc = 0; iCalc < statisticsSweep.getNumCalculations(); iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      auto chemPot = calcStat.chemicalPotential;
      auto temp = calcStat.temperature;
      double x = particle.getDnde(energy, temp, chemPot);
      for (int i : {0, 1, 2}) {
        // note: this is tuned for electrons
        // EDrift = e v dn/de
        operator()(iCalc, i, iBte) = x * vel(i);
      }
    }
  }
  mpi->allReduceSum(&data);
}

Vector0::Vector0(StatisticsSweep &statisticsSweep_,
                 BaseBandStructure &bandStructure_, SpecificHeat &specificHeat)
    : VectorBTE(statisticsSweep_, bandStructure_, 1) {

  Particle particle = bandStructure.getParticle();
#pragma omp parallel for default(none) shared(bandStructure,statisticsSweep,particle,specificHeat)
  for (int is : bandStructure.parallelIrrStateIterator()) {
    StateIndex isIdx(is);
    double energy = bandStructure.getEnergy(isIdx);
    int iBte = bandStructure.stateToBte(isIdx).get();
    for (int iCalc = 0; iCalc < statisticsSweep.getNumCalculations(); iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      double temp = calcStat.temperature;
      double chemPot = calcStat.chemicalPotential;
      double dnde = particle.getDnde(energy, temp, chemPot);
      // note dn/de = n(n+1)/T  (for bosons)
      auto c = specificHeat.get(iCalc);
      double x = -dnde / temp / c;
      operator()(iCalc, 0, iBte) = std::sqrt(x) * energy;
    }
  }
  mpi->allReduceSum(&data);
}
