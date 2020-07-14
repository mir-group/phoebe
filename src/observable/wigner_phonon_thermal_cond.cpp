#include "wigner_phonon_thermal_cond.h"

#include <iomanip>

#include "constants.h"

WignerPhononThermalConductivity::WignerPhononThermalConductivity(
    StatisticsSweep &statisticsSweep_, Crystal &crystal_,
    BaseBandStructure &bandStructure_, VectorBTE &relaxationTimes)
    : PhononThermalConductivity(statisticsSweep_, crystal_, bandStructure_),
      linewidths(relaxationTimes) {
  int numCalcs = statisticsSweep.getNumCalcs();

  // since we pass in input the relaxation times, we invert them here
  for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
    for (int is = 0; is < bandStructure.getNumStates(); is++) {
      linewidths.data(iCalc, is) = 1. / linewidths.data(iCalc, is);
    }
  }

  wignerCorrection =
      Eigen::Tensor<double, 3>(numCalcs, dimensionality, dimensionality);
  wignerCorrection.setZero();

  auto particle = bandStructure.getParticle();

  int dimensionality = crystal.getDimensionality();

  Eigen::VectorXd norm(numCalcs);
  for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
    double temperature = calcStat.temperature;
    norm(iCalc) = 1. / bandStructure.getNumPoints(true) /
                  crystal.getVolumeUnitCell(dimensionality) / temperature /
                  temperature / 2.;
  }

  int numPoints = bandStructure.getNumPoints();
  for (int iq : mpi->divideWorkIter(numPoints)) {
    State state = bandStructure.getState(iq);
    auto velocities = state.getVelocities();
    auto energies = state.getEnergies();

    int numBands = energies.size();

    Eigen::MatrixXd bose(numCalcs, numBands);
    for (int ib1 = 0; ib1 < numBands; ib1++) {
      for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
        auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
        double temperature = calcStat.temperature;
        bose(iCalc, ib1) = particle.getPopulation(energies(ib1), temperature);
      }
    }

    for (int ib1 = 0; ib1 < numBands; ib1++) {
      for (int ib2 = 0; ib2 < numBands; ib2++) {
        if (ib1 == ib2) continue;

        int is1 = bandStructure.getIndex(WavevectorIndex(iq), BandIndex(ib1));
        int is2 = bandStructure.getIndex(WavevectorIndex(iq), BandIndex(ib2));

        for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
          for (int ic1 = 0; ic1 < dimensionality; ic1++) {
            for (int ic2 = 0; ic2 < dimensionality; ic2++) {
              double num =
                  energies(ib1) * bose(iCalc, ib1) * (bose(iCalc, ib1) + 1.) +
                  energies(ib2) * bose(iCalc, ib2) * (bose(iCalc, ib2) + 1.);
              double vel =
                  (velocities(ib1, ib2, ic1) * velocities(ib2, ib1, ic2))
                      .real();
              double den =
                  4. * pow(energies(ib1) - energies(ib2), 2) +
                  pow(linewidths.data(iCalc, is1) + linewidths.data(iCalc, is2),
                      2);
              wignerCorrection(iCalc, ic1, ic2) +=
                  (energies(ib1) + energies(ib2)) * vel * num / den *
                  (linewidths.data(iCalc, is1) + linewidths.data(iCalc, is2)) *
                  norm(iCalc);
            }
          }
        }
      }
    }
  }
  mpi->allReduceSum(&wignerCorrection);
}

// copy constructor
WignerPhononThermalConductivity::WignerPhononThermalConductivity(
    const WignerPhononThermalConductivity &that)
    : PhononThermalConductivity(that),
      linewidths(that.linewidths),
      wignerCorrection(that.wignerCorrection) {}

// copy assigmnent
WignerPhononThermalConductivity &WignerPhononThermalConductivity::operator=(
    const WignerPhononThermalConductivity &that) {
  PhononThermalConductivity::operator=(that);
  if (this != &that) {
    linewidths = that.linewidths;
    wignerCorrection = that.wignerCorrection;
  }
  return *this;
}

void WignerPhononThermalConductivity::calcFromPopulation(VectorBTE &n) {
  PhononThermalConductivity::calcFromPopulation(n);
  tensordxd += wignerCorrection;
}

void WignerPhononThermalConductivity::calcVariational(VectorBTE &af,
                                                      VectorBTE &f,
                                                      VectorBTE &scalingCG) {
  PhononThermalConductivity::calcVariational(af, f, scalingCG);
  tensordxd += wignerCorrection;
}

void WignerPhononThermalConductivity::calcFromRelaxons(
    SpecificHeat &specificHeat, VectorBTE &relaxonV,
    VectorBTE &relaxationTimes) {
  PhononThermalConductivity::calcFromRelaxons(specificHeat, relaxonV,
                                              relaxationTimes);
  tensordxd += wignerCorrection;
}

void WignerPhononThermalConductivity::print() {
  if (!mpi->mpiHead()) return;  // debugging now

  std::string units;
  if (dimensionality == 1) {
    units = "W m / K";
  } else if (dimensionality == 2) {
    units = "W / K";
  } else {
    units = "W / m / K";
  }

  std::cout << "\n";
  std::cout << "Wigner Thermal Conductivity (" << units << ")\n";

  for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStat.temperature;

    std::cout << std::fixed;
    std::cout.precision(2);
    std::cout << "Temperature: " << temp * temperatureAuToSi << " (K)\n";
    std::cout.precision(5);
    for (long i = 0; i < dimensionality; i++) {
      std::cout << "  " << std::scientific;
      for (long j = 0; j < dimensionality; j++) {
        std::cout << " " << std::setw(13) << std::right;
        std::cout << tensordxd(iCalc, i, j) * thConductivityAuToSi;
      }
      std::cout << "\n";
    }
    std::cout << std::endl;
  }
}
