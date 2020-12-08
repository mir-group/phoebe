#include "wigner_phonon_thermal_cond.h"

#include <iomanip>

#include "constants.h"

WignerPhononThermalConductivity::WignerPhononThermalConductivity(
    Context &context_, StatisticsSweep &statisticsSweep_, Crystal &crystal_,
    BaseBandStructure &bandStructure_, VectorBTE &relaxationTimes)
    : PhononThermalConductivity(context_, statisticsSweep_, crystal_,
                                bandStructure_),
      smaRelTimes(relaxationTimes) {

  int numCalcs = statisticsSweep.getNumCalcs();

  wignerCorrection =
      Eigen::Tensor<double, 3>(numCalcs, dimensionality, dimensionality);
  wignerCorrection.setZero();

  auto particle = bandStructure.getParticle();
  int dimensionality = crystal.getDimensionality();

  Eigen::VectorXd norm(numCalcs);
  for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
    double temperature = calcStat.temperature;
    norm(iCalc) = 1. / context.getQMesh().prod() /
                  crystal.getVolumeUnitCell(dimensionality) / temperature /
                  temperature / 2.;
  }

  for (long iq = 0; iq < bandStructure.getNumPoints(); iq++) {
    auto iqIdx = WavevectorIndex(iq);
    auto velocities = bandStructure.getVelocities(iqIdx);
    auto energies = bandStructure.getEnergies(iqIdx);
    int numBands = energies.size();

    // calculate bose factors
    Eigen::MatrixXd bose(numCalcs, numBands);
    for (int ib1 = 0; ib1 < numBands; ib1++) {
      for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
        auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
        double temperature = calcStat.temperature;
        bose(iCalc, ib1) = particle.getPopulation(energies(ib1), temperature);

        // exclude acoustic phonons, cutoff at 0.1 cm^-1
        // setting this to zero here causes acoustic ph
        // contribution below to evaluate to zero
        if (energies(ib1) < 0.1 / ryToCmm1) {
          bose(iCalc, ib1) = 0.0;
        }
      }
    }

    // calculate wigner correction
    for (int ib1 = 0; ib1 < numBands; ib1++) {
      for (int ib2 = 0; ib2 < numBands; ib2++) {
        if (ib1 == ib2)
          continue;

        long is1 = bandStructure.getIndex(iqIdx, BandIndex(ib1));
        long is2 = bandStructure.getIndex(iqIdx, BandIndex(ib2));
        auto is1Idx = StateIndex(is1);
        auto is2Idx = StateIndex(is2);
        long ind1 = bandStructure.stateToBte(is1Idx).get();
        long ind2 = bandStructure.stateToBte(is2Idx).get();

        for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
          for (int ic1 = 0; ic1 < dimensionality; ic1++) {
            for (int ic2 = 0; ic2 < dimensionality; ic2++) {

              double num =
                  energies(ib1) * bose(iCalc, ib1) * (bose(iCalc, ib1) + 1.) +
                  energies(ib2) * bose(iCalc, ib2) * (bose(iCalc, ib2) + 1.);
              double vel =
                  (velocities(ib1, ib2, ic1) * velocities(ib2, ib1, ic2))
                      .real();
              double den = 4. * pow(energies(ib1) - energies(ib2), 2) +
                           pow(1. / smaRelTimes(iCalc, 0, ind1) +
                                   1. / smaRelTimes(iCalc, 0, ind2), 2);

              wignerCorrection(iCalc, ic1, ic2) +=
                  (energies(ib1) + energies(ib2)) * vel * num / den *
                  (1. / smaRelTimes(iCalc, 0, ind1) +
                   1. / smaRelTimes(iCalc, 0, ind2)) *
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
    : PhononThermalConductivity(that), smaRelTimes(that.smaRelTimes),
      wignerCorrection(that.wignerCorrection) {}

// copy assigmnent
WignerPhononThermalConductivity &WignerPhononThermalConductivity::operator=(
    const WignerPhononThermalConductivity &that) {
  PhononThermalConductivity::operator=(that);
  if (this != &that) {
    smaRelTimes = that.smaRelTimes;
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
    Context &context, StatisticsSweep &statisticsSweep,
    BaseBandStructure &bandStructure, ParallelMatrix<double> &eigenvectors,
    PhScatteringMatrix &scatteringMatrix, const Eigen::VectorXd &eigenvalues) {
  PhononThermalConductivity::calcFromRelaxons(context, statisticsSweep,
                                              bandStructure, eigenvectors,
                                              scatteringMatrix, eigenvalues);
  tensordxd += wignerCorrection;
}

void WignerPhononThermalConductivity::print() {

  if (!mpi->mpiHead())
    return; // debugging now

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
