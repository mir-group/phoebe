#include "wigner_electron.h"
#include "constants.h"
#include <iomanip>

WignerElCoefficients::WignerElCoefficients(StatisticsSweep &statisticsSweep_,
                                           Crystal &crystal_,
                                           BaseBandStructure &bandStructure_,
                                           Context &context_,
                                           VectorBTE &relaxationTimes)
    : OnsagerCoefficients(statisticsSweep_, crystal_, bandStructure_, context_),
      smaRelTimes(relaxationTimes) {

  correctionLEE.resize(numCalculations, dimensionality, dimensionality);
  correctionLTE.resize(numCalculations, dimensionality, dimensionality);
  correctionLET.resize(numCalculations, dimensionality, dimensionality);
  correctionLTT.resize(numCalculations, dimensionality, dimensionality);
  correctionLEE.setZero();
  correctionLTE.setZero();
  correctionLET.setZero();
  correctionLTT.setZero();

  auto particle = bandStructure.getParticle();

  double norm = spinFactor / context.getKMesh().prod() /
                crystal.getVolumeUnitCell(dimensionality) / 2.;

  Eigen::Tensor<std::complex<double>, 4> fE, fT;

  for (int ik : mpi->divideWorkIter(bandStructure.getNumPoints())) {
    WavevectorIndex ikIdx(ik);
    auto velocities = bandStructure.getVelocities(ikIdx);
    auto energies = bandStructure.getEnergies(ikIdx);
    int numBands = energies.size();

    Eigen::Vector3d k = bandStructure.getWavevector(ikIdx);
    auto t = bandStructure.getRotationToIrreducible(k, Points::cartesianCoordinates);
    int ikIrr = std::get<0>(t);

    // we do the calculation in two steps

    //-----------------------------------------------------------------------
    // Step 1: compute the off-diagonal population term (at given wavevector)

    fE.resize(numBands, numBands, dimensionality, numCalculations);
    fT.resize(numBands, numBands, dimensionality, numCalculations);
    fE.setZero();
    fT.setZero();

    for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      double chemPot = calcStat.chemicalPotential;
      double temp = calcStat.temperature;
      Eigen::VectorXd fermi(numBands);
      Eigen::VectorXd dfdt(numBands);
      for (int ib1 = 0; ib1 < numBands; ib1++) {
        fermi(ib1) = particle.getPopulation(energies(ib1), temp, chemPot);
        dfdt(ib1) = particle.getDndt(energies(ib1), temp, chemPot);
      }

      for (int ib1 = 0; ib1 < numBands; ib1++) {
        for (int ib2 = 0; ib2 < numBands; ib2++) {
          if (ib1 == ib2) {
            continue;
          }
          int is1 = bandStructure.getIndex(WavevectorIndex(ikIrr), BandIndex(ib1));
          int is2 = bandStructure.getIndex(WavevectorIndex(ikIrr), BandIndex(ib2));
          StateIndex is1Idx(is1);
          StateIndex is2Idx(is2);
          int iBte1 = bandStructure.stateToBte(is1Idx).get();
          int iBte2 = bandStructure.stateToBte(is2Idx).get();

          std::complex<double> xC = {1. / smaRelTimes(iCalc, 0, iBte1) +
                                         1. / smaRelTimes(iCalc, 0, iBte2),
                                     2. * (energies(ib1) - energies(ib2))};

          for (int ic1 = 0; ic1 < dimensionality; ic1++) {
            fE(ib1, ib2, ic1, iCalc) = -2. * velocities(ib1, ib2, ic1) / xC *
                                       (fermi(ib1) - fermi(ib2)) /
                                       (energies(ib1) - energies(ib2));
            fT(ib1, ib2, ic1, iCalc) =
                2. * velocities(ib1, ib2, ic1) / xC * (dfdt(ib1) + dfdt(ib2));
          }
        }
      }
    }

    //---------------------------------------------------------------------
    // Step 2: now compute the anti-commutator for the transport coefficient

    for (int ib1 = 0; ib1 < numBands; ib1++) {
      for (int ib2 = 0; ib2 < numBands; ib2++) {
        if (ib1 == ib2)
          continue;
        for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
          double chemicalPotential =
              statisticsSweep.getCalcStatistics(iCalc).chemicalPotential;
          for (int ic1 = 0; ic1 < dimensionality; ic1++) {
            for (int ic2 = 0; ic2 < dimensionality; ic2++) {
              double xE = std::real(
                  velocities(ib1, ib2, ic1) * fE(ib2, ib1, ic2, iCalc) +
                  velocities(ib2, ib1, ic1) * fE(ib1, ib2, ic2, iCalc));
              double xT = std::real(
                  velocities(ib1, ib2, ic1) * fT(ib2, ib1, ic2, iCalc) +
                  velocities(ib2, ib1, ic1) * fT(ib1, ib2, ic2, iCalc));
              correctionLEE(iCalc, ic1, ic2) += norm * xE;
              correctionLET(iCalc, ic1, ic2) -= norm * xT;
              correctionLTE(iCalc, ic1, ic2) -=
                  norm * (energies(ib1) - chemicalPotential) * xE;
              correctionLTT(iCalc, ic1, ic2) -=
                  norm * (energies(ib1) - chemicalPotential) * xT;
            }
          }
        }
      }
    }
  }
  mpi->allReduceSum(&correctionLEE);
  mpi->allReduceSum(&correctionLTE);
  mpi->allReduceSum(&correctionLET);
  mpi->allReduceSum(&correctionLTT);
}

// copy constructor
WignerElCoefficients::WignerElCoefficients(const WignerElCoefficients &that)
    : OnsagerCoefficients(that), smaRelTimes(that.smaRelTimes),
      correctionLEE(that.correctionLEE), correctionLTE(that.correctionLTE),
      correctionLET(that.correctionLET), correctionLTT(that.correctionLTT) {}

// copy assignment
WignerElCoefficients &
WignerElCoefficients::operator=(const WignerElCoefficients &that) {
  OnsagerCoefficients::operator=(that);
  if (this != &that) {
    smaRelTimes = that.smaRelTimes;
    correctionLEE = that.correctionLEE;
    correctionLTE = that.correctionLTE;
    correctionLET = that.correctionLET;
    correctionLTT = that.correctionLTT;
  }
  return *this;
}

void WignerElCoefficients::calcFromPopulation(VectorBTE &nE, VectorBTE &nT) {
  OnsagerCoefficients::calcFromPopulation(nE, nT);
  LEE += correctionLEE;
  LTE += correctionLTE;
  LET += correctionLET;
  LTT += correctionLTT;
  // calcTransportCoefficients is called twice, also in base calcFromPopulation.
  // Could this be improved?
  calcTransportCoefficients();
}

void WignerElCoefficients::print() {
  if (!mpi->mpiHead())
    return;
  std::cout << "Estimates with the Wigner transport equation.\n";
  OnsagerCoefficients::print();
}
