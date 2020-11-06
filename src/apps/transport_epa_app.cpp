#include "transport_epa_app.h"
#include "bandstructure.h"
#include "constants.h"
#include "context.h"
#include "delta_function.h"
#include "eigen.h"
#include "exceptions.h"
#include "interaction_epa.h"
#include "io.h"
#include "onsager.h"
#include "particle.h"
#include "qe_input_parser.h"
#include "statistics_sweep.h"
#include "utilities.h"
#include "vector_bte.h"
#include <math.h>

void TransportEpaApp::run(Context &context) {
  // parse QE-xml file
  auto t1 = QEParser::parseElHarmonicFourier(context);
  auto crystal = std::get<0>(t1);
  auto electronH0 = std::get<1>(t1);

  //--------------------------------
  // Setup energy grid

  double minEnergy = context.getFermiLevel() - context.getEpaEnergyRange();
  double maxEnergy = context.getFermiLevel() + context.getEpaEnergyRange();
  double energyStep = context.getEpaEnergyStep();
  // in principle, we should add 1 to account for ends of energy interval
  // i will not do that, because will work with the centers of energy steps
  long numEnergies = long((maxEnergy - minEnergy) / energyStep);
  // energies at the centers of energy steps
  Eigen::VectorXd energies(numEnergies);
  for (long i = 0; i < numEnergies; ++i) {
    // add 0.5 to be in the middle of the energy step
    energies(i) = (i + 0.5) * energyStep + minEnergy;
  }

  // Read and setup k-point mesh for interpolating bandstructure
  FullPoints fullPoints(crystal, context.getKMesh());
  bool withVelocities = true;
  bool withEigenvectors = false;

  if (mpi->mpiHead()) {
    std::cout << "\nBuilding electronic bandstructure" << std::endl;
  }

  electronH0.trimBands(context, minEnergy, maxEnergy);

  // Fourier interpolation of the electronic band structure
  FullBandStructure bandStructure =
      electronH0.populate(fullPoints, withVelocities, withEigenvectors);
  // set temperatures, chemical potentials and carrier concentrations
  StatisticsSweep statisticsSweep(context, &bandStructure);

  Particle particle = bandStructure.getParticle();

  if (mpi->mpiHead()) {
    std::cout << "\nStarting EPA with " << numEnergies << " energies and "
              << bandStructure.getNumStates() << " states" << std::endl;
  }

  //--------------------------------
  // set up tetrahedron method
  TetrahedronDeltaFunction tetrahedrons(bandStructure);

  //--------------------------------
  // Calculate EPA scattering rates
  BaseVectorBTE scatteringRates = getScatteringRates(
      context, statisticsSweep, bandStructure, energies, tetrahedrons);

  //--------------------------------
  // calc EPA velocities
  auto energyProjVelocity =
      calcEnergyProjVelocity(context, bandStructure, energies, tetrahedrons);

  //--------------------------------
  // compute transport coefficients
  if (mpi->mpiHead()) {
    std::cout << "\nComputing transport coefficients" << std::endl;
  }
  OnsagerCoefficients transCoeffs(statisticsSweep, crystal, bandStructure,
                                  context);

  transCoeffs.calcFromEPA(scatteringRates, energyProjVelocity, energies,
                          energyStep, particle);

  transCoeffs.calcTransportCoefficients();
  transCoeffs.print();
}

void TransportEpaApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getEpaFileName(), "epaFileName");
  throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");
  throwErrorIfUnset(context.getQMesh(), "kMesh");

  throwErrorIfUnset(context.getElectronFourierCutoff(),
                    "electronFourierCutoff");
  throwErrorIfUnset(context.getEpaEnergyStep(), "epaElectronStep");
  throwErrorIfUnset(context.getEpaEnergyRange(), "epaElectronRange");
  throwErrorIfUnset(context.getTemperatures(), "temperatures");
  if (context.getDopings().size() == 0 &&
      context.getChemicalPotentials().size() == 0) {
    Error e("Either chemical potentials or dopings must be set");
  }
}

void foldWithinBounds(int &idx, const int &numBins) {
  if (idx < 0) {
    idx = 0;
  }
  if (idx > numBins) {
    idx = numBins - 1;
  }
}

Eigen::Tensor<double, 3> TransportEpaApp::calcEnergyProjVelocity(
    Context &context, FullBandStructure &bandStructure,
    const Eigen::VectorXd &energies, TetrahedronDeltaFunction &tetrahedrons) {

  int numEnergies = energies.size();
  long numStates = bandStructure.getNumStates();
  int numPoints = bandStructure.getNumPoints(true);
  int dim = context.getDimensionality();

  Eigen::Tensor<double, 3> energyProjVelocity(dim, dim, numEnergies);
  energyProjVelocity.setZero();

  if (mpi->mpiHead()) {
    std::cout << "\nCalculating energy projected velocity tensor" << std::endl;
  }

  LoopPrint loopPrint("calculating energy projected velocity",
                      "energies", mpi->divideWorkIter(numEnergies).size());
#pragma omp parallel
  {
    Eigen::Tensor<double, 3> privateVel(dim, dim, numEnergies);
    privateVel.setZero();

#pragma omp for nowait
    for (long iEnergy : mpi->divideWorkIter(numEnergies)) {
      loopPrint.update();
      for (long iState = 0; iState != numStates; ++iState) {
        auto isIndex = StateIndex(iState);
        double deltaFunction =
            tetrahedrons.getSmearing(energies(iEnergy), isIndex);
        Eigen::Vector3d velocity = bandStructure.getGroupVelocity(iState);
        for (int j = 0; j < dim; ++j) {
          for (int i = 0; i < dim; ++i) {
            privateVel(i, j, iEnergy) +=
                velocity(i) * velocity(j) * deltaFunction / numPoints;
          }
        }
      }
    }
#pragma omp critical
    for (long iEnergy : mpi->divideWorkIter(numEnergies)) {
      for (int j = 0; j < dim; ++j) {
        for (int i = 0; i < dim; ++i) {
          energyProjVelocity(i, j, iEnergy) += privateVel(i, j, iEnergy);
        }
      }
    }
  }
  mpi->allReduceSum(&energyProjVelocity);
  loopPrint.close();
  return energyProjVelocity;
}

BaseVectorBTE TransportEpaApp::getScatteringRates(
    Context &context, StatisticsSweep &statisticsSweep,
    FullBandStructure &fullBandStructure, Eigen::VectorXd &energies,
    TetrahedronDeltaFunction &tetrahedrons) {

  long numStates = fullBandStructure.getNumStates();

  /*If constant relaxation time is specified in input, we don't need to
  calculate EPA lifetimes*/
  double constantRelaxationTime = context.getConstantRelaxationTime();
  if (constantRelaxationTime > 0.) {
    BaseVectorBTE crtRate(statisticsSweep, numStates, 1);
    crtRate.setConst(1. / constantRelaxationTime);
    return crtRate;
  }

  auto hasSpinOrbit = context.getHasSpinOrbit();
  int spinFactor = 2;
  if (hasSpinOrbit) {
    spinFactor = 1;
  }

  auto particle = Particle(Particle::electron);
  auto phParticle = Particle(Particle::phonon);

  if (particle.isPhonon())
    Error e("Electronic bandstructure has to be provided");

  long numCalcs = statisticsSweep.getNumCalcs();

  long numEnergies = energies.size();
  double energyStep = context.getEpaEnergyStep();

  // calculate the density of states at the energies in energies vector
  Eigen::VectorXd dos(numEnergies);
  {
    LoopPrint loopPrint1("calculating DoS", "energies",
                         mpi->divideWorkIter(numEnergies).size());
    dos.setZero();
#pragma omp parallel for
    for (long i : mpi->divideWorkIter(numEnergies)) {
      loopPrint1.update();
      dos(i) = tetrahedrons.getDOS(energies(i));
    }
    mpi->allReduceSum(&dos);
    loopPrint1.close();
  }

  // get vector containing averaged phonon frequencies per mode
  InteractionEpa couplingEpa = InteractionEpa::parseEpaCoupling(context);

  Eigen::VectorXd phEnergies = couplingEpa.getPhEnergies();
  int numPhEnergies = phEnergies.size();

  // phJump describes how bin-jumps the electron does after scattering
  // as a double
  Eigen::VectorXd phJump(numPhEnergies);
#pragma omp parallel for
  for (auto i = 0; i != phEnergies.size(); ++i) {
    phJump(i) = phEnergies(i) / energyStep;
  }

  Eigen::VectorXd elphEnergies = couplingEpa.getElEnergies();
  double minElphEnergy = elphEnergies(0);
  int numElphBins = elphEnergies.size();
  double binSize = 1.;
  if (numElphBins > 1) {
    binSize = elphEnergies(1) - elphEnergies(0);
  }

  LoopPrint loopPrint("calculation of EPA scattering rates", "energies",
                      mpi->divideWorkIter(numEnergies).size());

  BaseVectorBTE epaRate(statisticsSweep, numEnergies, 1);

  // loop over temperatures and chemical potentials
  // loop over energies
#pragma omp parallel
  {
    Eigen::MatrixXd privateRates(numCalcs, numEnergies);
    privateRates.setZero();

#pragma omp for nowait
    for (long iEnergy : mpi->divideWorkIter(numEnergies)) {
      loopPrint.update();

      for (long iCalc = 0; iCalc < numCalcs; ++iCalc) {
        double temp = statisticsSweep.getCalcStatistics(iCalc).temperature;
        double chemPot =
            statisticsSweep.getCalcStatistics(iCalc).chemicalPotential;

        // loop over phonon frequencies
        for (int iPhFreq = 0; iPhFreq < numPhEnergies; ++iPhFreq) {

          // Avoid some index out of bound errors
          if (iEnergy + phJump(iPhFreq) + 1 >= numEnergies ||
              iEnergy - phJump(iPhFreq) - 1 < 0) {
            continue;
          }

          // population of phonons, electron after emission/absorption
          double nBose = phParticle.getPopulation(phEnergies(iPhFreq), temp);
          double nFermiAbsorption = particle.getPopulation(
              energies[iEnergy] + phEnergies(iPhFreq), temp, chemPot);
          double nFermiEmission = particle.getPopulation(
              energies[iEnergy] - phEnergies(iPhFreq), temp, chemPot);

          // compute the dos for electron in the final state for the two
          // scatterings mechanisms
          // Note: we do a linear interpolation
          int iJump = (int)phJump(iPhFreq);
          double iInterp = phJump(iPhFreq) - (double)iJump;
          double dosAbsorption = dos(iEnergy + iJump) * (1. - iInterp) +
                                 dos(iEnergy + iJump + 1) * iInterp;
          double dosEmission = dos(iEnergy - iJump - 1) * iInterp +
                               dos(iEnergy - iJump) * (1. - iInterp);

          // find index of the energy in the bins of the elph energies
          int intBinPos =
              int(std::round((energies(iEnergy) - minElphEnergy) / binSize));
          int iAbsInt = int(std::round(
              (energies(iEnergy) + phEnergies(iPhFreq) - minElphEnergy) /
              binSize));
          int iEmisInt = int(std::round(
              (energies(iEnergy) - phEnergies(iPhFreq) - minElphEnergy) /
              binSize));
          // check and fold within bounds:
          foldWithinBounds(intBinPos, numElphBins);
          foldWithinBounds(iAbsInt, numElphBins);
          foldWithinBounds(iEmisInt, numElphBins);

          //------------------------------------
          // estimate strength of el-ph coupling |g|^2

          double gAbsorption =
              couplingEpa.getCoupling(iPhFreq, iAbsInt, intBinPos);
          double gEmission =
              couplingEpa.getCoupling(iPhFreq, iEmisInt, intBinPos);

          //-----------------------------
          // finally, the scattering rate

          privateRates(iCalc, iEnergy) +=
              twoPi / spinFactor * gAbsorption * (nBose + nFermiAbsorption) *
                  dosAbsorption +
              gEmission * (nBose + 1 - nFermiEmission) * dosEmission;
        }
      }
    }

#pragma omp critical
    {
      for (long iEnergy = 0; iEnergy < numEnergies; ++iEnergy) {
        for (long iCalc = 0; iCalc < numCalcs; ++iCalc) {
          epaRate.data(iCalc, iEnergy) += privateRates(iCalc, iEnergy);
        }
      }
    }
  }
  mpi->allReduceSum(&epaRate.data);
  loopPrint.close();

  return epaRate;
}
