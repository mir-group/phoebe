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
#include <cmath>
#include <nlohmann/json.hpp>
#include <vector>

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
  int numEnergies = int((maxEnergy - minEnergy) / energyStep);
  // energies at the centers of energy steps
  Eigen::VectorXd energies(numEnergies);
  for (int i = 0; i < numEnergies; ++i) {
    // add 0.5 to be in the middle of the energy step
    energies(i) = (double(i) + 0.5) * energyStep + minEnergy;
  }

  // Read and setup k-point mesh for interpolating band structure
  Points fullPoints(crystal, context.getKMesh());
  bool withVelocities = true;
  bool withEigenvectors = false;

  if (mpi->mpiHead()) {
    std::cout << "\nBuilding electronic band structure" << std::endl;
  }

  // filter to only the bands relevant to transport
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
      context, statisticsSweep, bandStructure, energies, tetrahedrons, crystal);
  outputToJSON("epa_relaxation_times.json", scatteringRates, statisticsSweep,
               numEnergies, energies);

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
  transCoeffs.outputToJSON("epa_onsager_coefficients.json");
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
    Error("Either chemical potentials or dopings must be set");
  }
}

// TODO what is this function doing
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
  int numStates = bandStructure.getNumStates();
  int numPoints = bandStructure.getNumPoints(true);
  int dim = context.getDimensionality();

  Eigen::Tensor<double, 3> energyProjVelocity(dim, dim, numEnergies);
  energyProjVelocity.setZero();

  if (mpi->mpiHead()) {
    std::cout << "\nCalculating energy projected velocity tensor" << std::endl;
  }

  LoopPrint loopPrint("calculating energy projected velocity", "energies",
                      mpi->divideWorkIter(numEnergies).size());
#pragma omp parallel default(none)                                             \
    shared(mpi, energyProjVelocity, dim, numEnergies, numStates,               \
           bandStructure, tetrahedrons, numPoints, loopPrint, energies)
  {
    Eigen::Tensor<double, 3> privateVel(dim, dim, numEnergies);
    privateVel.setZero();

#pragma omp for nowait
    for (int iEnergy : mpi->divideWorkIter(numEnergies)) {
#pragma omp critical
      {
        loopPrint.update(); // loop print not omp thread safe
      }
      for (int iState = 0; iState != numStates; ++iState) {
        StateIndex isIdx(iState);
        double deltaFunction =
            tetrahedrons.getSmearing(energies(iEnergy), isIdx);
        Eigen::Vector3d velocity = bandStructure.getGroupVelocity(isIdx);
        for (int j = 0; j < dim; ++j) {
          for (int i = 0; i < dim; ++i) {
            // TODO is this correctly doing outer product -- NO it's for sure
            // not! it's taking the product of specific elements, when this is
            // actually a cross product!
            // TODO is tetrahedron delta function doing it's job
            privateVel(i, j, iEnergy) +=
                velocity(i) * velocity(j) * deltaFunction / double(numPoints);
          }
        }
      }
    }
// TODO why are we copying here
#pragma omp critical
    for (int iEnergy : mpi->divideWorkIter(numEnergies)) {
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
    TetrahedronDeltaFunction &tetrahedrons, Crystal &crystal) {

  int numStates = fullBandStructure.getNumStates();

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

  // TODO this seems strangely redundant -- just ask the bands which aprticle we
  // have
  auto particle = Particle(Particle::electron);
  auto phParticle = Particle(Particle::phonon);

  if (particle.isPhonon())
    Error("Electronic band structure has to be provided");

  int numCalcs = statisticsSweep.getNumCalculations();
  int numEnergies = energies.size();
  double energyStep = context.getEpaEnergyStep();

  // calculate the density of states at the energies in energies vector
  Eigen::VectorXd dos(numEnergies);
  {
    LoopPrint loopPrint1("calculating DoS", "energies",
                         mpi->divideWorkIter(numEnergies).size());
    dos.setZero();

#pragma omp parallel for default(none)                                         \
    shared(numEnergies, loopPrint1, dos, tetrahedrons, energies, mpi)
    for (int i : mpi->divideWorkIter(numEnergies)) {
#pragma omp critical
      {
        loopPrint1.update(); // loop print not omp thread safe
      }
      dos(i) = tetrahedrons.getDOS(energies(i));
    }
    mpi->allReduceSum(&dos);
    loopPrint1.close();
  }

  // get vector containing averaged phonon frequencies per mode
  InteractionEpa couplingEpa = InteractionEpa::parseEpaCoupling(context);
  Eigen::VectorXd phEnergies = couplingEpa.getPhEnergies();
  int numPhEnergies = phEnergies.size();

  // phJump, a double, describes # of bin-jumps the electron does after
  // scattering
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

// loop over temperatures and chemical potentials, then loop over energies
#pragma omp parallel default(none)                                             \
    shared(mpi, loopPrint, statisticsSweep, numCalcs, numPhEnergies, phJump,   \
           phParticle, particle, energies, phEnergies, numEnergies, dos,       \
           binSize, minElphEnergy, numElphBins, couplingEpa, twoPi, fullBandStructure, spinFactor, crystal, epaRate)
  {
    Eigen::MatrixXd privateRates(numCalcs, numEnergies);
    privateRates.setZero();

#pragma omp for nowait
    for (int iEnergy : mpi->divideWorkIter(numEnergies)) {
#pragma omp critical
      {
        loopPrint.update(); // loop print not omp thread safe
      }

      // get statistics
      for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
        double temp = statisticsSweep.getCalcStatistics(iCalc).temperature;
        double chemPot =
            statisticsSweep.getCalcStatistics(iCalc).chemicalPotential;

        // loop over phonon frequencies
        for (int iPhFreq = 0; iPhFreq < numPhEnergies; iPhFreq++) {

          // Avoid some index out of bound errors
          if (double(iEnergy) + phJump(iPhFreq) + 1. >= double(numEnergies) ||
              double(iEnergy) - phJump(iPhFreq) - 1. < 0.) {
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
          // TODO the coupling is squared already here, right?
          privateRates(iCalc, iEnergy) +=
              gAbsorption * (nBose + nFermiAbsorption) * dosAbsorption +
              gEmission * (nBose + 1 - nFermiEmission) * dosEmission;
        }
      }
    }
    // TODO do we need this
    privateRates = (twoPi / spinFactor) * privateRates /
                   double(fullBandStructure.getNumPoints(true)) *
                   crystal.getVolumeUnitCell(crystal.getDimensionality());

// TODO again seems like unnecessary copying
#pragma omp critical
    {
      for (int iEnergy = 0; iEnergy < numEnergies; ++iEnergy) {
        for (int iCalc = 0; iCalc < numCalcs; ++iCalc) {
          epaRate.data(iCalc, iEnergy) += privateRates(iCalc, iEnergy);
        }
      }
    }
  }
  mpi->allReduceSum(&epaRate.data);
  loopPrint.close();

  return epaRate;
}

/* helper function to output scattering rates at each energy to JSON */
void TransportEpaApp::outputToJSON(const std::string &outFileName,
                                   BaseVectorBTE &scatteringRates,
                                   StatisticsSweep &statisticsSweep,
                                   int &numEnergies,
                                   Eigen::VectorXd &energiesEPA) {

  if (!mpi->mpiHead())
    return;

  std::string particleType = "electron";
  double energyConversion = energyRyToEv;
  std::string energyUnit = "eV";
  double energyToTime = timeRyToFs;

  // need to store as a vector format with dimensions
  // iCalc, ik. ib, iDim (where iState is unfolded into
  // ik, ib) for the velocities and lifetimes, no dim for energies
  std::vector<std::vector<double>> outTimes;
  std::vector<std::vector<double>> outLinewidths;
  std::vector<std::vector<double>> energies;
  std::vector<double> temps;
  std::vector<double> chemPots;

  for (int iCalc = 0; iCalc < statisticsSweep.getNumCalculations(); iCalc++) {
    auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStatistics.temperature;
    double chemPot = calcStatistics.chemicalPotential;
    temps.push_back(temp * temperatureAuToSi);
    chemPots.push_back(chemPot * energyConversion);

    // containers to hold data in std vectors
    std::vector<double> tempT;
    std::vector<double> tempL;
    std::vector<double> tempE;
    // loop over energy values on which the calculation was done
    for (int iEnergy = 0; iEnergy < numEnergies; ++iEnergy) {

      double ene = energiesEPA(iEnergy);
      double tau = 1. / scatteringRates.data(iCalc, iEnergy);
      double linewidth = 1. / tau;
      tempE.push_back(ene * energyConversion);
      tempT.push_back(tau * energyToTime);
      tempL.push_back(linewidth * energyConversion);
    }
    outTimes.push_back(tempT);
    outLinewidths.push_back(tempL);
    energies.push_back(tempE);
  }

  // output to json
  nlohmann::json output;
  output["temperatures"] = temps;
  output["temperatureUnit"] = "K";
  output["chemicalPotentials"] = chemPots;
  output["linewidths"] = outLinewidths;
  output["linewidthsUnit"] = energyUnit;
  output["relaxationTimes"] = outTimes;
  output["relaxationTimeUnit"] = "fs";
  output["energies"] = energies;
  output["energyUnit"] = energyUnit;
  output["particleType"] = particleType;
  std::ofstream o(outFileName);
  o << std::setw(3) << output << std::endl;
  o.close();
}
