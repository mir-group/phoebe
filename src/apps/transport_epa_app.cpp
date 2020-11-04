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

void TransportEpaApp::run(Context &context) {

  double fermiLevel = context.getFermiLevel();
  if (std::isnan(fermiLevel)) {
    Error e("Fermi energy must be provided for EPA calculation");
  }
  // Read necessary input: xml file of QE.
  // name of xml file should be provided in the input
  // electronFourierCutoff should be provided in the input (should it be the
  // same as encut in DFT?) Crystal crystal(directUnitCell, atomicPositions,
  // atomicSpecies, speciesNames, speciesMasses, dimensionality);
  // ElectronH0Fourier electronH0(crystal, coarsePoints, coarseBandStructure,
  // fourierCutoff);

  auto t1 = QEParser::parseElHarmonicFourier(context);
  auto crystal = std::get<0>(t1);
  auto electronH0 = std::get<1>(t1);

  // Read and setup k-point mesh for interpolating bandstructure
  FullPoints fullPoints(crystal, context.getKMesh());
  bool withVelocities = true;
  bool withEigenvectors = true;

  // Fourier interpolation of the electronic band structure
  FullBandStructure bandStructure =
      electronH0.populate(fullPoints, withVelocities, withEigenvectors);

  Particle particle = bandStructure.getParticle();

  // set temperatures, chemical potentials and carrier concentrations
  StatisticsSweep statisticsSweep(context, &bandStructure);

  //--------------------------------
  // Setup energy grid

  double minEnergy = fermiLevel - context.getEnergyRange();
  double maxEnergy = fermiLevel + context.getEnergyRange();
  double energyStep = context.getEnergyStep();
  // in principle, we should add 1 to account for ends of energy interval
  // i will not do that, because will work with the centers of energy steps
  long numEnergies = long((maxEnergy - minEnergy) / energyStep);
  if (mpi->mpiHead()) {
    std::cout << "Num energies: " << numEnergies << std::endl;
  }
  // energies at the centers of energy steps
  Eigen::VectorXd energies(numEnergies);
  for (long i = 0; i < numEnergies; ++i) {
    // add 0.5 to be in the middle of the energy step
    energies(i) = (i + 0.5) * energyStep + minEnergy;
  }

  //--------------------------------
  // Calculate EPA scattering rates
  BaseVectorBTE scatteringRates =
      getScatteringRates(context, statisticsSweep, bandStructure, energies);

  //--------------------------------
  // calc EPA velocities
  auto energyProjVelocity =
      calcEnergyProjVelocity(context, bandStructure, energies);

  //--------------------------------
  // compute transport coefficients
  OnsagerCoefficients transCoeffs(statisticsSweep, crystal, bandStructure,
                                  context);

  transCoeffs.calcFromEPA(scatteringRates, energyProjVelocity, energies,
                          energyStep, particle);

  transCoeffs.calcTransportCoefficients();
  transCoeffs.print();
}

// void TransportEpaApp::checkRequirements(Context & context) {
//    throwErrorIfUnset(context.getEpaEFileName(), "epaEFileName");
//    throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");
//    throwErrorIfUnset(context.getQMesh(), "kMesh");
//    throwErrorIfUnset(context.getDosMinEnergy(), "dosMinEnergy");
//    throwErrorIfUnset(context.getDosMaxEnergy(), "dosMaxEnergy");
//    throwErrorIfUnset(context.getDosDeltaEnergy(), "dosDeltaEnergy");
//    throwErrorIfUnset(context.getElectronFourierCutoff(),
//                      "electronFourierCutoff");
//}

Eigen::Tensor<double, 3>
TransportEpaApp::calcEnergyProjVelocity(Context &context,
                                        BaseBandStructure &bandStructure,
                                        const Eigen::VectorXd &energies) {

  int numEnergies = energies.size();
  long numStates = bandStructure.getNumStates();
  int numPoints = bandStructure.getNumPoints(true);
  int dim = context.getDimensionality();

  Eigen::Tensor<double, 3> energyProjVelocity(dim, dim, numEnergies);
  energyProjVelocity.setZero();

  TetrahedronDeltaFunction tetrahedra(bandStructure);

  if (mpi->mpiHead()) {
    std::cout << "Calculating energy projected velocity tensor" << std::endl;
  }

  for (long iEnergy = 0; iEnergy != numEnergies; ++iEnergy) {
    for (long iState = 0; iState != numStates; ++iState) {
      auto t = bandStructure.getIndex(iState);
      int ik = std::get<0>(t).get();
      int ib = std::get<1>(t).get();
      Eigen::Vector3d velocity = bandStructure.getGroupVelocity(iState);
      double deltaFunction = tetrahedra.getSmearing(energies(iEnergy), ik, ib);
      for (int j = 0; j < dim; ++j) {
        for (int i = 0; i < dim; ++i) {
          energyProjVelocity(i, j, iEnergy) +=
              velocity(i) * velocity(j) * deltaFunction / numPoints;
        }
      }
    }
  }

  return energyProjVelocity;
}

BaseVectorBTE TransportEpaApp::getScatteringRates(
    Context &context, StatisticsSweep &statisticsSweep,
    FullBandStructure &fullBandStructure, Eigen::VectorXd &energies) {

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
  if (hasSpinOrbit)
    spinFactor = 1;

  auto particle = Particle(Particle::electron);
  auto phParticle = Particle(Particle::phonon);

  if (particle.isPhonon())
    Error e("Electronic bandstructure has to be provided");

  long numCalcs = statisticsSweep.getNumCalcs();

  std::cout << "\nCalculate electronic density of states." << std::endl;
  TetrahedronDeltaFunction tetrahedra(fullBandStructure);

  long numEnergies = energies.size();
  double energyStep = context.getEnergyStep();

  // in principle, we should add 1 to account for ends of energy interval
  // i will not do that, because will work with the centers of energy steps
  //    long numEnergies = (long) (maxEnergy-minEnergy)/energyStep;

  // calculate the density of states at the energies in energies vector
  Eigen::VectorXd dos(numEnergies);
  for (long i = 0; i != numEnergies; ++i) {
    dos(i) = tetrahedra.getDOS(energies(i));
  }

  // get vector containing averaged phonon frequencies per mode
  InteractionEpa couplingEpa = InteractionEpa::parseEpaCoupling(context);

  Eigen::VectorXd phFreqAverage = couplingEpa.getPhFreqAverage();
  int numPhFreq = phFreqAverage.size();

  // phJump - contains the values of phFreqAverage/energyStep
  // defines in which step of electron energy grid the electron energy will go
  // after phonon absorption/emission
  Eigen::VectorXd phJump(phFreqAverage.size());
  for (auto i = 0; i != phFreqAverage.size(); ++i) {
    phJump(i) = phFreqAverage(i) / energyStep;
  }

  int numBandGroups = couplingEpa.getNumBandGroups();
  Eigen::VectorXd extrema = couplingEpa.getBandExtrema();
  Eigen::VectorXi numBins = couplingEpa.getNumBins();
  int numBinsMax = numBins.maxCoeff();
  Eigen::Tensor<double, 4> elPhMatElements = couplingEpa.getElPhMatAverage();
  Eigen::VectorXd binSize = couplingEpa.getBinSize();

  LoopPrint loopPrint("to calculate EPA scattering rates",
                      "pairs of temperatures and chemical potentials",
                      numCalcs);

  BaseVectorBTE epaRate(statisticsSweep, numEnergies, 1);

  // loop over temperatures and chemical potentials
  for (long iCalc = 0; iCalc < numCalcs; ++iCalc) {
    loopPrint.update();
    double temperature = statisticsSweep.getCalcStatistics(iCalc).temperature;
    double chemPotential =
        statisticsSweep.getCalcStatistics(iCalc).chemicalPotential;

    // loop over energies
    for (long iEnergy = 0; iEnergy < numEnergies; ++iEnergy) {

      // iTest: make sure that by phonon absorption or emission we will not go
      // outside of the considered energy range
      int iTest = (int)phJump(numPhFreq);
      if (iEnergy < iTest || iEnergy > numEnergies - iTest) {
        continue;
      }

      double scatRateTemp = 0.0;

      // loop over phonon frequencies
      for (int iPhFreq = 0; iPhFreq < numPhFreq; ++iPhFreq) {

        // population of phonons, electron after emission/absorption
        double nBose =
            phParticle.getPopulation(phFreqAverage(iPhFreq), temperature);
        double nFermiAbsorption =
            particle.getPopulation(energies[iEnergy] + phFreqAverage(iPhFreq),
                                   temperature, chemPotential);
        double nFermiEmission =
            particle.getPopulation(energies[iEnergy] - phFreqAverage(iPhFreq),
                                   temperature, chemPotential);

        int iJump = (int)phJump(iPhFreq);
        double iInterp = phJump(iPhFreq) - (double)iJump;
        double dosAbsorption = dos(iEnergy + iJump) * (1.0 - iInterp) +
                               dos(iEnergy + iJump + 1) * iInterp;
        double dosEmission = dos(iEnergy - iJump - 1) * iInterp +
                             dos(iEnergy - iJump) * (1.0 - iInterp);

        // this integer selects either valence or conduction bands
        int iBandGroup;
        if (energies(iEnergy) <= extrema.sum() / numBandGroups) {
          iBandGroup = 0;
        } else {
          iBandGroup = 1;
        }

        double iBinPos =
            (energies(iEnergy) - extrema(iBandGroup)) / binSize(iBandGroup);

        iBinPos = std::max(iBinPos, 1.0e-12);
        iBinPos = std::min(iBinPos, numBins(iBandGroup) - 1.0e-12);
        int intBinPos = (int)iBinPos;

        //------------------------------------
        // estimate strength of el-ph coupling

        double gAbsorption, gEmission;

        if (numBins(iBandGroup) == 1) {
          gAbsorption = elPhMatElements(iPhFreq, 0, 0, iBandGroup);
          gEmission = elPhMatElements(iPhFreq, 0, 0, iBandGroup);
        } else {
          Eigen::VectorXd elPhAvTemp(numBinsMax);
          for (int i = 0; i < numBinsMax; ++i) {
            elPhAvTemp(i) = elPhMatElements(iPhFreq, i, intBinPos, iBandGroup);
          }

          double iAbs = (energies(iEnergy) + phFreqAverage(iPhFreq) -
                         extrema(iBandGroup)) /
                        binSize(iBandGroup);
          iAbs = std::max(iAbs, 1.0e-12);
          iAbs = std::min(iAbs, numBins(iBandGroup) - 1.0e-12);
          int iAbsInt = (int)iAbs;

          gAbsorption = elPhAvTemp(iAbsInt);

          double iEmis = (energies[iEnergy] - phFreqAverage(iPhFreq) -
                          extrema(iBandGroup)) /
                         binSize(iBandGroup);
          iEmis = std::max(iEmis, 1.0e-12);
          iEmis = std::min(iEmis, numBins(iBandGroup) - 1.0e-12);
          int iEmisInt = (int)iEmis;

          gEmission = elPhAvTemp(iEmisInt);
        }

        //-----------------------------
        // finally, the scattering rate

        scatRateTemp +=
            gAbsorption * (nBose + nFermiAbsorption) * dosAbsorption +
            gEmission * (nBose + 1 - nFermiEmission) * dosEmission;
      }
      // scattering rate in Rydbergs
      epaRate.data(iCalc, iEnergy) = twoPi * scatRateTemp / spinFactor;
    }
  }
  loopPrint.close();
  return epaRate;
}
