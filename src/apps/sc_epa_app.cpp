#include "sc_epa_app.h"
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
#include <nlohmann/json.hpp>
#include <vector>

void ScEpaApp::run(Context &context) {

  // parse QE-xml file
  auto t1 = QEParser::parseElHarmonicFourier(context);
  auto crystal = std::get<0>(t1);
  auto electronH0 = std::get<1>(t1);

  //--------------------------------
  // Setup energy grid
  // TODO run with a tiny energy range and 1 bin
  double minEnergy = context.getFermiLevel() - context.getEpaEnergyRange();
  double maxEnergy = context.getFermiLevel() + context.getEpaEnergyRange();
  double energyStep = context.getEpaEnergyStep();

  // in principle, we should add 1 to account for ends of energy interval
  // I will not do that, because will work with the centers of energy steps
  int numEnergies = int((maxEnergy - minEnergy) / energyStep);
  // energies at the centers of energy steps
  Eigen::VectorXd energies(numEnergies);
  for (int i = 0; i < numEnergies; ++i) {
    // add 0.5 to be in the middle of the energy step
    energies(i) = (double(i) + 0.5) * energyStep + minEnergy;
  }

  if(mpi->mpiHead()) std::cout << "numEnergies " << numEnergies << std::endl;

  // Read and setup k-point mesh for interpolating band structure
  Points fullPoints(crystal, context.getKMesh());
  fullPoints.setIrreduciblePoints();

  if (mpi->mpiHead()) {
    std::cout << "\nBuilding electronic band structure" << std::endl;
  }

  // filter to only the bands relevant to transport
  electronH0.trimBands(context, minEnergy, maxEnergy);

  // Fourier's interpolation of the electronic band structure ----------------
  bool withVelocities = false; // TODO in this case we don't need velocities (for now)
  bool withEigenvectors = false;
  FullBandStructure bandStructure =
      electronH0.populate(fullPoints, withVelocities, withEigenvectors);

  if(mpi->mpiHead()) std::cout << "numBands post trim " << bandStructure.getNumBands() << std::endl;

  // set temperatures, chemical potentials and carrier concentrations
  StatisticsSweep statisticsSweep(context, &bandStructure);

  // TODO for now let's start with 1 calculation
  if(statisticsSweep.getNumCalcs() != 1) Error("scEPA is only for 1 calc at a time right now!");
  double efermi = statisticsSweep.getCalcStatistics().chemicalPotential;
  TetrahedronDeltaFunction tetrahedrons(bandStructure);
  double NF = tetrahedrons.getDOS(efermi);

  // TODO throw an error if mu isn't inside the slim energy range?

  // interpolate the phonon band structure here ----------------------
  auto tup = Parser::parsePhHarmonic(context);
  auto crystal = std::get<0>(tup);
  auto phononH0 = std::get<1>(tup);

  // throw an error if k/q points arent the same somehow

  if (mpi->mpiHead()) {
    std::cout << "\nComputing phonon band structure." << std::endl;
  }
  auto tup1 = ActiveBandStructure::builder(context, phononH0, fullPoints);
  auto phBandStructure = std::get<0>(tup1);
  int numModes = phononH0.getNumBands();

  // print some info about state number reduction
  if (mpi->mpiHead()) {
    if(bandStructure.hasWindow() != 0) {
        std::cout << "Window selection reduced phonon band structure from "
                << fullPoints.getNumPoints()*phononH0.getNumBands() << " to "
                << phBandStructure.getNumStates() << " states."  << std::endl;
    }
    if(context.getUseSymmetries()) {
      std::cout << "Symmetries reduced phonon band structure from "
          << phBandStructure.getNumStates() << " to "
          << phBandStructure.irrStateIterator().size() << " states." << std::endl;
    }
    std::cout << "Done computing phonon band structure.\n" << std::endl;
  }

  // set up a delta function for integration
  GaussianDeltaFunction delta(context);

  // TODO we need to set this to a correct factor
  double norm = twoPi / spinFactor * crystal.getVolumeUnitCell(crystal.getDimensionality());

  // get vector containing averaged phonon frequencies per mode
  InteractionEpa couplingEpa = InteractionEpa::parseEpaCoupling(context);

  LoopPrint loopPrint("calculation of lambda values", "phonon modes",

  Eigen::MatrixXd lambda(numModes, numPoints);

  // TODO add parallel

  // loop over phonon modes
  for (int iMode = 0; iMode < numModes; iModes++) {

    // get mode-dept g coupling TODO make sure the first index is mode, and that we can trick it with Ef
    // TODO check that this is sq
    double gvSq = couplingEpa.getCoupling(iMode, efermi, efermi);

    #pragma OMP critical
    { loopPrint.update(); }

    // loop over qpoints
    for (int iq = 0; iq < numPoints; iq++) {

      // phonon energy
      double omega = phBandStructure.getEnergy(WavectorIndex(iq), BandIndex(iMode));

      // joint density of states integral
      double jdos = 0;

      // do the phase space integral here ------------------
      for (size_t ik = 0; ik < numPoints; ik++) {

        // get energies at k
        WavevectorIndex kIdx = WavevectorIndex(ik);
        Eigen::VectorXd ek = bandStructure.getEnergies(kIdx);
        // k+q in cartesian, look up k+q idx, get energies at k+q
        Eigen::VectorXd kpq =
                bandStructure.getWavevector(kIdx) + phBandStructure.getWavevector(WavevectorIndex(iq));
        // TODO might need to fold this wavevector
        kpq = fullPoints.cartesianToCrystal(kpq);
        int ikq = fullPoints.getIndex(kpq);
        kqIdx = WavevectorIndex(ikq);
        Eigen::VectorXd ekpq = bandStructure.getEnergies(kIdx);

        // loop over band index
        for(int n = 0; n < ek.size(); n++) {
          for (int m = 0; m < ekpq.size(); m++) {
            jdos += delta(ek) + delta(ekpq);
          }
        }
      } // close phase space integral

      // normalize here -- TODO need also NF
      if(omega > 1e-8) { // is there a better way to handle this zero gamma info?
        lambda(iMode,iq) = 1./omega * gvSq * jdos/numPoints;
      }

    } // close point loop
  } // close mode
  loopPrint.close();


}

void ScEpaApp::checkRequirements(Context &context) {
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
    Error("Either chemical potentials or doping must be set");
  }
}

std::tuple<Eigen::Tensor<double, 3>, Eigen::VectorXd>
ScEpaApp::calcEnergyProjVelocity(Context &context,
                                        FullBandStructure &bandStructure,
                                        const Eigen::VectorXd &energies) {

  //--------------------------------
  // set up tetrahedron method
  TetrahedronDeltaFunction tetrahedrons(bandStructure);

  int numEnergies = int(energies.size());
  int numPoints = std::get<0>(bandStructure.getPoints().getMesh()).prod();

  if (mpi->mpiHead()) {
    std::cout << "\nCalculating energy projected velocity tensor" << std::endl;
  }

  auto crystal = bandStructure.getPoints().getCrystal();
  int dim = context.getDimensionality();
  double norm = pow(twoPi, dim) / crystal.getVolumeUnitCell(dim) / numPoints;
  double normDos = 1. / crystal.getVolumeUnitCell(dim) / numPoints;

  Eigen::Tensor<double, 3> energyProjVelocity(3, 3, numEnergies);
  energyProjVelocity.setZero();
  Eigen::VectorXd dos(numEnergies);
  dos.setZero();

  LoopPrint loopPrint("calculating energy projected velocity", "states",
      numEnergies);
  std::vector<size_t> iEnergies = mpi->divideWorkIter(numEnergies);
  size_t niEnergies = iEnergies.size();
#pragma omp parallel for
  for (size_t iiEnergy = 0; iiEnergy < niEnergies; iiEnergy++) {
    int iEnergy = iEnergies[iiEnergy];
#pragma omp critical
    { loopPrint.update(); }
    for (int iState : bandStructure.irrStateIterator()) {
      StateIndex isIdx(iState);
      auto rotations = bandStructure.getRotationsStar(isIdx);
      Eigen::Vector3d velIrr = bandStructure.getGroupVelocity(isIdx);
      double deltaFunction = tetrahedrons.getSmearing(energies(iEnergy), isIdx);

      // integrate DOS
      auto degeneracy = double(rotations.size());
      dos(iEnergy) += deltaFunction * degeneracy * normDos;

      for (const Eigen::Matrix3d &r : rotations) {
        Eigen::Vector3d velocity = r * velIrr;
        for (int j  = 0; j < 3; j++) {
          for (int i = 0; i < 3; i++) {
            energyProjVelocity(i, j, iEnergy) +=
              velocity(i) * velocity(j) * deltaFunction * norm;
          }
        }
      }
    }
  }
  mpi->allReduceSum(&energyProjVelocity);
  mpi->allReduceSum(&dos);
  loopPrint.close();
  return std::make_tuple(energyProjVelocity, dos);
}

VectorEPA ScEpaApp::getScatteringRates(
    Context &context, StatisticsSweep &statisticsSweep,
    const Eigen::VectorXd &energies, Crystal &crystal,
    const Eigen::VectorXd &dos) {

  int numEnergies = int(energies.size());

  /*If constant relaxation time is specified in input, we don't need to
  calculate EPA lifetimes*/
  double constantRelaxationTime = context.getConstantRelaxationTime();
  if (constantRelaxationTime > 0.) {
    VectorEPA crtRate(statisticsSweep, numEnergies, 1);
    crtRate.setConst(twoPi / constantRelaxationTime);
    return crtRate;
  }

  int spinFactor = 2;
  if (context.getHasSpinOrbit()) {
    spinFactor = 1;
  }

  auto particle = Particle(Particle::electron);

  int numCalculations = statisticsSweep.getNumCalculations();
  double energyStep = context.getEpaEnergyStep();

  // get vector containing averaged phonon frequencies per mode
  InteractionEpa couplingEpa = InteractionEpa::parseEpaCoupling(context);

  Eigen::VectorXd phEnergies = couplingEpa.getPhEnergies();
  int numModes = int(phEnergies.size());

  // precompute bose-einstein populations
  Eigen::MatrixXd bose(numModes, numCalculations);
  {
    auto phParticle = Particle(Particle::phonon);
    for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
      double temp = statisticsSweep.getCalcStatistics(iCalc).temperature;
      for (int iPhFreq = 0; iPhFreq < numModes; iPhFreq++) {
        bose(iPhFreq, iCalc) =
            phParticle.getPopulation(phEnergies(iPhFreq), temp);
      }
    }
  }

  LoopPrint loopPrint("calculation of EPA scattering rates", "phonon modes",
                      int(mpi->divideWorkIter(numEnergies).size()));

  VectorEPA epaRate(statisticsSweep, numEnergies, 1);

  double norm = twoPi / spinFactor *
                crystal.getVolumeUnitCell(crystal.getDimensionality());

  std::vector<size_t> iEnergies = mpi->divideWorkIter(numEnergies);
  size_t niEnergies = iEnergies.size();

#pragma omp parallel for
  for (size_t iiEnergy = 0; iiEnergy < niEnergies; iiEnergy++) {
    int iEnergy = iEnergies[iiEnergy];

#pragma omp critical
    { loopPrint.update(); }

    // loop over phonon frequencies
    for (int iPhFreq = 0; iPhFreq < numModes; iPhFreq++) {

      // note: phEnergies(iPhFreq)/energyStep =
      // # of bin-jumps the electron does after scattering

      // compute the dos for electron in the final state for the two
      // scatterings mechanisms, and do a linear interpolation
      int iJump = int(phEnergies(iPhFreq) / energyStep);
      double iInterp = phEnergies(iPhFreq) / energyStep - double(iJump);
      auto largeIndex = int(iEnergy + iJump + 1);
      auto smallIndex = int(iEnergy - iJump - 1);
      double dosEmission = 0.;
      double dosAbsorption = 0.;
      // Avoid some index out of bound errors
      if (smallIndex >= 0) {
        dosEmission =
            dos(smallIndex) * iInterp + dos(iEnergy - iJump) * (1. - iInterp);
      }
      if (largeIndex < dos.size()) {
        dosAbsorption =
            dos(iEnergy + iJump) * (1. - iInterp) + dos(largeIndex) * iInterp;
      }

      //------------------------------------
      // estimate strength of el-ph coupling |g|^2
      double en = energies(iEnergy);
      double enP = energies(iEnergy) + phEnergies(iPhFreq);
      double enM = energies(iEnergy) - phEnergies(iPhFreq);
      double gAbsorption = couplingEpa.getCoupling(iPhFreq, en, enP);
      double gEmission = couplingEpa.getCoupling(iPhFreq, en, enM);

      // population of phonons, electron after emission/absorption
      for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
        auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
        double temp = calcStat.temperature;
        double chemPot = calcStat.chemicalPotential;

        double nFermiAbsorption = particle.getPopulation(
            energies[iEnergy] + phEnergies(iPhFreq), temp, chemPot);
        double nFermiEmission = particle.getPopulation(
            energies[iEnergy] - phEnergies(iPhFreq), temp, chemPot);

        //-----------------------------
        // now the scattering rate
        // Note that g (the coupling) is already squared

        double rate = gAbsorption * (bose(iPhFreq, iCalc) + nFermiAbsorption) *
                          dosAbsorption +
                      gEmission * (bose(iPhFreq, iCalc) + 1 - nFermiEmission) *
                          dosEmission;

        epaRate.data(iCalc, iEnergy) += norm * rate;
      }
    }
  }
  mpi->allReduceSum(&epaRate.data);
  loopPrint.close();

  return epaRate;
}

/* helper function to output scattering rates at each energy to JSON */
void ScEpaApp::outputToJSON(const std::string &outFileName,
                                   VectorEPA &scatteringRates,
                                   StatisticsSweep &statisticsSweep,
                                   Eigen::VectorXd &energiesEPA, Context &context) {

  if (!mpi->mpiHead()) {
    return;
  }

  std::string particleType = "electron";
  double energyConversion = energyRyToEv;
  std::string energyUnit = "eV";
  double energyToTime = energyRyToFs;

  int dimensionality = context.getDimensionality();

  // need to store as a vector format with dimensions
  // iCalc, ik. ib, iDim (where iState is unfolded into
  // ik, ib) for the velocities and lifetimes, no dim for energies
  std::vector<std::vector<double>> outTimes;
  std::vector<std::vector<double>> outLinewidths;
  std::vector<std::vector<double>> energies;
  std::vector<double> temps;
  std::vector<double> chemPots;
  std::vector<double> dopings;

  auto numEnergies = int(energiesEPA.size());

  for (int iCalc = 0; iCalc < statisticsSweep.getNumCalculations(); iCalc++) {
    auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStatistics.temperature;
    double chemPot = calcStatistics.chemicalPotential;
    double doping = calcStatistics.doping;
    temps.push_back(temp * temperatureAuToSi);
    chemPots.push_back(chemPot * energyConversion);
    dopings.push_back(doping);

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
  output["chemicalPotentialUnit"] = "eV";
  output["dopingConcentrations"] = dopings;
  output["dopingConcentrationUnit"] = "cm$^{-" + std::to_string(dimensionality) + "}$";
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
