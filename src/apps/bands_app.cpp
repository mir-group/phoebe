#include "bands_app.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "eigen.h"
#include "constants.h"
#include "qe_input_parser.h"
#include "path_points.h"
#include "mpiHelper.h"
#include <nlohmann/json.hpp>

// forward declare this helper function, so we can leave the important
// run functions at the top
void outputBandsToJSON(FullBandStructure& fullBandStructure,
          Context& context, PathPoints& pathPoints,
          std::string particleType, std::string outFileName,
          std::string energyUnit, double energyConversion,
          double chemicalPotential);

/* --------------------- PhononBandsApp ------------------------------ */
void PhononBandsApp::run(Context &context) {

  if ( mpi->mpiHead()) {
    std::cout << "Starting phonon bands calculation" << std::endl;
  }

  // Read the necessary input files
  auto tup = QEParser::parsePhHarmonic(context);
  auto crystal = std::get<0>(tup);
  auto phononH0 = std::get<1>(tup);

  // first we make compute the band structure on the fine grid
  PathPoints pathPoints(crystal, context.getPathExtrema(),
                        context.getDeltaPath());
  bool withVelocities = false;
  bool withEigenvectors = false;
  FullBandStructure fullBandStructure =
      phononH0.populate(pathPoints, withVelocities, withEigenvectors);

  // arguments: bandStructure, context, pathPoints, bandsType, outputFileName,
  // energyUnits, energyConversionFactor, chemicalPotential
  outputBandsToJSON(fullBandStructure, context, pathPoints, "phonon",
        "phonon_bands.json", "cm$^{-1}$", ryToCmm1, 0.0);

  if ( mpi->mpiHead()) {
    std::cout << "Finishing phonon bands calculation" << std::endl;
  }
}

/* --------------------- ElectronWannierBandsApp --------------------- */
void ElectronWannierBandsApp::run(Context &context) {
  if ( mpi->mpiHead()) {
    std::cout << "Starting electron (Wannier) bands calculation" << std::endl;
  }

  // Read the necessary input files
  auto tup = QEParser::parseElHarmonicWannier(context);
  auto crystal = std::get<0>(tup);
  auto electronH0 = std::get<1>(tup);

  // first we make compute the band structure on the fine grid
  PathPoints pathPoints(crystal, context.getPathExtrema(),
                        context.getDeltaPath());

  bool withVelocities = false;
  bool withEigenvectors = false;
  FullBandStructure fullBandStructure =
      electronH0.populate(pathPoints, withVelocities, withEigenvectors);

  // Use statisticsSweep to get the chemical potential
  // TODO do we want to do this, or would we prefer to use whatever was read in
  // by context (the value provided by QE)
  Eigen::VectorXd dummyZero(1);
  dummyZero(0) = 0.0; // set both temperature and doping to zero
  context.setTemperatures(dummyZero);
  context.setDopings(dummyZero);
  StatisticsSweep statisticsSweep(context,&fullBandStructure);
  auto stats = statisticsSweep.getCalcStatistics(0);

  // arguments: bandStructure, context, pathPoints, bandsType, outputFileName,
  // energyUnits, energyConversionFactor, chemicalPotential
  outputBandsToJSON(fullBandStructure, context, pathPoints, "electron",
        "electron_bands.json", "eV", energyRyToEv, stats.chemicalPotential);

  if ( mpi->mpiHead()) {
    std::cout << "Finishing electron (Wannier) bands calculation" << std::endl;
  }
}

/* --------------------- ElectronFourierBandsApp --------------------- */
void ElectronFourierBandsApp::run(Context &context) {
  if ( mpi->mpiHead()) {
    std::cout << "Starting electron (Fourier) bands calculation" << std::endl;
  }

  // Read the necessary input files
  auto tup = QEParser::parseElHarmonicFourier(context);
  auto crystal = std::get<0>(tup);
  auto electronH0 = std::get<1>(tup);

  // first we make compute the band structure on the fine grid
  PathPoints pathPoints(crystal, context.getPathExtrema(),
                        context.getDeltaPath());

  bool withVelocities = false;
  bool withEigenvectors = false;
  FullBandStructure fullBandStructure =
      electronH0.populate(pathPoints, withVelocities, withEigenvectors);

  // Use statisticsSweep to get the chemical potential
  // TODO do we want to do this, or would we prefer to use whatever was read in
  // by context (the value provided by QE)
  Eigen::VectorXd dummyZero(1);
  dummyZero(0) = 0.0; // set both temperature and doping to zero
  context.setTemperatures(dummyZero);
  context.setDopings(dummyZero);
  StatisticsSweep statisticsSweep(context,&fullBandStructure);
  auto stats = statisticsSweep.getCalcStatistics(0);

  // arguments: bandStructure, context, pathPoints, particleType, outputFileName,
  // energyUnits, energyConversionFactor, chemicalPotential
  outputBandsToJSON(fullBandStructure, context, pathPoints, "electron",
        "electron_bands.json", "eV", energyRyToEv, stats.chemicalPotential);

  if ( mpi->mpiHead()) {
    std::cout << "Finishing electron (Fourier) bands calculation" << std::endl;
  }
}

/* helper function to output bands to a json file */
void outputBandsToJSON(FullBandStructure& fullBandStructure,
                 Context& context, PathPoints& pathPoints,
                 std::string particleType, std::string outFileName,
                 std::string energyUnit, double energyConversion,
                 double chemicalPotential) {

  if (!mpi->mpiHead()) return;

  std::vector<std::vector<double>> outEnergies;
  std::vector<int> wavevectorIndices;
  std::vector<double> tempEns;
  std::vector<std::vector<double>> pathCoords;
  int numBands = fullBandStructure.getNumBands();

  // determine path extrema and their ik indices, to output to json
  std::vector<std::vector<double>> extremaCoords;
  std::vector<int> pathLabelIndices;
  Eigen::Tensor<double, 3> pathExtrema = context.getPathExtrema();
  auto numExtrema = pathExtrema.dimensions();
  for (int pe = 0; pe < numExtrema[0]; pe++) {
    // store coordinates of the extrema
    extremaCoords.push_back({pathExtrema(pe,0,0),pathExtrema(pe,0,1),pathExtrema(pe,0,2)});
    extremaCoords.push_back({pathExtrema(pe,1,0),pathExtrema(pe,1,1),pathExtrema(pe,1,2)});
    // determine the indices of the extrema and save for plotting
    //Eigen::Vector3d tempCoords1 = {pathExtrema(pe,0,0),pathExtrema(pe,0,1),pathExtrema(pe,0,2)};
    //pathLabelIndices.push_back(pathPoints.getIndex(tempCoords1));
    //Eigen::Vector3d tempCoords2 = {pathExtrema(pe,1,0),pathExtrema(pe,1,1),pathExtrema(pe,1,2)};
    //pathLabelIndices.push_back(pathPoints.getIndex(tempCoords2));
  }

  unsigned int extremaCount = 0; //keeps track of how many high sym points we've found
  for (long ik = 0; ik < pathPoints.getNumPoints(); ik++) {
    // store wavevector indices
    wavevectorIndices.push_back(ik);
    auto ikIndex = WavevectorIndex(ik);

    // store the path coordinates
    auto coord = pathPoints.getPointCoords(ik);
    pathCoords.push_back({coord[0],coord[1],coord[2]});

    // check if this point is one of the high sym points,
    // and if it is, save the index
    if(coord[0] == extremaCoords[extremaCount][0] &&
          coord[1] == extremaCoords[extremaCount][1] &&
          coord[2] == extremaCoords[extremaCount][2]) {
      pathLabelIndices.push_back(ik);
      // if the next extrema is the same as this one, add the index
      //  twice and increment extrema count twice.
      // we have to do this because the point won't occur a second time
      // in the path list.
      if(extremaCount+1 < extremaCoords.size()) { //check bounds
        if(extremaCoords[extremaCount] == extremaCoords[extremaCount+1])  {
          pathLabelIndices.push_back(ik);
          extremaCount += 2;
        } else {
          extremaCount++;
        }
      }
    }

    // store the energies
    Eigen::VectorXd energies = fullBandStructure.getEnergies(ikIndex);
    for (int ib = 0; ib < numBands; ib++) {
      tempEns.push_back(energies(ib) * energyConversion);
    }
    outEnergies.push_back(tempEns);
    tempEns.clear();
  }

  // output to json
  nlohmann::json output;
  output["wavevectorIndices"] = wavevectorIndices;
  output["wavevectorCoordinates"] = pathCoords;
  output["highSymLabels"] = context.getPathLabels();
  output["highSymIndices"] = pathLabelIndices;
  output["highSymCoordinates"] = extremaCoords;
  output["numBands"] = numBands;
  output["energies"] = outEnergies;
  output["chemicalPotential"] = chemicalPotential*energyConversion;
  output["particleType"] = particleType;
  output["energyUnit"] = energyUnit;
  output["coordsType"] = "lattice";
  std::ofstream o(outFileName);
  o << std::setw(3) << output << std::endl;
  o.close();
}

void PhononBandsApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getPhD2FileName(), "PhD2FileName");
  throwErrorIfUnset(context.getPathExtrema(), "points path extrema");
  throwErrorIfUnset(context.getDeltaPath(), "deltaPath");
  throwWarningIfUnset(context.getSumRuleD2(), "sumRuleD2");
}

void ElectronWannierBandsApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");
  throwErrorIfUnset(context.getPathExtrema(), "points path extrema");
  throwErrorIfUnset(context.getDeltaPath(), "deltaPath");

  std::string crystalMsg = "crystal structure";
  throwErrorIfUnset(context.getInputAtomicPositions(), crystalMsg);
  throwErrorIfUnset(context.getInputSpeciesNames(), crystalMsg);
  throwErrorIfUnset(context.getInputAtomicSpecies(), crystalMsg);
}

void ElectronFourierBandsApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");
  throwErrorIfUnset(context.getPathExtrema(), "points path extrema");
  throwErrorIfUnset(context.getDeltaPath(), "deltaPath");
  throwErrorIfUnset(context.getElectronFourierCutoff(),
                    "electronFourierCutoff");
}
