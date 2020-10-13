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
          std::string bandsType, std::string outFileName,
          std::string energyUnit, double energyConversion,
          double chemicalPotential);

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
        "phonon_bands.json", "cmm-1", ryToCmm1, 0.0);

  if ( mpi->mpiHead()) {
    std::cout << "Finishing phonon bands calculation" << std::endl;
  }
/*
 *
  // Save phonon band structure to file
  std::vector<std::vector<double>> outEnergies; 
  std::vector<int> wavevectorIndices;
  std::vector<double> tempEns;  
  std::vector<std::vector<double>> pathCoords; 
  if ( mpi->mpiHead()) {
    int numBands = phononH0.getNumBands();
    for (long ik = 0; ik < pathPoints.getNumPoints(); ik++) {
      // store wavevector indices 
      wavevectorIndices.push_back(ik);
      auto ikIndex = WavevectorIndex(ik);
      Eigen::VectorXd energies = fullBandStructure.getEnergies(ikIndex);

      // store the path coordinates 
      auto coord = pathPoints.getPointCoords(ik); 
      pathCoords.push_back({coord[0],coord[1],coord[2]});

      // store the energies 
      Eigen::VectorXd energies = fullBandStructure.getEnergies(ikIndex);
      for (int ib = 0; ib < numBands; ib++) {
        tempEns.push_back(energies(ib)*ryToCmm1); 
      }
      outEnergies.push_back(tempEns);
      tempEns.clear(); 
    }
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
      Eigen::Vector3d tempCoords1 = {pathExtrema(pe,0,0),pathExtrema(pe,0,1),pathExtrema(pe,0,2)};
      pathLabelIndices.push_back(pathPoints.getIndex(tempCoords1));
      Eigen::Vector3d tempCoords2 = {pathExtrema(pe,1,0),pathExtrema(pe,1,1),pathExtrema(pe,1,2)};
      pathLabelIndices.push_back(pathPoints.getIndex(tempCoords2));
    }

    // output to json 
    nlohmann::json output;
    output["pathIndices"] = wavevectorIndices;
    output["pathCoordinates"] = pathCoords;
    output["highSymLabels"] = context.getPathLabels();
    output["highSymIndices"] = pathLabelIndices; 
    output["highSymCoordinates"] = extremaCoords; 
    output["numBands"] = numBands;
    output["energies"] = outEnergies;
    output["mu"] = 0.0;
    output["bandsType"] = "phonon";
    output["energyUnits"] = "cmm-1"; 
    output["coordsType"] = "cartesian";
    std::ofstream o("phonon_bands.json");
    o << std::setw(3) << output << std::endl;
*/
    // arguments: numBands, context, pathPoints, bandsType, outputFileName, 
    // energyUnits, energyConversionFactor, chemicalPotential 
   // outputBandsToJSON(numBands, context, pathPoints, "phonon", "phonon_bands.json", 
    //                    "cmm-1", ryToCmm1, 0.0);

   // std::cout << "Finishing phonon bands calculation" << std::endl;
  //}  
}

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

  // Save phonon band structure to file
  if ( mpi->mpiHead()) {
    int numBands = electronH0.getNumBands();
    std::ofstream outfile("./electron_bands.dat");
    outfile << "# Electron bands: path index, Bands[eV]" << std::endl;
    for (long ik = 0; ik < pathPoints.getNumPoints(); ik++) {
      outfile << ik;
      auto ikIndex = WavevectorIndex(ik);
      Eigen::VectorXd energies = fullBandStructure.getEnergies(ikIndex);
      for (int ib = 0; ib < numBands; ib++) {
        outfile << "\t" << energies(ib) * energyRyToEv;
      }
      outfile << std::endl;
    }
    std::cout << "Finishing electron (Wannier) bands calculation" << std::endl;
  }
}

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
  Eigen::VectorXd dummyZero(1);
  dummyZero(0) = 0.0; // set both temperature and doping to zero
  context.setTemperatures(dummyZero); 
  context.setDopings(dummyZero); 
  StatisticsSweep statisticsSweep(context,&fullBandStructure);
  // in this case, there should only be one object containing calculated 
  // statistics to reference at index 0
  auto stats = statisticsSweep.getCalcStatistics(0);

  // arguments: bandStructure, context, pathPoints, bandsType, outputFileName, 
  // energyUnits, energyConversionFactor, chemicalPotential 
  outputBandsToJSON(fullBandStructure, context, pathPoints, "electron",
        "electron_bands.json", "eV", energyRyToEv, stats.chemicalPotential);

  if ( mpi->mpiHead()) {
    std::cout << "Finishing electron (Fourier) bands calculation" << std::endl;
  }

  // Save phonon band structure to file
/*  if ( mpi->mpiHead()) {
    int numBands = electronH0.getNumBands();
    std::ofstream outfile("./electron_bands.dat");
    outfile << "# Electron bands: path index, Bands[eV]" << std::endl;
    for (long ik = 0; ik < pathPoints.getNumPoints(); ik++) {
      outfile << ik;
      auto ikIndex = WavevectorIndex(ik);
      Eigen::VectorXd energies = fullBandStructure.getEnergies(ikIndex);
      for (int ib = 0; ib < numBands; ib++) {
        outfile << "\t" << energies(ib) * energyRyToEv;
      }
      outfile << std::endl;
    }*/
}
/* helper function to output bands to a json file */ 
void outputBandsToJSON(FullBandStructure& fullBandStructure, 
                 Context& context, PathPoints& pathPoints,
                 std::string bandsType, std::string outFileName,
                 std::string energyUnit, double energyConversion,
                 double chemicalPotential) {

  if ( mpi->mpiHead()) {
    std::vector<std::vector<double>> outEnergies;
    std::vector<int> wavevectorIndices;
    std::vector<double> tempEns;
    std::vector<std::vector<double>> pathCoords;
    int numBands = fullBandStructure.getNumBands();

    for (long ik = 0; ik < pathPoints.getNumPoints(); ik++) {
      // store wavevector indices 
      wavevectorIndices.push_back(ik);
      auto ikIndex = WavevectorIndex(ik);
      
      // store the path coordinates 
      auto coord = pathPoints.getPointCoords(ik);
      pathCoords.push_back({coord[0],coord[1],coord[2]});
      
      // store the energies 
      Eigen::VectorXd energies = fullBandStructure.getEnergies(ikIndex);
      for (int ib = 0; ib < numBands; ib++) {
        tempEns.push_back(energies(ib)*energyConversion);
      }
      outEnergies.push_back(tempEns);
      tempEns.clear();
    }

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
      Eigen::Vector3d tempCoords1 = {pathExtrema(pe,0,0),pathExtrema(pe,0,1),pathExtrema(pe,0,2)};
      pathLabelIndices.push_back(pathPoints.getIndex(tempCoords1));
      Eigen::Vector3d tempCoords2 = {pathExtrema(pe,1,0),pathExtrema(pe,1,1),pathExtrema(pe,1,2)};
      pathLabelIndices.push_back(pathPoints.getIndex(tempCoords2));
    }

    // output to json 
    nlohmann::json output;
    output["pathIndices"] = wavevectorIndices;
    output["pathCoordinates"] = pathCoords;
    output["highSymLabels"] = context.getPathLabels();
    output["highSymIndices"] = pathLabelIndices;
    output["highSymCoordinates"] = extremaCoords;
    output["numBands"] = numBands; 
    output["energies"] = outEnergies;
    output["chemicalPotential"] = chemicalPotential*energyConversion;
    output["bandsType"] = bandsType;
    output["energyUnits"] = energyUnit;
    output["coordsType"] = "lattice";
    std::ofstream o(outFileName);
    o << std::setw(3) << output << std::endl;
    o.close();
  }
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
