#include "bands_app.h"
#include <iostream>
#include <fstream>
#include "eigen.h"
#include "constants.h"
#include "qe_input_parser.h"
#include "path_points.h"
#include "mpiHelper.h"
#include <nlohmann/json.hpp>

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

  // Save phonon band structure to file
  std::vector<std::vector<double>> outEnergies; 
  std::vector<int> wavevectorIndices;
  std::vector<double> tempEns;  
  std::vector<std::vector<double>> pathCoords; 
  if ( mpi->mpiHead()) {
    int numBands = phononH0.getNumBands();
    std::ofstream outfile("./phonon_bands.dat");
    outfile << "# Phonon bands: path index, Bands[cmm1]" << std::endl;
    for (long ik = 0; ik < pathPoints.getNumPoints(); ik++) {
      outfile << ik;
      wavevectorIndices.push_back(ik);
      auto ikIndex = WavevectorIndex(ik);
      Eigen::VectorXd energies = fullBandStructure.getEnergies(ikIndex);
      auto coord = pathPoints.getPointCoords(ik); 
      pathCoords.push_back({coord[0],coord[1],coord[2]});
      for (int ib = 0; ib < numBands; ib++) {
        outfile << "\t" << energies(ib) * ryToCmm1;
        tempEns.push_back(energies(ib)); 
      }
      outEnergies.push_back(tempEns);
      tempEns.clear(); 
      outfile << std::endl;
    }
    // determine path extrema to output to json
    std::vector<std::vector<double>> outExtrema; 
    std::vector<std::string> pathLabels = context.getPathLabels(); 
    Eigen::Tensor<double, 3> pathExtrema = context.getPathExtrema(); 
    auto numExtrema = pathExtrema.dimensions(); 
    for (int pe = 0; pe < numExtrema[0]; pe++) {
      outExtrema.push_back({pathExtrema(pe,0,0),pathExtrema(pe,0,1),pathExtrema(pe,0,2)});
    }

    // output to json 
    nlohmann::json output;
    output["wavevectorIndices"] = wavevectorIndices; 
    output["pathCoordinates"] = pathCoords; 
    output["pathLabels"] = pathLabels; 
    output["numBands"] = numBands; 
    output["energies"] = outEnergies;
    output["mu"] = 0.0; 
    output["bandsType"] = "phonon";
    std::ofstream o("phonon.json");
    o << output << std::endl;

    std::cout << "Finishing phonon bands calculation" << std::endl;
  }  
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
    std::cout << "Finishing electron (Fourier) bands calculation" << std::endl;
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
