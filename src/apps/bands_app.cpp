#include "bands_app.h"
#include <iostream>
#include <fstream>
#include "eigen.h"
#include "constants.h"
#include "qe_input_parser.h"
#include "path_points.h"

void PhononBandsApp::run(Context &context) {
    std::cout << "Starting phonon bands calculation" << std::endl;

    // Read the necessary input files
    auto tup = QEParser::parsePhHarmonic(context);
 auto crystal = std::get<0>(tup);
 auto phononH0 = std::get<1>(tup);

    // first we make compute the band structure on the fine grid
    PathPoints pathPoints(crystal, context.getPathExtrema(),
            context.getDeltaPath());

    bool withVelocities = false;
    bool withEigenvectors = false;
    FullBandStructure fullBandStructure = phononH0.populate(pathPoints,
            withVelocities, withEigenvectors);

    // Save phonon band structure to file
    long numPoints = pathPoints.getNumPoints();
    long numBands = phononH0.getNumBands();

    std::ofstream outfile("./phonon_bands.dat");
    outfile << "# Phonon bands: path index, Bands[cmm1]" << std::endl;

    for (long ik = 0; ik < numPoints; ik++) {
        outfile << ik;

        auto p = fullBandStructure.getPoint(ik);
        auto s = fullBandStructure.getState(p);
        Eigen::VectorXd energies = s.getEnergies();
        for (long ib = 0; ib < numBands; ib++) {
            outfile << "\t" << energies(ib) * ryToCmm1;
        }
        outfile << std::endl;
    }
    std::cout << "Finishing phonon bands calculation" << std::endl;
}

void ElectronWannierBandsApp::run(Context &context) {
    std::cout << "Starting electron (Wannier) bands calculation" << std::endl;

    // Read the necessary input files
    auto tup = QEParser::parseElHarmonicWannier(context);
 auto crystal = std::get<0>(tup);
 auto electronH0 = std::get<1>(tup);

    // first we make compute the band structure on the fine grid
    PathPoints pathPoints(crystal, context.getPathExtrema(),
            context.getDeltaPath());

    bool withVelocities = false;
    bool withEigenvectors = false;
    FullBandStructure fullBandStructure = electronH0.populate(pathPoints,
            withVelocities, withEigenvectors);

    // Save phonon band structure to file
    long numPoints = pathPoints.getNumPoints();
    long numBands = electronH0.getNumBands();

    std::ofstream outfile("./electron_bands.dat");
    outfile << "# Electron bands: path index, Bands[eV]" << std::endl;

    for (long ik = 0; ik < numPoints; ik++) {
        outfile << ik;
        auto p = fullBandStructure.getPoint(ik);
        auto s = fullBandStructure.getState(p);
        Eigen::VectorXd energies = s.getEnergies();
        for (long ib = 0; ib < numBands; ib++) {
            outfile << "\t" << energies(ib) * energyRyToEv;
        }
        outfile << std::endl;
    }
    std::cout << "Finishing electron (Wannier) bands calculation" << std::endl;
}

void ElectronFourierBandsApp::run(Context &context) {
    std::cout << "Starting electron (Fourier) bands calculation" << std::endl;

    // Read the necessary input files
    auto tup = QEParser::parseElHarmonicFourier(context);
 auto crystal = std::get<0>(tup);
 auto electronH0 = std::get<1>(tup);

    // first we make compute the band structure on the fine grid
    PathPoints pathPoints(crystal, context.getPathExtrema(),
            context.getDeltaPath());

    bool withVelocities = false;
    bool withEigenvectors = false;
    FullBandStructure fullBandStructure = electronH0.populate(pathPoints,
            withVelocities, withEigenvectors);

    // Save phonon band structure to file
    long numPoints = pathPoints.getNumPoints();
    long numBands = electronH0.getNumBands();

    std::ofstream outfile("./electron_bands.dat");
    outfile << "# Electron bands: path index, Bands[eV]" << std::endl;

    for (long ik = 0; ik < numPoints; ik++) {
        outfile << ik;
        auto p = fullBandStructure.getPoint(ik);
        auto s = fullBandStructure.getState(p);
        Eigen::VectorXd energies = s.getEnergies();
        for (long ib = 0; ib < numBands; ib++) {
            outfile << "\t" << energies(ib) * energyRyToEv;
        }
        outfile << std::endl;
    }
    std::cout << "Finishing electron (Fourier) bands calculation" << std::endl;
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

