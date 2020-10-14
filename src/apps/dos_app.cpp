#include "dos_app.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

#include "constants.h"
#include "delta_function.h"
#include "electron_h0_fourier.h"
#include "exceptions.h"
#include "full_points.h"
#include "mpiHelper.h"
#include "qe_input_parser.h"
#include "utilities.h"
#include <nlohmann/json.hpp>

// Forward declare this helper function so that 
// app::run functions can remain at the top 
std::tuple<std::vector<double>, std::vector<double>> calcDOS(
      Context& context, FullBandStructure& fullBandStructure);

void outputDOSToJSON(std::vector<double> energies, std::vector<double> dos,
                 std::string particleType, std::string outFileName,
                 std::string energyUnit, std::string dosUnit,  
                 double energyConversion, double chemicalPotential);

/* ------------------- PhononDoSApp --------------------*/
// Compute the DOS with the tetrahedron method
void PhononDosApp::run(Context &context) {
  if (mpi->mpiHead()) {
    std::cout << "Starting phonon DoS calculation" << std::endl;
  }

  // Read the necessary input files
  auto tup = QEParser::parsePhHarmonic(context);
  auto crystal = std::get<0>(tup);
  auto phononH0 = std::get<1>(tup);

  // first we make compute the band structure on the fine grid
  FullPoints fullPoints(crystal, context.getQMesh());
  bool withVelocities = false;
  bool withEigenvectors = false;
  FullBandStructure fullBandStructure =
      phononH0.populate(fullPoints, withVelocities, withEigenvectors);

  // calculate DOS via the tetrahedron method 
  auto tup1 = calcDOS(context, fullBandStructure);
  std::vector<double> energies = std::get<0>(tup1);
  std::vector<double> dos = std::get<1>(tup1);

  // save dos to an output file
  // arguments are: energies, dos, particleType, outFileName,
  //      energyUnit, dosUnit, energyConversion, chemicalPotential
  outputDOSToJSON(energies, dos, "phonon", "phonon_dos.json",
           "cm$^{-1}$", "1/(cm$^{-1}$)", ryToCmm1, 0.0);

  // Save phonon DOS to file
  if (mpi->mpiHead()) {
    std::cout << "Phonon DoS computed." << std::endl;
  }
}

/* ---------------- ElectronWannierDos --------------------*/
// Compute the Electron DOS with tetrahedron method and Wannier interpolation
void ElectronWannierDosApp::run(Context &context) {
  if (mpi->mpiHead())
    std::cout << "Starting electronic (Wannier) DoS calculation" << std::endl;

  // Read the necessary input files
  auto tup = QEParser::parseElHarmonicWannier(context);
  auto crystal = std::get<0>(tup);
  auto h0 = std::get<1>(tup);

  // first we make compute the band structure on the fine grid
  FullPoints fullPoints(crystal, context.getKMesh());
  bool withVelocities = false;
  bool withEigenvectors = false;
  FullBandStructure fullBandStructure =
      h0.populate(fullPoints, withVelocities, withEigenvectors);

  // calculate DOS via the tetrahedron method 
  auto tup1 = calcDOS(context, fullBandStructure);
  std::vector<double> energies = std::get<0>(tup1);
  std::vector<double> dos = std::get<1>(tup1);

  // Use statisticsSweep to get the chemical potential 
  // TODO do we want to do this, or would we prefer to use whatever was read in
  // by context (the value provided by QE)
  Eigen::VectorXd dummyZero(1);
  dummyZero(0) = 0.0; // set both temperature and doping to zero
  context.setTemperatures(dummyZero);
  context.setDopings(dummyZero);
  StatisticsSweep statisticsSweep(context,&fullBandStructure);
  auto stats = statisticsSweep.getCalcStatistics(0);

  // save dos to an output file
  // arguments are: energies, dos, particleType, outFileName,
  //      energyUnit, dosUnit, energyConversion, chemicalPotential
  outputDOSToJSON(energies, dos, "electron", "electron_dos.json",
           "eV", "1/eV", energyRyToEv, stats.chemicalPotential);

  // Merge different ranks and save DOS to file (mpi head only)
  if (mpi->mpiHead()) {
    std::cout << "Electronic (Wannier) DoS computed" << std::endl;
  }
}

/* ---------------- ElectronFourierDos --------------------*/
// Compute the Electron DOS with tetrahedron method and Fourier interpolation
void ElectronFourierDosApp::run(Context &context) {
  if (mpi->mpiHead())
    std::cout << "Starting electronic (Fourier) DoS calculation" << std::endl;

  // Read the necessary input files
  auto tup = QEParser::parseElHarmonicFourier(context);
  auto crystal = std::get<0>(tup);
  auto h0 = std::get<1>(tup);

  // first we make compute the band structure on the fine grid
  FullPoints fullPoints(crystal, context.getKMesh());
  bool withVelocities = false;
  bool withEigenvectors = false;
  FullBandStructure fullBandStructure =
      h0.populate(fullPoints, withVelocities, withEigenvectors);

  // calculate DOS via the tetrahedron method 
  auto tup1 = calcDOS(context, fullBandStructure); 
  std::vector<double> energies = std::get<0>(tup1); 
  std::vector<double> dos = std::get<1>(tup1);

  // Use statisticsSweep to get the chemical potential 
  // TODO do we want to do this, or would we prefer to use whatever was read in
  // by context (the value provided by QE)
  Eigen::VectorXd dummyZero(1);
  dummyZero(0) = 0.0; // set both temperature and doping to zero
  context.setTemperatures(dummyZero);
  context.setDopings(dummyZero);
  StatisticsSweep statisticsSweep(context,&fullBandStructure);
  auto stats = statisticsSweep.getCalcStatistics(0);

  // save dos to an output file
  // arguments are: energies, dos, particleType, outFileName,
  //      energyUnit, dosUnit, energyConversion, chemicalPotential
  outputDOSToJSON(energies, dos, "electron", "electron_dos.json",
           "eV", "1/eV", energyRyToEv, stats.chemicalPotential);

  // save dos to an output file
  if (mpi->mpiHead()) {
    std::cout << "Electronic (Fourier) DoS computed" << std::endl;
  }  
}

/* 
 * Generic helper function to calculate DoS via tetrahedron method 
 */ 
std::tuple<std::vector<double>, std::vector<double>> calcDOS(
        Context& context, FullBandStructure& fullBandStructure) { 

  // Form tetrahedra and fill them with eigenvalues
  TetrahedronDeltaFunction tetrahedra(fullBandStructure);

  double minEnergy = context.getDosMinEnergy();
  double maxEnergy = context.getDosMaxEnergy();
  double deltaEnergy = context.getDosDeltaEnergy();
  long numEnergies = (maxEnergy - minEnergy) / deltaEnergy + 1;

  // create instructions about how to divide up the work
  auto divs = mpi->divideWork(numEnergies);
  // the amount of values this process has to take care of
  int start = divs[0];
  int stop = divs[1];
  int workFraction = stop - start;

  std::vector<double> energies(workFraction);
  #pragma omp parallel for
  for (long i = start; i < stop; i++) {
    energies[i - start] = i * deltaEnergy + minEnergy;
  }

  // Calculate phonon density of states (DOS) [1/Ry]
  std::vector<double> dos(workFraction, 0.);  // DOS initialized to zero
  #pragma omp parallel for
  for (long i = 0; i < workFraction; i++) {
    dos[i] += tetrahedra.getDOS(energies[i]);
  }

  std::vector<double> dosTotal;
  std::vector<double> eneTotal;
  // root allocates memory to receive the gathered data
  if (mpi->mpiHead()) {
    dosTotal.resize(numEnergies);
    eneTotal.resize(numEnergies);
  }
  // all processes send their data to be gathered into dos/eneTotal
  mpi->gatherv(&energies, &eneTotal);
  mpi->gatherv(&dos, &dosTotal);
  return {eneTotal, dosTotal}; 

}

/* 
 * helper function to output dos to a json file 
 */
void outputDOSToJSON(std::vector<double> energies, std::vector<double> dos, 
                 std::string particleType, std::string outFileName, 
                 std::string energyUnit, std::string dosUnit, 
                 double energyConversion, double chemicalPotential) {

  if ( mpi->mpiHead()) {

    for (unsigned int i = 0; i < energies.size(); i++) {
      energies[i] *= energyConversion; 
      dos[i] /= energyConversion; 
    }

    // output to json 
    nlohmann::json output;
    output["energies"] = energies;
    output["dos"] = dos; 
    output["chemicalPotential"] = chemicalPotential*energyConversion;
    output["particleType"] = particleType;
    output["energyUnit"] = energyUnit;
    output["dosUnit"] = dosUnit; 
    std::ofstream o(outFileName);
    o << std::setw(3) << output << std::endl;
    o.close();
  }
}

void PhononDosApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getPhD2FileName(), "PhD2FileName");
  throwErrorIfUnset(context.getQMesh(), "qMesh");
  throwErrorIfUnset(context.getDosMinEnergy(), "dosMinEnergy");
  throwErrorIfUnset(context.getDosMaxEnergy(), "dosMaxEnergy");
  throwErrorIfUnset(context.getDosDeltaEnergy(), "dosDeltaEnergy");
  throwWarningIfUnset(context.getSumRuleD2(), "sumRuleD2");
}

void ElectronWannierDosApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");
  throwErrorIfUnset(context.getQMesh(), "kMesh");
  throwErrorIfUnset(context.getDosMinEnergy(), "dosMinEnergy");
  throwErrorIfUnset(context.getDosMaxEnergy(), "dosMaxEnergy");
  throwErrorIfUnset(context.getDosDeltaEnergy(), "dosDeltaEnergy");
  std::string crystalMsg = "crystal structure";
  throwErrorIfUnset(context.getInputAtomicPositions(), crystalMsg);
  throwErrorIfUnset(context.getInputSpeciesNames(), crystalMsg);
  throwErrorIfUnset(context.getInputAtomicSpecies(), crystalMsg);
}

void ElectronFourierDosApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");
  throwErrorIfUnset(context.getQMesh(), "kMesh");
  throwErrorIfUnset(context.getDosMinEnergy(), "dosMinEnergy");
  throwErrorIfUnset(context.getDosMaxEnergy(), "dosMaxEnergy");
  throwErrorIfUnset(context.getDosDeltaEnergy(), "dosDeltaEnergy");
  throwErrorIfUnset(context.getElectronFourierCutoff(),
                    "electronFourierCutoff");
}
