#include "dos_app.h"

#include <fstream>
#include <iostream>
#include <string>

#include "constants.h"
#include "delta_function.h"
#include "electron_h0_fourier.h"
#include "exceptions.h"
#include "full_points.h"
#include "mpiHelper.h"
#include "qe_input_parser.h"
#include "utilities.h"

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

  // Form tetrahedra and fill them with eigenvalues
  TetrahedronDeltaFunction tetrahedra(fullBandStructure);

  double minEnergy = context.getDosMinEnergy();
  double maxEnergy = context.getDosMaxEnergy();
  double deltaEnergy = context.getDosDeltaEnergy();
  long numEnergies = (maxEnergy - minEnergy) / deltaEnergy + 1;

  // create instructions about how to divide up the work
  mpi->divideWork(numEnergies);
  // the amount of values this process has to take care of
  int start = mpi->workHead();
  int stop = mpi->workTail();
  int workFraction = stop - start;

  std::vector<double> energies(workFraction);
  for (long i = start; i < stop; i++) {
    energies[i - start] = i * deltaEnergy;
  }

  // Calculate phonon density of states (DOS) [1/Ry]
  std::vector<double> dos(workFraction, 0.);  // phonon DOS initialized to zero
  for (long i = 0; i < workFraction; i++) {
    dos[i] += tetrahedra.getDOS(energies[i]);
  }

  // containers for the final result, only allocated for head process
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

  // Save phonon DOS to file
  if (mpi->mpiHead()) {
    std::ofstream outfile("./phonon_dos.dat");
    outfile << "# Phonon density of states: frequency[Cmm1], Dos[1/Ry]\n";
    for (long i = 0; i < numEnergies; i++) {
      outfile << eneTotal[i] * ryToCmm1 << "\t" << dosTotal[i] << "\n";
    }
    std::cout << "Phonon DoS computed" << std::endl;
  }
}

// Compute the Electron DOS with tetrahedron method and Fourier interpolation
void ElectronWannierDosApp::run(Context &context) {
  if (mpi->mpiHead())
    std::cout << "Starting electronic DoS calculation" << std::endl;

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

  // Form tetrahedra and fill them with eigenvalues
  TetrahedronDeltaFunction tetrahedra(fullBandStructure);

  double minEnergy = context.getDosMinEnergy();
  double maxEnergy = context.getDosMaxEnergy();
  double deltaEnergy = context.getDosDeltaEnergy();
  long numEnergies = (maxEnergy - minEnergy) / deltaEnergy + 1;

  // create instructions about how to divide up the work
  mpi->divideWork(numEnergies);
  // the amount of values this process has to take care of
  int start = mpi->workHead();
  int stop = mpi->workTail();
  int workFraction = stop - start;

  std::vector<double> energies(workFraction);
// loop over the energies that this process should work on
#pragma omp parallel for
  for (int i = start; i < stop; i++) {
    energies[i - start] = i * deltaEnergy + minEnergy;
  }

  // Calculate density of states (DOS) [1/Ry], for energies on this process
  std::vector<double> dos(workFraction, 0.);  // DOS initialized to zero
#pragma omp parallel for
  for (int i = 0; i < workFraction; i++) {
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

  // Merge different ranks and save DOS to file (mpi head only)
  if (mpi->mpiHead()) {
    std::ofstream outfile("./electron_dos.dat");
    outfile << "# Electronic density of states: energy[eV], Dos[1/Ry]\n";
    for (long i = 0; i < numEnergies; i++) {
      outfile << eneTotal[i] * energyRyToEv << "\t"
              << dosTotal[i] / energyRyToEv << "\n";
    }
    std::cout << "Electronic DoS computed" << std::endl;
  }
}

// Compute the Electron DOS with tetrahedron method and Fourier interpolation
void ElectronFourierDosApp::run(Context &context) {
  if (mpi->mpiHead())
    std::cout << "Starting electronic DoS calculation" << std::endl;

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

  // Form tetrahedra and fill them with eigenvalues
  TetrahedronDeltaFunction tetrahedra(fullBandStructure);

  double minEnergy = context.getDosMinEnergy();
  double maxEnergy = context.getDosMaxEnergy();
  double deltaEnergy = context.getDosDeltaEnergy();
  long numEnergies = (maxEnergy - minEnergy) / deltaEnergy + 1;

  // create instructions about how to divide up the work
  mpi->divideWork(numEnergies);
  // the amount of values this process has to take care of
  int start = mpi->workHead();
  int stop = mpi->workTail();
  int workFraction = stop - start;

  std::vector<double> energies(workFraction);
  for (long i = start; i < stop; i++) {
    energies[i - start] = i * deltaEnergy + minEnergy;
  }

  // Calculate phonon density of states (DOS) [1/Ry]
  std::vector<double> dos(workFraction, 0.);  // DOS initialized to zero
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

  // save dos to an output file
  if (mpi->mpiHead()) {
    std::ofstream outfile("./electron_dos.dat");
    outfile << "# Electronic density of states: energy[eV], Dos[1/Ry]\n";
    for (long i = 0; i < numEnergies; i++) {
      outfile << eneTotal[i] * energyRyToEv << "\t"
              << dosTotal[i] / energyRyToEv << "\n";
    }
    std::cout << "Electronic DoS computed" << std::endl;
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
