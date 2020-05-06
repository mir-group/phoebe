#include <string>
#include "dos_app.h"
#include "exceptions.h"
#include "constants.h"

#include "delta_function.h"

// Compute the DOS with the tetrahedron method
void DosApp::run(Context & context) {
	std::cout << "Starting DoS calculation" << std::endl;

	// Read the necessary input files
	auto [crystal, phononH0] = setupPhononH0(context);
//	std::tuple<Crystal, PhononH0> t = setupPhononH0(context);
//	Crystal crystal = std::get<0>(t);
//	PhononH0 phononH0 = std::get<1>(t);

	// first we make compute the band structure on the fine grid
	FullPoints fullPoints(crystal, context.getQMesh());
	bool withVelocities = false;
	bool withEigenvectors = false;
	auto fullBandStructure = buildFullBandStructure(fullPoints, phononH0,
			withVelocities, withEigenvectors);

	// Form tetrahedra and fill them with eigenvalues
	Tetrahedra tetrahedra(fullPoints, fullBandStructure);
	// Hard coded limits of energy. Later change them to user input?
	// Array of uniform frequencies to sample

	double minEnergy = context.getDosMinEnergy();
	double maxEnergy = context.getDosMaxEnergy();
	double deltaEnergy = context.getDosDeltaEnergy();
	long numEnergies = (maxEnergy - minEnergy) / deltaEnergy + 1;
	std::vector<double> energies(numEnergies);
	for ( long i=0; i<numEnergies; i++ ) {
		energies[i] = i * deltaEnergy;
	}

	// Calculate phonon density of states (DOS) [1/Ry]
	std::vector<double> dos(numEnergies, 0.); // phonon DOS initialized to zero
	for ( long i=0; i<numEnergies; i++ ) {
		dos[i] += tetrahedra.getDOS(energies[i]);
	}

	// Save phonon DOS to file
	std::ofstream outfile("./phonon_dos.dat");
	outfile << "# Phonon density of states: frequency[Cmm1], Dos[1/Ry]\n";
	for ( long i=0; i<numEnergies; i++ ) {
		outfile << energies[i] * ryToCmm1 << "\t" << dos[i] << "\n";
	}
	std::cout << "DoS computed" << std::endl;
}
