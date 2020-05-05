#include <string>
#include "dos_app.h"
#include "exceptions.h"
#include "constants.h"

#include "delta_function.h"

// Compute the DOS with the tetrahedron method
void DosApp::run(Context & context) {

	// Read the necessary input files
	auto [crystal, phononH0] = setupPhononH0(context);

	// first we make compute the band structure on the fine grid
	bool withVelocities=false;
	bool withEigenvectors = false;
	auto [fullPoints,fullBandStructure] = buildFullBandStructure(crystal,
			context.getQMesh(), phononH0, withVelocities, withEigenvectors);

	// Free some memory
	phononH0.~PhononH0();

	// Form tetrahedra and fill them with eigenvalues
	Tetrahedra tetrahedra(fullPoints, fullBandStructure);

	// Hard coded limits of energy. Later change them to user input?
	// Array of uniform frequencies to sample
	long numFreq = 401;
	std::vector<double> freq(numFreq, 0.);
	double del = 1. / ryToCmm1; //frequency spacing, THz
	for ( long i=0; i<numFreq; i++ ) {
		freq[i] = i * del;
	}

	// Calculate phonon density of states (DOS) [1/Ry]
	std::vector<double> dos(numFreq, 0.); // phonon DOS initialized to zero
	for ( long i=0; i<numFreq; i++ ) {
		dos[i] += tetrahedra.getTotalWeight(freq[i]);
	}

	// Save phonon DOS to file
	std::ofstream outfile("./phonon_dos.dat");
	outfile << "# Phonon density of states: frequency[Cmm1], Dos[a.u.]";
	for ( long i=0; i<numFreq; i++ ) {
		outfile << freq[i] * ryToCmm1 << "\t" << dos[i] << "\n";
	}
}
