#include <string>
#include "transport_app.h"
#include "qe_input_parser.h"
#include "context.h"
#include "constants.h"
#include "exceptions.h"
#include "points.h"
#include "io.h"
#include "window.h"
#include "statistics.h"

TransportApp::TransportApp(int argc, char** argv) {
	IO io(argc, argv);

	std::string inputFileName = io.getInputFileName();

	std::cout << "Reading from input file: " << inputFileName << endl;

	// Read user input file

	Context context;
	context.setupFromInput(inputFileName);

	// Read the necessary input files

//	Eigen::VectorXd q(3);
//	q << 0.,0.,0.1;
//	auto [energies, eigenvectors] = phononH0.diagonalize(q);
//	std::cout << energies.transpose() * ryToCmm1 << std::endl;

//	Eigen::Vector3i mesh;
//	mesh << 4, 4, 4;
//	Points ps(crystal, mesh);

//
//	//	TODO: 'elph' shoudn't be a string, we should use a dictionary
//	//	and store which are the allowed values of calculations.
//

	QEParser qeParser;

	if ( context.getCalculation() == "electron-phonon" ) {
		auto [crystalEl, electronH0Fourier] =
				qeParser.parseElHarmonicFourier(context);
	}

	auto [crystalPh, phononH0] =
			qeParser.parsePhHarmonic(context.getPhD2FileName());
	phononH0.setAcousticSumRule(context.getSumRuleD2());
	phononH0.setAcousticSumRule(context.getSumRuleD2());

	// Now, we build the harmonic phonon properties
	Statistics phStatistics(Statistics::bose);

	// first we make compute the band structure on the fine grid
	FullPoints fullQPoints(crystalPh, context.getKMesh());
	bool withVelocities=true, withEigenvectors=true;
	FullBandStructure fullPhBandStructure(phononH0.getNumBands(), phStatistics,
			withVelocities, withEigenvectors, &fullQPoints);
	fullPhBandStructure.populate(phononH0);

	// then we apply a filter to retain only useful energies
	Window phononWindow(context, phStatistics);
	ActiveBandStructure phBandStructure(phStatistics);
	phBandStructure.buildAsPostprocessing(phononWindow, fullPhBandStructure);
};
