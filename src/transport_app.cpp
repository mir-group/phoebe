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

	QEParser qeParser;

	if ( context.getCalculation() == "electron-phonon" ) {
		auto [crystalEl, electronH0Fourier] =
				qeParser.parseElHarmonicFourier(context);
	}

	auto [crystalPh, phononH0] =
			qeParser.parsePhHarmonic(context.getPhD2FileName());
	phononH0.setAcousticSumRule(context.getSumRuleD2());

	// Now, we build the harmonic phonon properties
	Statistics phStatistics(Statistics::bose);

	// first we make compute the band structure on the fine grid
	FullPoints fullQPoints(crystalPh, context.getQMesh());
	bool withVelocities=true, withEigenvectors=true;
	FullBandStructure fullPhBandStructure(phononH0.getNumBands(), phStatistics,
			withVelocities, withEigenvectors, &fullQPoints);
	fullPhBandStructure.populate(phononH0);

	// then we apply a filter to retain only useful energies
	Window phononWindow(context, phStatistics);
	ActiveBandStructure phBandStructure(phStatistics);
	ActivePoints phQPoints = phBandStructure.buildAsPostprocessing(phononWindow,
			fullPhBandStructure);
	// note: phQPoints should not go out of scope
	// instead window, phStatistics, FullPhBandStructure and FullPoints
	// can all be left go out of scope

	/** out of this TransportApp(), we have to store somewhere:
	 * phQPoints, phBandStructure
	 * CrystalEl, CrystalPh: need both to be referenced
	 * context
	 * phononH0, ElectronH0 are not strictly needed anymore
	 * (note: we need phononH0 and ElectronH0,
	 * for example, for the U matrices or Z charges in the polar correction)
	 *
	 */


};
