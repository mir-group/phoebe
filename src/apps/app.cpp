#include <string>
#include "app.h"
#include "qe_input_parser.h"
#include "context.h"
#include "constants.h"
#include "exceptions.h"
#include "points.h"
#include "io.h"
#include "window.h"
#include "statistics.h"
#include "dos_app.h"
#include "bands_app.h"
#include "phonon_transport_app.h"
#include "utilities.h"
#include "active_bandstructure.h"
#include "bandstructure.h"

// app factory
std::unique_ptr<App> App::loadApp(std::string & choice) {
	if ( choice == "phononTransport" ) {
		return std::unique_ptr<App> (new PhononTransportApp);
	} else if ( choice == "phononDos" ) {
		return std::unique_ptr<App> (new PhononDosApp);
	} else if ( choice == "electronWannierDos" ) {
		return std::unique_ptr<App> (new ElectronWannierDosApp);
	} else if ( choice == "electronFourierDos" ) {
		return std::unique_ptr<App> (new ElectronFourierDosApp);
	} else if ( choice == "phononBands" ) {
		return std::unique_ptr<App> (new PhononBandsApp);
	} else if ( choice == "electronWannierBands" ) {
		return std::unique_ptr<App> (new ElectronWannierBandsApp);
	} else {
		return std::unique_ptr<App> (nullptr);
	}
}

void App::run(Context & context) {
	(void) context; // suppress unused variable compiler warnings
	Error e("Base class app doesn't have a run()");
}

std::tuple<Crystal, PhononH0> App::setupPhononH0(Context & context) {
	auto [crystal, phononH0] =
			parser.parsePhHarmonic(context);
	phononH0.setAcousticSumRule(context.getSumRuleD2());
	return {crystal, phononH0};
}

std::tuple<ActivePoints, ActiveBandStructure> App::restrictBandStructure(
		Context & context, FullBandStructure<FullPoints> & fullBandStructure) {

	Statistics statistics = fullBandStructure.getStatistics();

	// we create the window object
	Window window(context, statistics);

	// initialize activebandstructure
	ActiveBandStructure activeBandStructure(statistics);
	ActivePoints activePoints = activeBandStructure.buildAsPostprocessing(
			window, fullBandStructure);
	// note: activePoints should not go out of scope
	return {activePoints, activeBandStructure};
};
