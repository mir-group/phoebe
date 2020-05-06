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
#include "phonon_transport_app.h"

// app factory
std::unique_ptr<App> App::loadApp(std::string choice) {
	if ( choice == "phononTransport" ) {
		return std::unique_ptr<App> (new PhononTransportApp);
	} else if ( choice == "phononDos" ) {
		return std::unique_ptr<App> (new PhononDosApp);
	} else if ( choice == "electronDos" ) {
		return std::unique_ptr<App> (new ElectronDosApp);
	}
}

void App::run(Context & context) {
}

std::tuple<Crystal, PhononH0> App::setupPhononH0(Context & context) {
	auto [crystal, phononH0] =
			qeParser.parsePhHarmonic(context);
	phononH0.setAcousticSumRule(context.getSumRuleD2());
	return {crystal, phononH0};
}

std::tuple<Crystal, ElectronH0Fourier> App::setupElectronH0Fourier(Context & context) {
	return qeParser.parseElHarmonicFourier(context);
}

ElectronH0Wannier App::setupElectronH0Wannier(Context & context) {
	return qeParser.parseElHarmonicWannier(context);
}

std::tuple<ActivePoints, ActiveBandStructure> App::restrictBandStructure(
		Context & context, FullBandStructure & fullBandStructure) {

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
