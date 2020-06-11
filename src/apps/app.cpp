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
#include "polarization_app.h"
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
	} else if ( choice == "electronFourierBands" ) {
		return std::unique_ptr<App> (new ElectronFourierBandsApp);
	} else if ( choice == "electronPolarization" ) {
		return std::unique_ptr<App> (new ElectronPolarizationApp);
	} else {
		return std::unique_ptr<App> (nullptr);
	}
}

void App::run(Context & context) {
	(void) context; // suppress unused variable compiler warnings
	Error e("Base class app doesn't have a run()");
}
