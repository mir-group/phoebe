#ifndef APP_H
#define APP_H

#include <string>
#include <memory>
#include "context.h"
#include "phonon_h0.h"
#include "electron_h0_fourier.h"
#include "electron_h0_wannier.h"
#include "qe_input_parser.h"
#include "statistics.h"
#include "bandstructure.h"
#include "active_bandstructure.h"

/** Base class for apps
 *
 */
class App {
public:
	// technically, this is called a factory method
	static std::unique_ptr<App> loadApp(std::string & choice);
	// note: auto_ptr transfers the ownership to the returned pointer itself
	// in this way, the returned app is correctly destructed when out of scope
	virtual void run(Context & context);
};

#endif
