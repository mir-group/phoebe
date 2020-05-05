#ifndef APP_H
#define APP_H

#include <string>
#include "context.h"
#include "phonon_h0.h"
#include "electron_h0_fourier.h"
#include "electron_h0_wannier.h"
#include "qe_input_parser.h"
#include "statistics.h"

/** Base class for apps
 *
 */
class App {
public:
	// technically, this is called a factory method
	static App * loadApp(std::string choice);
	virtual void run(Context & context);
protected:
	std::tuple<Crystal, PhononH0> setupPhononH0(Context & context);
	std::tuple<Crystal, ElectronH0Fourier> setupElectronH0Fourier(Context & context);
	ElectronH0Wannier setupElectronH0Wannier(Context & context);
	std::tuple<ActivePoints, ActiveBandStructure> restrictBandStructure(
			Context & context, FullBandStructure & fullBandStructure);
	std::tuple<FullPoints, FullBandStructure> buildFullBandStructure(
			Crystal & crystal, Eigen::Vector3i mesh, HarmonicHamiltonian & h0,
			bool & withVelocities, bool & withEigenvectors);

	QEParser qeParser;
};

#endif
