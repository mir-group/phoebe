#ifndef APP_H
#define APP_H

#include "active_bandstructure.h"
#include "bandstructure.h"
#include "context.h"
#include "electron_h0_fourier.h"
#include "electron_h0_wannier.h"
#include "phonon_h0.h"
#include "qe_input_parser.h"
#include "statistics.h"
#include <memory>
#include <string>

/** Base class for apps
 *
 */
class App {
public:
  // technically, this is called a factory method
  static std::unique_ptr<App> loadApp(std::string &choice);
  // note: auto_ptr transfers the ownership to the returned pointer itself
  // in this way, the returned app is correctly destructed when out of scope
  virtual void run(Context &context);

protected:
  std::tuple<Crystal, PhononH0> setupPhononH0(Context &context);
  //	std::tuple<ActivePoints, ActiveBandStructure> restrictBandStructure(
  //			Context & context,
  //			FullBandStructure<FullPoints> & fullBandStructure);
  QEParser parser;
};

#endif
