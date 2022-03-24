#ifndef EW_TRANSPORT_APP_H
#define EW_TRANSPORT_APP_H

#include <string>
#include "app.h"
#include "el_scattering.h"

/** Main driver for the transport calculation
 */
class ElectronWannierTransportApp: public App {
public:
  void run(Context &context) override;
  void checkRequirements(Context &context) override;

private:
  /** Method for running the variational solver of the electron BTE
   */
  static void runVariationalMethod(Context &context,
                            Crystal &crystal,
                            StatisticsSweep &statisticsSweep,
                            ActiveBandStructure &bandStructure,
                            ElScatteringMatrix &scatteringMatrix);
};

// TODO how can i have the test access this wihtout making it public
void unfoldLinewidths(Context& context, ElScatteringMatrix& oldMatrix,
                         ActiveBandStructure& bandStructure,
                         StatisticsSweep& statisticsSweep, HarmonicHamiltonian& electronH0);


#endif
