#ifndef EWTRANSPORTAPP_H
#define EWTRANSPORTAPP_H

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
  void runVariationalMethod(Context &context,
                            Crystal &crystal,
                            StatisticsSweep &statisticsSweep,
                            ActiveBandStructure &bandStructure,
                            ElScatteringMatrix &scatteringMatrix);
};

#endif
