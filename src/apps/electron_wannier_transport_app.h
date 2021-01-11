#ifndef EWTRANSPORTAPP_H
#define EWTRANSPORTAPP_H

#include <string>
#include "app.h"

/** Main driver for the transport calculation
 */
class ElectronWannierTransportApp: public App {
public:
  void run(Context &context) override;
  void checkRequirements(Context &context) override;
};

#endif
