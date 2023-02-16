#ifndef EL_EL_PLOT_APP_H
#define EL_El_PLOT_APP_H

#include <string>
#include "app.h"

/** Main driver for the calculation
 */
class ElElCouplingPlotApp: public App {
public:
  void run(Context &context) override;
  void checkRequirements(Context &context) override;
};

#endif
