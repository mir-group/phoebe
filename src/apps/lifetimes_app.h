#ifndef LIFETIMES_APP_H
#define LIFETIMES_APP_H

#include <string>
#include "app.h"

/** Driver for computing electronic lifetimes on a path
 */
class ElectronLifetimesApp: public App {
public:
  void run(Context &context) override;
  void checkRequirements(Context &context) override;
};

/** Driver for computing phonon lifetimes on a path
 */
class PhononLifetimesApp: public App {
public:
  void run(Context &context) override;
  void checkRequirements(Context &context) override;
};

#endif
