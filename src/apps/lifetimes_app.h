#ifndef LIFETIMESAPP_H
#define LIFETIMESAPP_H

#include <string>
#include "app.h"

/** Driver for computing electronic lifetimes on a path
 */
class ElLifetimesApp: public App {
public:
  void run(Context &context);
  void checkRequirements(Context &context);
};

/** Driver for computing phonon lifetimes on a path
 */
class PhLifetimesApp: public App {
public:
  void run(Context &context);
  void checkRequirements(Context &context);
};

#endif
