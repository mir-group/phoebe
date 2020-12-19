#ifndef DOS_APP_H
#define DOS_APP_H

#include "app.h"
#include <string>

class PhononDosApp : public App {
public:
  void run(Context &context) override;
  void checkRequirements(Context &context) override;
};

class ElectronWannierDosApp : public App {
public:
  void run(Context &context) override;
  void checkRequirements(Context &context) override;
};

class ElectronFourierDosApp : public App {
public:
  void run(Context &context) override;
  void checkRequirements(Context &context) override;
};

#endif
