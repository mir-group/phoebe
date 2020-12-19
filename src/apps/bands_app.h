#ifndef BANDS_APP_H
#define BANDS_APP_H

#include "app.h"
#include "path_points.h"

class PhononBandsApp : public App {
public:
  void run(Context &context) override;
  void checkRequirements(Context &context) override;
};

class ElectronWannierBandsApp : public App {
public:
  void run(Context &context) override;
  void checkRequirements(Context &context) override;
};

class ElectronFourierBandsApp : public App {
public:
  void run(Context &context) override;
  void checkRequirements(Context &context) override;
};

void outputBandsToJSON(FullBandStructure &fullBandStructure, Context &context,
                       PathPoints &pathPoints, std::string outFileName);

#endif
