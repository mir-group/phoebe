#ifndef BANDSAPP_H
#define BANDSAPP_H

#include "app.h"
#include "path_points.h"

class PhononBandsApp: public App {
public:
    void run(Context &context);
    void checkRequirements(Context &context);
};

class ElectronWannierBandsApp: public App {
public:
    void run(Context &context);
    void checkRequirements(Context &context);
};

class ElectronFourierBandsApp: public App {
public:
    void run(Context &context);
    void checkRequirements(Context &context);
};

void outputBandsToJSON(FullBandStructure& fullBandStructure,
                       Context& context, PathPoints& pathPoints,
                       std::string outFileName);

#endif
