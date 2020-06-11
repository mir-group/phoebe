#ifndef BANDSAPP_H
#define BANDSAPP_H

#include "app.h"

class PhononBandsApp : public App {
public:
	void run(Context & context);
	void checkRequirements(Context & context);
};

class ElectronWannierBandsApp : public App {
public:
	void run(Context & context);
	void checkRequirements(Context & context);
};

class ElectronFourierBandsApp : public App {
public:
	void run(Context & context);
	void checkRequirements(Context & context);
};

#endif
