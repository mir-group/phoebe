#ifndef DOSAPP_H
#define DOSAPP_H

#include <string>
#include "app.h"

class PhononDosApp : public App {
public:
	void run(Context & context);
	void checkRequirements(Context & context);
};

class ElectronWannierDosApp : public App {
public:
	void run(Context & context);
	void checkRequirements(Context & context);
};

class ElectronFourierDosApp : public App {
public:
	void run(Context & context);
	void checkRequirements(Context & context);
};

#endif
