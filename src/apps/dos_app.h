#include <string>
#include "app.h"

class PhononDosApp : public App {
public:
	void run(Context & context);
};

class ElectronWannierDosApp : public App {
public:
	void run(Context & context);
};

class ElectronFourierDosApp : public App {
public:
	void run(Context & context);
};
