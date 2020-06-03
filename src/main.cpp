#include "app.h"
#include "context.h"
#include "io.h"

int main(int argc, char** argv) {

	// here launch parallel environment

	// setup input/output

	IO io(argc, argv);
	io.welcome();

	// Read user input file

	Context context;
	context.setupFromInput(io.getInputFileName());

	// decide which app to use

	std::string appName = context.getAppName();
	std::unique_ptr<App> app = App::loadApp(appName);
	if ( app != nullptr ) {
		// launch it

		app->run(context);
	} else {
		std::cout << "No app to launch found." << std::endl;
	}

	// exiting program

	io.goodbye();

	// here close parallel environment

	return(0);
}
