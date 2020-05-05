#include "app.h"
#include "context.h"
#include "io.h"

int main(int argc, char** argv) {

	// here launch parallel environment

	IO io(argc, argv);

	// Read user input file

	Context context;
	context.setupFromInput(io.getInputFileName());

	// decide which app to use

	std::string appName = context.getAppName();
	App * app = App::loadApp(appName);

	// then launch it

	app->run(context);

	// here we should close the parallel environment

	return(0);
}
