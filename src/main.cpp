#include "app.h"
#include "context.h"
#include "io.h"
#include "mpi/mpiHelper.h"

int main(int argc, char** argv) {

	// here launch parallel environment
	// Call proxy function from MPI Helper, which makes mpi object 
	// globally available. 
        initMPI(); 

	// setup input/output
	IO io(argc, argv);
	if(mpi->mpiHead()) io.welcome();

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
	if(mpi->mpiHead()) io.goodbye();

	// here close parallel environment
        mpi->barrier(); 
        mpi->finalize(); 

	return(0);
}
