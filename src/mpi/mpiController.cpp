#include <vector>
#include <complex>
#include <chrono>
#include <iostream>
#include "mpiController.h"

#ifdef MPI_AVAIL 
#include <mpi.h>
#endif

// default constructor
MPIcontroller::MPIcontroller(){

        fprintf(stdout,"inside the controller");	
	#ifdef MPI_AVAIL
	// start the MPI environment
	MPI_Init(NULL, NULL);
        fprintf(stdout,"initialized the env");  
	// get the number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	// get rank of current process
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	// start a timer
	startTime = MPI_Wtime();
	fprintf(stdout,"started the timer"); 
	// set this so that MPI returns errors and lets us handle them, rather
	// than using the default, MPI_ERRORS_ARE_FATAL
	MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
	#else
 
	startTime = std::chrono::steady_clock::now();
	#endif
        fprintf(stdout, "started the timer and finished");
}

// default destructor
MPIcontroller::~MPIcontroller(){
	#ifdef MPI_AVAIL
	if (!MPI::Is_finalized()) finalize();
	#endif
}

// TODO: any other stats would like to output here? 
void MPIcontroller::finalize() const {
	#ifdef MPI_AVAIL
	fprintf(stdout, "Final time for rank %3d: %3f\n ", rank, MPI_Wtime() - startTime );
	MPI_Finalize();
	#else
        std::cout << "Final time for rank 0 : " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count() << " mu.secs" << std::endl;
	#endif
}

// Utility functions  -----------------------------------

// get the error string and print it to stderr before returning
void MPIcontroller::errorReport(int errCode) const{
	#ifdef MPI_AVAIL
	char errString[BUFSIZ];
	int lengthOfString;
	
	MPI_Error_string(errCode, errString, &lengthOfString);
	fprintf(stderr, "Error from rank %3d: %s\n", rank, errString);
	MPI_Abort(MPI_COMM_WORLD, errCode);
	#else 
	// TODO: how are we throwing non-mpi errors? 
	#endif
}

void MPIcontroller::time() const{
	#ifdef MPI_AVAIL
	fprintf(stdout, "Time for rank %3d : %3f\n", rank, MPI_Wtime() - startTime );
	#else
	std::cout << "Time for rank 0 :" << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - startTime).count() << " secs" << std::endl;
	#endif
}

// Labor division functions -----------------------------------------i
void MPIcontroller::divideWork(size_t numTasks) {
        // each should be nranks long
        workDivisionHeads.resize(size);  
        workDivisionTails.resize(size); 
        //size_t numDivisons = numTasks/size + (numTasks % size != 0); // ceiling of work/ranks
        for(int r = 0; r<size; r++) {
                workDivisionHeads[r] = (numTasks * r)/size;
                workDivisionTails[r] = (numTasks * (r+1))/size;
        }
}

size_t MPIcontroller::workHead() { return workDivisionHeads[rank]; }
size_t MPIcontroller::workTail() { return workDivisionTails[rank]; } 

