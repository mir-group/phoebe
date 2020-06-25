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

	#ifdef MPI_AVAIL
	// start the MPI environment
        MPI_Init(NULL, NULL);

        // set this so that MPI returns errors and lets us handle them, rather
        // than using the default, MPI_ERRORS_ARE_FATAL
        MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

	// get the number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	// get rank of current process
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	// start a timer
	startTime = MPI_Wtime();

	#else
        // To maintain consistency when running in serial
        size = 1;
        rank = 0;  
	startTime = std::chrono::steady_clock::now();
	#endif
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
	if(rank==0) fprintf(stdout, "Final time for MPI head: %3f\n ", MPI_Wtime() - startTime );
	MPI_Finalize();
	#else
        std::cout << "Total runtime: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count()*1e-6 << " secs" << std::endl;
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
// Asynchronous support functions -----------------------------------------
void MPIcontroller::barrier() const{ 
        #ifdef MPI_AVAIL
        int errCode; 
        errCode = MPI_Barrier(MPI_COMM_WORLD); 
        if(errCode != MPI_SUCCESS) {  errorReport(errCode); }
        #endif
}

// Labor division functions -----------------------------------------
std::vector<int> MPIcontroller::divideWork(size_t numTasks) {
        // return a vector of the start and stop points for task division
        std::vector<int> divs(2); 
        divs[0] = (numTasks * rank)/size;
        divs[1] = (numTasks * (rank+1))/size;
        return divs; 
}

std::vector<int> MPIcontroller::divideWorkIter(size_t numTasks) {
        // return a vector of the start and stop points for task division
        std::vector<int> divs;
        int start = (numTasks * rank)/size;
        int stop = (numTasks * (rank+1))/size;
        for( int i = start; i <= stop; i++ ) divs.push_back(i);
        return divs;
}
