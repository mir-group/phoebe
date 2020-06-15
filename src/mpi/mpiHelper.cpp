#include "mpiHelper.h"
#include <stdio.h>

#ifdef OMP_AVAIL
#include "omp.h"
#endif

MPIcontroller* mpi = 0;

// A function to set up the mpi env by creating the controller object. 
void initMPI(){

        mpi = new MPIcontroller();
        return;  
}

void parallelInfo() { 
        fprintf(stdout,"\tIntialized with: \n");
        #ifdef MPI_AVAIL
        fprintf(stdout,"\tMPI Processes\t %d \n", mpi->getSize());
        #endif
        #ifdef OMP_AVAIL 
        fprintf(stdout,"\tOMP Threads  \t %d \n", omp_get_max_threads());
        #endif
}
