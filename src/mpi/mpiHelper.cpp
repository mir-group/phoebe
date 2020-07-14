#include "mpiHelper.h"
#include <stdio.h>
#include <iostream>

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
  std::cout << "Initialized with:\n";
#ifdef MPI_AVAIL
  std::cout << "MPI Processes " << mpi->getSize() << std::endl;
#endif
#ifdef OMP_AVAIL 
  std::cout << "OMP Threads " << omp_get_max_threads() << std::endl;
#endif
}
