#ifndef MPIASSIST_H
#define MPIASSIST_H

#include "mpiController.h"

// Currently, just using this file to make the mpiController globally available. 
// Other contents may be added in the future. 
extern MPIcontroller* mpi; 

// Function to initialize MPI environment
void initMPI(); 

// Prints info about omp num threads and mpi processes
void parallelInfo(); 

#endif
