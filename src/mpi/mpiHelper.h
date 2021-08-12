#ifndef MPI_HELPER_H
#define MPI_HELPER_H

#include "mpiController.h"

// Currently, just using this file to make the mpiController globally available. 
// Other contents may be added in the future. 
extern MPIcontroller* mpi;

// Function to initialize MPI environment
void initMPI(int argc, char *argv[]);

// Function to delete MPI environment
void deleteMPI();

// Prints info about omp num threads and mpi processes
void parallelInfo(); 

#endif
