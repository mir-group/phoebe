#ifndef MPIHELPER_H
#define MPIHELPER_H

#include "mpiController.h"

// Currently, just using this file to make the mpiController globally available. 
// Other contents may be added in the future. 
extern MPIcontroller* mpi; 

// Function to initialize MPI environment
void initEnv(); 


#endif
