#include "mpiHelper.h"

MPIcontroller* mpi = 0;

// A function to set up the mpi env by creating the controller object. 
void initMPI(){

        mpi = new MPIcontroller();
        return;  
}

