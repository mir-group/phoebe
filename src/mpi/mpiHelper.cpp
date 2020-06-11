#include "mpiHelper.h"

MPIcontroller* mpi = 0;

// A function to set up the mpi env by creating the controller object. 
void initEnv(){

        std::cout << "trying to create mpi controller" << std::endl;
        mpi = new MPIcontroller();
        std::cout << "created mpi controller" << std::endl;
        return;  
}

