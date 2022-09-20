#include "mpiHelper.h"
#include <iostream>

#ifdef OMP_AVAIL
#include "omp.h"
#endif

MPIcontroller *mpi = nullptr;

// A function to set up the mpi env by creating the controller object.
void initMPI(int argc, char *argv[]) {
  mpi = new MPIcontroller(argc, argv);
}

void deleteMPI() {
  mpi->finalize();
  delete mpi;
}

void parallelInfo() {
  if (mpi->mpiHead()) {
    std::cout << "Initialized with:\n";
#ifdef MPI_AVAIL
    if (mpi->hasPools()) {
      std::cout << "MPI Processes " << mpi->getSize() << " (" << mpi->getSize(mpi->intraPoolComm) << " processes per pool)" << std::endl;
    } else {
      std::cout << "MPI Processes " << mpi->getSize() << std::endl;
    }
#endif
#ifdef OMP_AVAIL
    std::cout << "OMP Threads " << omp_get_max_threads() << std::endl;
#endif
  }
}
