#include "exceptions.h"
#include "mpiHelper.h"
#include <iostream>
#include <string>

Error::Error(const std::string &errMessage, const int &errCode) {
  if (errCode != 0) {
    if (mpi->mpiHead()) {
      std::cout << "Error!" << std::endl;
      std::cout << errMessage << std::endl;
    }
    mpi->barrier();
    mpi->finalize();
    exit(errCode);
  }
}

Warning::Warning(const std::string &errMessage) {
  if (!mpi->mpiHead())
    return;
  std::cout << "WARNING: " << errMessage << std::endl;
}
