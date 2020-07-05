#include "example.h"
#include "gtest/gtest.h"
#include "mpiHelper.h"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  initMPI();

  int errCode = RUN_ALL_TESTS();

  mpi->finalize();
  
  return errCode;
}
