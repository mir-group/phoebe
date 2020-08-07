#include "context.h"
#include "example.h"
#include "mpiHelper.h"
#include "gtest/gtest.h"
#include <Kokkos_Core.hpp>

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  initMPI();
  Kokkos::initialize(argc, argv);

  int errCode = RUN_ALL_TESTS();

  Kokkos::finalize();
  mpi->finalize();

  return errCode;
}
