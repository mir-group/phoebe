#include "example.h"
#include "mpiHelper.h"
#include "gtest/gtest.h"
#include <Kokkos_Core.hpp>
#include "io.h"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  initMPI(argc,argv);
  Kokkos::initialize(argc, argv);

  // mute messages from other than root rank
  ::testing::TestEventListeners& listeners =
    ::testing::UnitTest::GetInstance()->listeners();
  if (!mpi->mpiHead()) {
    delete listeners.Release(listeners.default_result_printer());
  }

  int errCode = RUN_ALL_TESTS();

  Kokkos::finalize();
  mpi->finalize();

  return errCode;
}
