#include "mpiHelper.h"
#include "gtest/gtest.h"
#include "common_kokkos.h"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  initMPI(argc,argv);
  initKokkos(argc, argv);

  // mute messages from other than root rank
  ::testing::TestEventListeners& listeners =
    ::testing::UnitTest::GetInstance()->listeners();
  if (!mpi->mpiHead()) {
    delete listeners.Release(listeners.default_result_printer());
  }

  int errCode = RUN_ALL_TESTS();

  deleteKokkos();
  mpi->finalize();

  return errCode;
}
