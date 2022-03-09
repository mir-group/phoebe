#include "common_kokkos.h"
#include "gtest/gtest.h"

/** Here I estimate the mass at the top of the valence band of silicon
 */
TEST(Kokkos, KokkosDeviceManager) {
  // check that the default device memory is 16Gb
  double x = kokkosDeviceMemory->getTotalMemory();
  ASSERT_EQ(x, 16.0e9);

  // at the start, available memory is approx all the total memory
  double y = kokkosDeviceMemory->getAvailableMemory();
  ASSERT_NEAR(x, y, 1.0e-2); // at the start, they should be about the same

  // check we remember we added memory
  double z = 2.0e9;
  kokkosDeviceMemory->addDeviceMemoryUsage(z);
  ASSERT_NEAR(kokkosDeviceMemory->getAvailableMemory(), y-z, 1.0e-2);

  // check that if we remove memory, we go back to where we started
  kokkosDeviceMemory->removeDeviceMemoryUsage(z);
  ASSERT_NEAR(kokkosDeviceMemory->getAvailableMemory(), y, 1.0e-2);

  std::vector<int> test = {0,1,2,3,4,5,6,7,8,9,10};
  int batchSize = 3;
  auto batches = kokkosDeviceMemory->splitToBatches(test, batchSize);

  // test: all indices inside batches should give test back
  std::vector<int> verify;
  for (auto batch : batches) {
    for (int b : batch) {
      verify.push_back(b);
    }
  }
  for (unsigned int i=0; i<test.size(); ++i) {
    ASSERT_EQ(verify[i], test[i]);
  }
}
