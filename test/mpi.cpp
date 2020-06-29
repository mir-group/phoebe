#include "gtest/gtest.h"
#include "mpiHelper.h"
#include "eigen.h"

TEST(MPITest, AllReduceSum) {

  int size = mpi->getSize();
  int rank = mpi->getRank();

  // test allReduceSum on a std::vector

  std::vector<double> x(size, 0.);
  x[rank] = 1.;
  // the vector has only one element different from zero
  mpi->allReduceSum(&x);
  // now the vector should contain 1 everywhere
  for ( double y : x ) {
    EXPECT_EQ(y, 1.);
  }

  // test allReduceSum on a Eigen::MatrixXd

  Eigen::MatrixXd x2(size,size);
  x2.setZero();
  for ( long i=0; i<size; i++ ) {
      x2(rank,i) = 1.;
  }
  // the vector has only one element different from zero
  mpi->allReduceSum(&x2);
  // now the vector should contain 1 everywhere
  for ( long i=0; i<size; i++ ) {
      for ( long j=0; j<size; j++ ) {
          EXPECT_EQ(x2(i,j), 1.);
      }
  }


}
