#include "constants.h"
#include "eigen.h"
#include "mpiHelper.h"
#include "gtest/gtest.h"

TEST(MPITest, AllReduceSum) {

  int size = mpi->getSize();
  int rank = mpi->getRank();

  // test allReduceSum on a std::vector -------------------------
  std::vector<double> x(size, 0.);
  x[rank] = 1.;
  // the vector has only one element different from zero
  mpi->allReduceSum(&x);
  // now the vector should contain 1 everywhere
  for (double y : x) {
    EXPECT_EQ(y, 1.);
  }

  // test allReduceSum on a Eigen::MatrixXd ---------------------
  Eigen::MatrixXd x2(size, size);
  x2.setZero();
  for (long i = 0; i < size; i++) {
    x2(rank, i) = 1.;
  }
  // the vector has only one element different from zero
  mpi->allReduceSum(&x2);
  // now the vector should contain 1 everywhere
  for (long i = 0; i < size; i++) {
    for (long j = 0; j < size; j++) {
      EXPECT_EQ(x2(i, j), 1.);
    }
  }

  // test allReduceSum on an Eigen::VectorXi ---------------------
  Eigen::VectorXi x3(size);
  x3.setZero();
  for (long i = 0; i < size; i++) {
    x3(rank, i) = 1;
  }
  // the vector has only one element different from zero
  mpi->allReduceSum(&x3);
  // now the vector should contain 1 everywhere
  for (long i = 0; i < size; i++) {
    EXPECT_EQ(x3(i), 1);
  }

  // std::vector test -------------------------------------------
  std::vector<std::complex<double>> x4(size, complexZero);
  x4[rank] = complexI;
  // the vector has only one element different from zero
  mpi->allReduceSum(&x4);
  // now the vector should contain 1 everywhere
  for (std::complex<double> y : x4) {
    EXPECT_EQ(y, complexI);
  }

  // Test allReduce Sum on an Eigen::Tensor ------------------------
  Eigen::Tensor<double, 5> x5(size, size, size, size, size);
  x5.setZero();

  // fill every element with one
  for (long i = 0; i < size; i++) {
    for (long j = 0; j < size; j++) {
      for (long k = 0; k < size; k++) {
        for (long m = 0; m < size; m++) {
          for (long n = 0; n < size; n++) {
            x5(i, j, k, m, n) = i * j * k * m * n;
          }
        }
      }
    }
  }
  // the vector has only one element different from zero
  mpi->allReduceSum(&x5);
  // now the vector should contain nRanks everywhere
  for (long i = 0; i < size; i++) {
    for (long j = 0; j < size; j++) {
      for (long k = 0; k < size; k++) {
        for (long m = 0; m < size; m++) {
          for (long n = 0; n < size; n++) {
            EXPECT_EQ(x5(i, j, k, m, n), i * k * j * m * n * mpi->getSize());
          }
        }
      }
    }
  }
}
