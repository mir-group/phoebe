#include "common_kokkos.h"
#include "eigen.h"

void kokkosZHEEV(ComplexView3D& A, DoubleView2D& W) {
  // kokkos people didn't implement the diagonalization of matrices.
  // So, we have to do a couple of dirty tricks

#if defined(KOKKOS_ENABLE_SERIAL) || defined(KOKKOS_ENABLE_OPENMP)
  // in this case, we use Eigen to solve the diagonalization problem
  // critically, here I assume host == device!
  // hence, no need for deep_copy or Kokkos lambdas

  int M = A.extent(0); // number of matrices
  int N = A.extent(1); // matrix size is NxN

#pragma omp parallel for
  for (int i = 0; i < M; ++i) {
    // extract the block that describes the H(numW,numW) at fixed k
    ComplexView2D H = Kokkos::subview(A, i, Kokkos::ALL, Kokkos::ALL);

    // this is a pointer to the storage, N values, viewing it as std::complex
    // this hopes that kokkos::complex uses the same pattern of std::complex
    auto* storage = reinterpret_cast<std::complex<double>*>(H.data());

    // now, I feed the data to Eigen
    // BEWARE: this statement doesn't do data copy, but points directly to the
    // memory array. No out-of-bounds checks are made!
    Eigen::Map<Eigen::MatrixXcd> thisH(storage, N, N);
    // also, MatrixXcd is col-major. But it doesn't matter since H hermitian

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigenSolver(thisH);
    Eigen::VectorXd energies = eigenSolver.eigenvalues();
    Eigen::MatrixXcd eigenvectors = eigenSolver.eigenvectors();
    for (int m=0; m<N; ++m) {
      W(i, m) = energies(m);
      for (int n = 0; n < N; ++n) {
        A(i, m, n) = eigenvectors(m, n);
      }
    }
  }
#else
  Error("Kokkos@Phoebe: implement diagonalization in this architecture");
#endif
}
