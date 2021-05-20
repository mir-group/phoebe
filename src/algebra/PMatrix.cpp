#ifdef MPI_AVAIL

#include "PMatrix.h"

#include "Blacs.h"
#include "constants.h"
#include "mpiHelper.h"
#include "utilities.h"

template <>
ParallelMatrix<double> ParallelMatrix<double>::prod(
    const ParallelMatrix<double>& that, const char& trans1,
    const char& trans2) {

  if(cols() != that.rows()) {
    Error("Cannot multiply matrices for which lhs.cols != rhs.rows.");
  }
  auto result = that;
  result.zeros();

  int m;
  if (trans1 == transN) {
    m = numRows_;
  } else {
    m = numCols_;
  }
  int n;
  if (trans2 == transN) {
    n = that.numCols_;
  } else {
    n = that.numRows_;
  }
  int k;
  if (trans1 == transN) {
    k = numCols_;
  } else {
    k = numRows_;
  }
  // check on k being consistent
  if (trans2 == transN) {
    assert(k == that.numRows_);
  } else {
    assert(k == that.numCols_);
  }
  double alpha = 1.;
  double beta = 0.;
  int one = 1;
  pdgemm_(&trans1, &trans2, &m, &n, &k, &alpha, mat, &one, &one, &descMat_[0],
          that.mat, &one, &one, &that.descMat_[0], &beta, result.mat, &one,
          &one, &result.descMat_[0]);
  return result;
}

template <>
ParallelMatrix<std::complex<double>> ParallelMatrix<std::complex<double>>::prod(
    const ParallelMatrix<std::complex<double>>& that, const char& trans1,
    const char& trans2) {
  ParallelMatrix<std::complex<double>> result(numRows_, numCols_,
                                              numBlocksRows_, numBlocksCols_);
  if(cols() != that.rows()) {
    Error("Cannot multiply matrices for which lhs.cols != rhs.rows.");
  }
  int m;
  if (trans1 == transN) {
    m = numRows_;
  } else {
    m = numCols_;
  }
  int n;
  if (trans2 == transN) {
    n = that.numCols_;
  } else {
    n = that.numRows_;
  }
  int k;
  if (trans1 == transN) {
    k = numCols_;
  } else {
    k = numRows_;
  }
  // check on k being consistent
  if (trans2 == transN) {
    assert(k == that.numRows_);
  } else {
    assert(k == that.numCols_);
  }
  std::complex<double> alpha = complexOne;
  std::complex<double> beta = complexZero;
  int one = 1;
  pzgemm_(&trans1, &trans2, &m, &n, &k, &alpha, mat, &one, &one, &descMat_[0],
          that.mat, &one, &one, &that.descMat_[0], &beta, result.mat, &one,
          &one, &result.descMat_[0]);
  return result;
}

template <>
std::tuple<std::vector<double>, ParallelMatrix<double>>
ParallelMatrix<double>::diagonalize() {
  if (numRows_ != numCols_) {
    Error("Cannot diagonalize non-square matrix");
  }
  if ( numBlasRows_ != numBlasCols_ ) {
    Error("Cannot diagonalize via scalapack with a non-square process grid!");
  }

  double *eigenvalues;
  allocate(eigenvalues, numRows_);

  // Make a new PMatrix to receive the output
  ParallelMatrix<double> eigenvectors(numRows_,numCols_,
                                      numBlocksRows_,numBlocksCols_);

  char jobz = 'V';  // also eigenvectors
  char uplo = 'U';  // upper triangolar
  int ia = 1;       // row index from which we diagonalize
  int ja = 1;       // row index from which we diagonalize

  int info = 0;

  // we will let pdseyv determine lwork for us. if we run it with
  // lwork = -1 and work of length 1, it will fill work with
  // an appropriate lwork number
  double *work;
  int lwork = -1;
  allocate(work, 1);
  pdsyev_(&jobz, &uplo, &numRows_, mat, &ia, &ja, &descMat_[0], eigenvalues,
          eigenvectors.mat, &ia, &ja, &descMat_[0], work, &lwork, &info);
  lwork=int(work[0]);
  delete[] work;
  allocate(work, lwork);

  // Here, we are tricking scalapack a little bit.
  // normally, one would provide the descMat for the eigenvectors matrix
  // as the 12th argument to pdsyev. However, this will result in an error,
  // because the blacsContext for eigenvectors is not the exact same reference
  // as the one for the input matrix. However, because eigenvectors is
  // the same size and block distribution as the input matrix,
  // there is no harm in using the same values for descMat.
  pdsyev_(&jobz, &uplo, &numRows_, mat, &ia, &ja, &descMat_[0], eigenvalues,
          eigenvectors.mat, &ia, &ja, &descMat_[0], work, &lwork, &info);

  if (info != 0) {
    Error("PDSYEV failed", info);
  }

  std::vector<double> eigenvalues_(numRows_);
  for (int i = 0; i < numRows_; i++) {
    eigenvalues_[i] = *(eigenvalues + i);
  }
  delete[] eigenvalues;
  delete[] work;
  // note that the scattering matrix now has different values

  return {eigenvalues_, eigenvectors};
}

template <>
std::tuple<std::vector<double>, ParallelMatrix<std::complex<double>>>
ParallelMatrix<std::complex<double>>::diagonalize() {
  if (numRows_ != numCols_) {
    Error("Can not diagonalize non-square matrix");
  }
  if ( numBlasRows_ != numBlasCols_ ) {
    Error("Cannot diagonalize via scalapack with a non-square process grid!");
  }
  double* eigenvalues = nullptr;
  eigenvalues = new double[numRows_];

  ParallelMatrix<std::complex<double>> eigenvectors(
      numRows_, numCols_, numBlocksRows_, numBlocksCols_);

  // find the value of lwork and lrwork. These are internal "scratch" arrays
  int NB = descMat_[5];
  int aZero = 0;
  int NN = std::max(std::max(numRows_, NB), 2);
  int NP0 = numroc_(&NN, &NB, &aZero, &aZero, &numBlasRows_);
  int NQ0 = numroc_(&NN, &NB, &aZero, &aZero, &numBlasCols_);
  int lwork = (NP0 + NQ0 + NB) * NB + 3 * numRows_ + numRows_ * numRows_;
  int lrwork = 2 * numRows_ + 2 * numRows_ - 2;

  // double work[lwork];
  std::complex<double>* work = nullptr;
  work = new std::complex<double>[lwork];
  std::complex<double>* rwork = nullptr;
  rwork = new std::complex<double>[lrwork];

  char jobz = 'V';  // also eigenvectors
  char uplo = 'U';  // upper triangolar
  int ia = 1;       // row index from which we diagonalize
  int ja = 1;       // row index from which we diagonalize
  int info = 0;
  pzheev_(&jobz, &uplo, &numRows_, mat, &ia, &ja, &descMat_[0], eigenvalues,
          eigenvectors.mat, &ia, &ja, &eigenvectors.descMat_[0], work, &lwork,
          rwork, &lrwork, &info);

  if (info != 0) {
    Error("PZHEEV failed", info);
  }

  std::vector<double> eigenvalues_(numRows_);
  for (int i = 0; i < numRows_; i++) {
    eigenvalues_[i] = *(eigenvalues + i);
  }
  delete[] eigenvalues;
  delete[] work;
  delete[] rwork;
  // note that the scattering matrix now has different values

  return {eigenvalues_, eigenvectors};
}

#endif  // MPI_AVAIL
