#include "Matrix.h"

// Explict specialization of BLAS matrix-matrix mult for Matrix<complex<double>>
template <>
Matrix<std::complex<double>> Matrix<std::complex<double>>::prod(
    const Matrix<std::complex<double>>& that, const char& trans1,
    const char& trans2) {
  assert(cols() == that.rows());
  Matrix<std::complex<double>> ret(rows(), that.cols());  // newly sized matrix
  // throw away variables
  std::complex<double> alpha(1.0, 0.0);
  std::complex<double> beta(0.0, 0.0);
  zgemm_(trans1, trans2, rows(), that.cols(), cols(), alpha, mat,
         rows(), that.mat, that.rows(), beta, ret.mat, rows());
  return ret;
}

// Explicit specializiation of BLAS matrix-matrix mult for Matrix<double>
template <>
Matrix<double> Matrix<double>::prod(const Matrix<double>& that,
                                    const char& trans1, const char& trans2) {
  assert(cols() == that.rows());
  Matrix<double> ret(rows(), that.cols());  // newly sized matrix
  // throw away variables
  double alpha = 1.0;
  double beta = 0.0;
  dgemm_(trans1, trans2, rows(), that.cols(), cols(), alpha, mat,
         rows(), that.mat, that.rows(), beta, ret.mat, rows());
  return ret;
}

// Diagonalize a complex double hermitian matrix
template <>
std::tuple<std::vector<double>, Matrix<std::complex<double>>>
Matrix<std::complex<double>>::diagonalize() {
  assert(nRows == nCols);  // needs to be square

  std::vector<double> eigvals(nRows);
  Matrix<std::complex<double>> eigvecs(nRows, nCols);

  // throw away variables
  Matrix<std::complex<double>> temp;
  temp = *this;  // copy

  char jobz = 'V';
  char uplo = 'U';
  int lwork = std::max(1, 2 * nRows - 1) * 2;
  std::vector<std::complex<double>> work(lwork);
  int rworkSize = std::max(1, 3 * nRows - 2);
  std::vector<double> rwork(rworkSize);
  int info;

  zheev_(&jobz, &uplo, &nRows, temp.mat, &nRows, &eigvals[0], &work[0], &lwork,
         &rwork[0], &info);

  assert(info == 0);  // if it doesn't =0, there was an error. Different errors
                      // for info< or > 0.

  return {eigvals, eigvecs};
}

// Diagonalize for real double symmetric matrix
template <>
std::tuple<std::vector<double>, Matrix<double>> Matrix<double>::diagonalize() {
  assert(nRows == nCols);  // needs to be square
  std::vector<double> eigvals(nRows);
  Matrix<double> eigvecs = *this;
  // throw away variables

  char jobz = 'V';
  char uplo = 'U';
  int lwork = std::max(1, 3 * nRows - 1) * 2;
  std::vector<double> work(lwork);
  int info;

  dsyev_(&jobz, &uplo, &nRows, eigvecs.mat, &nRows, &eigvals[0], &work[0],
         &lwork, &info);

  assert(info == 0);  // if it doesn't =0, there was an error.
  return {eigvals, eigvecs};
}

// Explicit specialization of norm for doubles
template <>
double Matrix<double>::norm() {
  char norm = 'F';  // tells lapack to give us Frobenius norm
  int nr = nRows;
  int nc = nCols;
  // TODO: should we allocate *work?
  std::vector<double> work(nRows);
  return dlange_(&norm, &nr, &nc, this->mat, &nr, &work[0]);
}

// Explicit specialization of norm for complex doubles
template <>
double Matrix<std::complex<double>>::norm() {
  char norm = 'F';  // tells lapack to give us Frobenius norm
  int nr = nRows;
  int nc = nCols;
  std::vector<double> work(nRows);
  return zlange_(&norm, &nr, &nc, this->mat, &nr, &work[0]);
}
