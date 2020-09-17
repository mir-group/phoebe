#include "Matrix.h"

// Explict specialization of BLAS matrix-matrix mult for Matrix<complex<double>>
template <>
Matrix<std::complex<double>> Matrix<std::complex<double>>::prod(
    const Matrix<std::complex<double>>& that, const char& trans1,
    const char& trans2) {
  Matrix<std::complex<double>> c(*this); // copy this matrix
  if(isDistributed) *(c.pmat) = pmat->prod( *(that.pmat),trans1,trans2);
  else{ *(c.mat) = mat->prod( *(that.mat),trans1,trans2); }
  return c;
}

// Explicit specializiation of BLAS matrix-matrix mult for Matrix<double>
template <>
Matrix<double> Matrix<double>::prod(const Matrix<double>& that,
                                    const char& trans1, const char& trans2) {
  Matrix<double> c(*this); // copy this matrix
  if(isDistributed) *(c.pmat) = pmat->prod( *(that.pmat),trans1,trans2);
  else{ *(c.mat) = mat->prod( *(that.mat),trans1,trans2); }
  return c;
}

// Diagonalize a complex double hermitian matrix
template <>
std::tuple<std::vector<double>, Matrix<std::complex<double>>>
Matrix<std::complex<double>>::diagonalize() {

  std::vector<double> eigvals;
  Matrix<std::complex<double>> eigvecs(*this);

  if(isDistributed) {
    auto tup = pmat->diagonalize();
    eigvals = std::get<0>(tup);
    eigvecs.pmat = &(std::get<1>(tup));
  }
  else{
    auto tup = mat->diagonalize();
    eigvals = std::get<0>(tup);
    eigvecs.mat = &(std::get<1>(tup));
  }
  return {eigvals,eigvecs};
}

// Diagonalize for real double symmetric matrix
template <>
std::tuple<std::vector<double>, Matrix<double>> Matrix<double>::diagonalize() {

  std::vector<double> eigvals;
  Matrix<double> eigvecs(*this);

  if(isDistributed) {
    auto tup = pmat->diagonalize();
    eigvals = std::get<0>(tup);
    eigvecs.pmat = &(std::get<1>(tup)); // returns a pmat, need the pointer to it
  }
  else{
    auto tup = mat->diagonalize();
    eigvals = std::get<0>(tup);
    eigvecs.mat = &(std::get<1>(tup));
  }
  return {eigvals,eigvecs};
}

// Explicit specialization of norm for doubles
template <>
double Matrix<double>::norm() {
  if(isDistributed) { return pmat->norm(); }
  else{ return mat->norm(); }
}

// Explicit specialization of norm for complex doubles
template <>
double Matrix<std::complex<double>>::norm() {
  if(isDistributed) {
    Error e(" No implemented specialization of PMatrix::complex<double>.norm()");
    return 0;
  }
  else{ return mat->norm(); }
}
