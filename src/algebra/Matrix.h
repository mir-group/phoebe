#ifndef MATRIX_H
#define MATRIX_H

// include statements
#include <assert.h>
#include <tuple>
#include <vector>
#include <cmath>
#include <complex>
#include <iostream>
#include <type_traits>

#include "Blas.h"

/** Matrix parent class, which can be used to define matrix classes of different
 *  types
 * brief General templated matrix class, with explicit specialization for
 * double and complex<double> types.
 */
template <typename T>
class Matrix {
  /// Class variables
  int nRows;
  int nCols;
  int numElements_;

  T* mat = nullptr;  // pointer to the internal array structure.

  /// Index from a 1D array to a position in a 2D array (matrix)
  long global2Local(const long& row, const long& col);
  std::tuple<long, long> local2Global(const long& k);

 public:
  static const char transN = 'N';  // no transpose nor adjoint
  static const char transT = 'T';  // transpose
  static const char transC = 'C';  // adjoint (for complex numbers)

  /** Matrix class constructor.
   * numBlocks* (ignored) are put for compatibility with ParallelMatrix.
   */
  Matrix(const int& numRows, const int& numCols, const int& numBlocksRows = 0,
         const int& numBlocksCols = 0);

  /** Default constructor
   */
  Matrix();

  /** Destructor, to delete raw buffer
   */
  ~Matrix();

  /** Copy constructor
   */
  Matrix(const Matrix<T>& that);

  /** Copy constructor
   */
  Matrix<T>& operator=(const Matrix<T>& that);

  /** Find the global indices of the matrix elements that are stored locally
   * by the current MPI process.
   */
  std::vector<std::tuple<long, long>> getAllLocalStates();

  /** Returns true if the global indices (row,col) identify a matrix element
   * stored by the MPI process.
   */
  bool indecesAreLocal(const int& row, const int& col);

  /** Find global number of rows
   */
  long rows() const;
  /** Return local number of rows */
  long localRows() const;
  /** Find global number of columns
   */
  long cols() const;
  /** Return local number of rows */
  long localCols() const;
  /** Find global number of matrix elements
   */
  long size() const;
  /** Returns a pointer to the raw matrix data buffer */ 
  T* data() const;

  /** Get and set operator
   */
  T& operator()(const int row, const int col);

  /** Const get and set operator
   */
  const T& operator()(const int row, const int col) const;

  /** Matrix-matrix multiplication.
   */
  Matrix<T> prod(const Matrix<T>& that, const char& trans1 = transN,
                 const char& trans2 = transN);

  /** Matrix-matrix addition.
   */
  Matrix<T> operator+=(const Matrix<T>& m1) {
    assert((m1.rows() == nRows) && (m1.cols() == nCols));
    for (int s = 0; s < size(); s++) mat[s] += m1.mat[s];
    return *this;
  }

  /** Matrix-matrix subtraction.
   */
  Matrix<T> operator-=(const Matrix<T>& m1) {
    assert((m1.rows() == nRows) && (m1.cols() == nCols));
    for (int s = 0; s < size(); s++) mat[s] -= m1.mat[s];
    return *this;
  }

  /** Matrix-scalar multiplication.
   */
  Matrix<T> operator*=(const T& that) {
    for (int s = 0; s < size(); s++) mat[s] *= that;
    return *this;
  }

  /** Matrix-scalar division.
   */
  Matrix<T> operator/=(const T& that) {
    for (int s = 0; s < size(); s++) mat[s] /= that;
    return *this;
  }

  /** Sets this matrix as the identity.
   * Deletes any previous content.
   */
  void eye();

  /** Diagonalize a complex-hermitian / real-symmetric matrix.
   * Nota bene: we don't check if it's hermitian/symmetric.
   */
  std::tuple<std::vector<double>, Matrix<T>> diagonalize();

  /** Computes the squared Frobenius norm of the matrix
   * (or Euclidean norm, or L2 norm of the matrix)
   */
  double squaredNorm();

  /** Computes the Frobenius norm of the matrix
   * (or Euclidean norm, or L2 norm of the matrix)
   */
  double norm();

  /** Computes a "scalar product" between two matrices A and B,
   * defined as \sum_ij A_ij * B_ij. For vectors, this reduces to the standard
   * scalar product.
   */
  T dot(const Matrix<T>& that);

  /** Unary negation
   */
  Matrix<T> operator-() const;
};

// A default constructor to build a dense matrix of zeros to be filled
template <typename T>
Matrix<T>::Matrix(const int& numRows, const int& numCols,
                  const int& numBlocksRows, const int& numBlocksCols) {
  (void) numBlocksRows;
  (void) numBlocksCols;
  nRows = numRows;
  nCols = numCols;
  numElements_ = nRows * nCols;
  mat = new T[nRows * nCols];
  for (int i = 0; i < numElements_; i++) mat[i] = 0;  // fill with zeroes
  assert(mat != nullptr);  // Memory could not be allocated, end program
}

// default constructor
template <typename T>
Matrix<T>::Matrix() {
  mat = nullptr;
  nRows = 0;
  nCols = 0;
  numElements_ = 0;
}

// copy constructor
template <typename T>
Matrix<T>::Matrix(const Matrix<T>& that) {
  nRows = that.nRows;
  nCols = that.nCols;
  numElements_ = that.numElements_;
  if (mat != nullptr) {
    delete[] mat;
    mat = nullptr;
  }
  mat = new T[numElements_];
  assert(mat != nullptr);
  for (long i = 0; i < numElements_; i++) {
    mat[i] = that.mat[i];
  }
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& that) {
  if (this != &that) {
    nRows = that.nRows;
    nCols = that.nCols;
    numElements_ = that.numElements_;
    // matrix allocation
    if (mat != nullptr) {
      delete[] mat;
    }
    mat = new T[numElements_];
    assert(mat != nullptr);
    for (long i = 0; i < numElements_; i++) {
      mat[i] = that.mat[i];
    }
  }
  return *this;
}

// destructor
template <typename T>
Matrix<T>::~Matrix() {
  delete[] mat;
}

/* ------------- Very basic operations -------------- */
template <typename T>
long Matrix<T>::rows() const {
  return nRows;
}
template <typename T>
long Matrix<T>::localRows() const {
  return nRows;
}
template <typename T>
long Matrix<T>::cols() const {
  return nCols;
}
template <typename T>
long Matrix<T>::localCols() const {
  return nCols;
}
template <typename T>
long Matrix<T>::size() const {
  return numElements_;
}
template <typename T>
T* Matrix<T>::data() const{ 
  return mat; 
}

// Get/set element
template <typename T>
T& Matrix<T>::operator()(const int row, const int col) {
  return mat[global2Local(row, col)];
}

template <typename T>
const T& Matrix<T>::operator()(const int row, const int col) const {
  return mat[global2Local(row, col)];
}

template <typename T>
bool Matrix<T>::indecesAreLocal(const int& row, const int& col) {
  (void) row;
  (void) col;
  return true;
}

template <typename T>
std::tuple<long, long> Matrix<T>::local2Global(const long& k) {
  // we convert this combined local index k into row / col indeces
  // k = j * nRows + i
  int j = k / nRows;
  int i = k - j * nRows;
  return {i, j};
}

// Indexing to set up the matrix in col major format
template <typename T>
long Matrix<T>::global2Local(const long& row, const long& col) {
  return nRows * col + row;
}

template <typename T>
std::vector<std::tuple<long, long>> Matrix<T>::getAllLocalStates() {
  std::vector<std::tuple<long, long>> x;
  for (long k = 0; k < numElements_; k++) {
    std::tuple<long, long> t = local2Global(k);  // bloch indices
    x.push_back(t);
  }
  return x;
}

// General unary negation
template <typename T>
Matrix<T> Matrix<T>::operator-() const {
  Matrix<T> c(nRows, nCols);
  for (int row = 0; row < nRows; row++) {
    for (int col = 0; col < nCols; col++) c(row, col) = -(*this)(row, col);
  }
  return c;
}

// Sets the matrix to the idenity matrix
template <typename T>
void Matrix<T>::eye() {
  assert(nRows == nCols);
  for (int row = 0; row < nRows; row++) (*this)(row, row) = (T)1.0;
}

// General function for the norm
template <typename T>
double Matrix<T>::norm() {
  T sumSq = 0;
  for (int row = 0; row < nRows; row++) {
    for (int col = 0; col < nCols; col++) {
      sumSq += ((*this)(row, col) * (*this)(row, col));
    }
  }
  return sqrt(sumSq);
}

template <typename T>
double Matrix<T>::squaredNorm() {
  double x = norm();
  return x * x;
}

template <typename T>
T Matrix<T>::dot(const Matrix<T>& that) {
  T scalar = (T)0.;
  for (int i = 0; i < numElements_; i++) {
    scalar += (*(mat + i)) * (*(that.mat + i));
  }
  return scalar;
}

#endif  // MATRIX_H
