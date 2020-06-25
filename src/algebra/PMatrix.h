#ifdef MPI_AVAIL

#ifndef PMATRIX_H
#define PMATRIX_H

#include <tuple>
#include <vector>
#include "bandstructure.h"

// double matrix (I skip the templates for now)
template <typename T>
class Matrix {  // P stands for parallel
 protected:
  /// Class variables
  int numRows_, numCols_;
  int numLocalRows_, numLocalCols_;
  int numLocalElements_;
  int numBlocksRows_, numBlocksCols_;
  int blockSizeRows_, blockSizeCols_;
  int descMat_[9];
  int numBlasRows_, numBlasCols_;
  int myBlasRow_, myBlasCol_;

  //  char transposition = 'N'; // status of the transposition
  const char transN = 'N';  // no transpose nor adjoint
  const char transT = 'T';  // transpose
  const char transC = 'C';  // adjoint (for complex numbers)

  std::tuple<long, long> local2Global(const int& k);

  // dummy values to return when accessing elements not available locally
  T dummyZero = 0;
  T const dummyConstZero = 0;

 public:
  T* mat = nullptr;

  // Construct using row, col numbers, and block distribution
  Matrix(const int& numRows, const int& numCols, const int& numBlocksRows = 0,
         const int& numBlocksCols = 0);
  Matrix(); // default constructor
  ~Matrix(); // deallocate pointers

  /** Copy constructor
   */
  Matrix(const Matrix<T>& that);

  /** Copy assignment
   */
  Matrix& operator=(const Matrix<T>& that);

  /** Find all the wavevector pairs (iq1,iq2) that should be computed by the
   * local MPI process. This method is specifically made for the scattering
   * matrix, which has rows spanned by Bloch states (iq,ib)
   */
  std::vector<std::tuple<long, long>> getAllLocalWavevectors(
      BaseBandStructure& bandStructure);

  /** Find the global indices of the matrix elements that are stored locally
   * by the current MPI process.
   */
  std::vector<std::tuple<long, long>> getAllLocalStates();

  // utilities for the global matrix size
  /** Find global number of rows
   */
  long rows() const;
  /** Find global number of columns
   */
  long cols() const;
  /** Find global number of matrix elements
   */
  long size() const;

  /** Get and set operator
   */
  T& operator()(const int row, const int col);

  /** Const get and set operator
   */
  const T& operator()(const int row, const int col) const;

  /** Matrix-matrix multiplication.
   */
  Matrix<T> prod(const Matrix<T>& that, const char& trans1, const char& trans2);
  /** Matrix-scalar multiplication.
   */
  Matrix<T> operator*=(const T& that);
  /** Matrix-scalar division.
   */
  Matrix<T> operator/=(const T& that);

  /** Matrix-matrix addition.
   */
  Matrix<T> operator+=(const Matrix<T>& that);

  /** Matrix-matrix subtraction.
   */
  Matrix<T> operator-=(const Matrix<T>& that);

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
  T squaredNorm();

  /** Computes the Frobenius norm of the matrix
   * (or Euclidean norm, or L2 norm of the matrix)
   */
  T norm();

  /** Computes a "scalar product" between two matrices A and B,
   * defined as \sum_ij A_ij * B_ij. For vectors, this reduces to the standard
   * scalar product.
   */
  T dot(const Matrix<T>& that);

  /** Unary negation
   */
  Matrix<T> operator-() const;
};

#include "PMatrix.cpp"

#endif // include safeguard

#endif // mpi_avail
