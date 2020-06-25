#ifndef PMATRIX_H
#define PMATRIX_H

#include <tuple>
#include <vector>
#include "bandstructure.h"

// double matrix (I skip the templates for now)
class PMatrix {  // P stands for parallel
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

  double dummyZero = 0.;
  double const dummyConstZero = 0.;

 public:
  double* mat = nullptr;

  PMatrix(const int& numRows, const int& numCols, const int& numBlocksRows = 0,
          const int& numBlocksCols = 0);  // Construct using row, col numbers
  PMatrix();                              // default constructor
  ~PMatrix();

  PMatrix(const PMatrix& that);            /// Copy constructor
  PMatrix& operator=(const PMatrix& that);  /// Copy constructor

  long rows() const;
  long cols() const;
  long size() const;

  // Get and set operators
  double& operator()(const int row, const int col);
  const double& operator()(const int row, const int col) const;

  // Generic matrix multiplication.
  //  PMatrix& operator*=(const PMatrix& that);
  PMatrix prod(const PMatrix& that, const char& trans1, const char& trans2);
  PMatrix operator*=(const double& that);
  PMatrix operator/=(const double& that);

  PMatrix operator+=(const PMatrix& that);
  PMatrix operator-=(const PMatrix& that);

  /** Sets this matrix as the identity
   */
  void eye();

  /** Diagonalize a matrix
   */
  std::tuple<std::vector<double>, PMatrix> diagonalize();

  std::vector<std::tuple<long, long>> getAllLocalWavevectors(
      BaseBandStructure& bandStructure);

  std::vector<std::tuple<long, long>> getAllLocalStates();

  // Matrix multiplication.
  PMatrix& operator*=(const PMatrix& that);

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
  double dot(const PMatrix& that);

  /** Unary negation
   */
  PMatrix operator-() const;
};

#endif
