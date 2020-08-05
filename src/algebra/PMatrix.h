#ifndef PMATRIX_H
#define PMATRIX_H

#include <tuple>
#include <vector>

#include "Blacs.h"
#include "Matrix.h"
#include "bandstructure.h"
#include "constants.h"
#include "mpiHelper.h"

#ifndef MPI_AVAIL

// alias template
template<typename T>
using ParallelMatrix = Matrix<T>;

#else

#include <utility>
#include <set>

/** Class for managing a matrix MPI-distributed in memory.
 *
 * This class uses the Scalapack library for matrix-matrix multiplication and
 * matrix diagonalization. For the time being we don't use other scalapack
 * functionalities.
 *
 * If the code is compiled without MPI, ParallelMatrix reduces to its serial
 * version Matrix.
 *
 * Template specialization only valid for double or complex<double>.
 */
template <typename T>
class ParallelMatrix {
 private:
  /// Class variables
  int numRows_, numCols_;
  int numLocalRows_, numLocalCols_;
  int numLocalElements_;
  int numBlocksRows_, numBlocksCols_;
  int blockSizeRows_, blockSizeCols_;
  int descMat_[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  int numBlasRows_, numBlasCols_;
  int myBlasRow_, myBlasCol_;

  // dummy values to return when accessing elements not available locally
  T dummyZero = 0;
  T const dummyConstZero = 0;

  T* mat = nullptr; // raw buffer

  /** Converts a local one-dimensional storage index (MPI-dependent) into the
   * row/column index of the global matrix.
   */
  std::tuple<long, long> local2Global(const long& k);

  /** Converts a global row/column index of the global matrix into a local
   * one-dimensional storage index (MPI-dependent),
   * with value ranging from  0 to numLocalElements_-1.
   * Returns -1 if the matrix element is not stored on the current MPI process.
   */
  long global2Local(const long& row, const long& col);

 public:
  static const char transN = 'N';  // no transpose nor adjoint
  static const char transT = 'T';  // transpose
  static const char transC = 'C';  // adjoint (for complex numbers)

  /** Default constructor of the matrix class.
   * Matrix elements are set to zero in the initialization.
   * @param numRows: global number of matrix rows
   * @param numCols: global number of matrix columns
   * @param numBLocksRows: row size of the block for Blacs distribution
   * @param numBLocksCols: column size of the block for Blacs distribution
   */
  ParallelMatrix(const int& numRows, const int& numCols,
                 const int& numBlocksRows = 0, const int& numBlocksCols = 0);

  /** Empty constructor
   */
  ParallelMatrix();

  /** Destructor, to deallocate pointers
   */
  ~ParallelMatrix();

  /** Copy constructor
   */
  ParallelMatrix(const ParallelMatrix<T>& that);

  /** Copy assignment
   */
  ParallelMatrix& operator=(const ParallelMatrix<T>& that);

  /** Find all the wavevector pairs (iq1,iq2) that should be computed by the
   * local MPI process. This method is specifically made for the scattering
   * matrix, which has rows spanned by Bloch states (iq,ib)
   */
  std::vector<std::pair<int, int>> getAllLocalWavevectors(
      BaseBandStructure& bandStructure);

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
  /** Find global number of columns
   */
  long cols() const;
  /** Find global number of matrix elements
   */
  long size() const;

  /** Get and set operator.
   * Returns the stored value if the matrix element (row,col) is stored in
   * memory by the MPI process, otherwise returns zero.
   */
  T& operator()(const int row, const int col);

  /** Const get and set operator
   * Returns the stored value if the matrix element (row,col) is stored in
   * memory by the MPI process, otherwise returns zero.
   */
  const T& operator()(const int row, const int col) const;

  /** Matrix-matrix multiplication.
   * Computes result = trans1(*this) * trans2(that)
   * where trans(1/2( can be "N" (matrix as is), "T" (transpose) or "C" adjoint
   * (these flags are used so that the transposition/adjoint operation is never
   * directly operated on the stored values)
   *
   * @param that: the matrix to multiply "this" with
   * @param trans2: controls transposition of "this" matrix
   * @param trans1: controls transposition of "that" matrix
   * @return result: a ParallelMatrix object.
   */
  ParallelMatrix<T> prod(const ParallelMatrix<T>& that,
                         const char& trans1 = transN,
                         const char& trans2 = transN);

  /** Matrix-matrix addition.
   */
  ParallelMatrix<T>& operator+=(const ParallelMatrix<T>& that);

  /** Matrix-matrix subtraction.
   */
  ParallelMatrix<T>& operator-=(const ParallelMatrix<T>& that);

  /** Matrix-scalar multiplication.
   */
  ParallelMatrix<T>& operator*=(const T& that);

  /** Matrix-scalar division.
   */
  ParallelMatrix<T>& operator/=(const T& that);

  /** Sets this matrix as the identity.
   * Deletes any previous content.
   */
  void eye();

  /** Diagonalize a complex-hermitian or real-symmetric matrix.
   * Nota bene: we don't check if the matrix is hermitian/symmetric or not.
   * By default, it operates on the upper-triangular part of the matrix.
   */
  std::tuple<std::vector<double>, ParallelMatrix<T>> diagonalize();

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
  T dot(const ParallelMatrix<T>& that);

  /** Unary negation
   */
  ParallelMatrix<T> operator-() const;
};

template <typename T>
ParallelMatrix<T>::ParallelMatrix(const int& numRows, const int& numCols,
                                  const int& numBlocksRows,
                                  const int& numBlocksCols) {

  // if blacs is not initalized, we need to start it. 
  mpi->initBlacs(); 
  
  // initialize number of rows and columns of the global matrix
  numRows_ = numRows;
  numCols_ = numCols;

  numBlasRows_ = mpi->getNumBlasRows();
  numBlasCols_ = mpi->getNumBlasCols();
  myBlasRow_ = mpi->getMyBlasRow();
  myBlasCol_ = mpi->getMyBlasCol();

  // determine the number of blocks (for parallel distribution) along rows/cols
  // by default, we set the number of blocks equal to the number of
  // rows in the blacs grid of processes
  if (numBlocksRows != 0) {
    numBlocksRows_ = numBlocksRows;
  } else {
    numBlocksRows_ = numBlasRows_;
  }
  if (numBlocksCols != 0) {
    numBlocksCols_ = numBlocksCols;
  } else {
    numBlocksCols_ = numBlasCols_;
  }

  // compute the block size (over which matrix is distributed)
  blockSizeRows_ = numRows_ / numBlocksRows_;
  if (numRows_ % numBlocksRows_ != 0) blockSizeRows_ += 1;
  blockSizeCols_ = numCols_ / numBlocksCols_;
  if (numCols_ % numBlocksCols_ != 0) blockSizeCols_ += 1;

  // determine the number of local rows and columns
  int iZero = 0;  // helper variable
  numLocalRows_ =
      numroc_(&numRows_, &blockSizeRows_, &myBlasRow_, &iZero, &numBlasRows_);
  numLocalCols_ =
      numroc_(&numCols_, &blockSizeCols_, &myBlasCol_, &iZero, &numBlasCols_);
  numLocalElements_ = numLocalRows_ * numLocalCols_;
  mat = new T[numLocalElements_];

  // Memory could not be allocated, end program
  assert(mat != nullptr);

  // fill the matrix with zeroes
  for (long i = 0; i < numLocalElements_; ++i) *(mat + i) = 0.;  // mat[i] = 0.;

  // Create descriptor for block cyclic distribution of matrix
  int info;  // error code
  int lddA =
      numLocalRows_ > 1 ? numLocalRows_ : 1;  // if mpA>1, ldda=mpA, else 1
  int blacsContext = mpi->getBlacsContext();
  descinit_(descMat_, &numRows_, &numCols_, &blockSizeRows_, &blockSizeCols_,
            &iZero, &iZero, &blacsContext, &lddA, &info);
  if (info != 0) {
    Error e("Something wrong calling descinit", info);
  }
}

template <typename T>
ParallelMatrix<T>::ParallelMatrix() {
  numRows_ = 0;
  numCols_ = 0;
  numLocalRows_ = 0;
  numLocalCols_ = 0;
  numLocalElements_ = 0;
  numBlocksRows_ = 0;
  numBlocksCols_ = 0;
  blockSizeRows_ = 0;
  blockSizeCols_ = 0;
  numBlasRows_ = 0;
  numBlasCols_ = 0;
  myBlasRow_ = 0;
  myBlasCol_ = 0;
  mpi->initBlacs(); 
}

template <typename T>
ParallelMatrix<T>::ParallelMatrix(const ParallelMatrix<T>& that) {
  numRows_ = that.numRows_;
  numCols_ = that.numCols_;
  numLocalRows_ = that.numLocalRows_;
  numLocalCols_ = that.numLocalCols_;
  numLocalElements_ = that.numLocalElements_;
  numBlocksRows_ = that.numBlocksRows_;
  numBlocksCols_ = that.numBlocksCols_;
  blockSizeRows_ = that.blockSizeRows_;
  blockSizeCols_ = that.blockSizeCols_;
  numBlasRows_ = that.numBlasRows_;
  numBlasCols_ = that.numBlasCols_;
  myBlasRow_ = that.myBlasRow_;
  myBlasCol_ = that.myBlasCol_;

  for (int i = 0; i < 9; i++) {
    descMat_[i] = that.descMat_[i];
  }
  if ( mat != nullptr ) {
    delete[] mat;
  }
  // matrix allocation
  mat = new T[numLocalElements_];
  // Memory could not be allocated, end program
  assert(mat != nullptr);
  for (long i = 0; i < numLocalElements_; i++) {
    mat[i] = that.mat[i];
  }
}

template <typename T>
ParallelMatrix<T>& ParallelMatrix<T>::operator=(const ParallelMatrix<T>& that) {
  if (this != &that) {
    numRows_ = that.numRows_;
    numCols_ = that.numCols_;
    numLocalRows_ = that.numLocalRows_;
    numLocalCols_ = that.numLocalCols_;
    numLocalElements_ = that.numLocalElements_;
    numBlocksRows_ = that.numBlocksRows_;
    numBlocksCols_ = that.numBlocksCols_;
    blockSizeRows_ = that.blockSizeRows_;
    blockSizeCols_ = that.blockSizeCols_;
    numBlasRows_ = that.numBlasRows_;
    numBlasCols_ = that.numBlasCols_;
    myBlasRow_ = that.myBlasRow_;
    myBlasCol_ = that.myBlasCol_;

    for (int i = 0; i < 9; i++) {
      descMat_[i] = that.descMat_[i];
    }
    if ( mat != nullptr ) {
      delete[] mat;
    }
    // matrix allocation
    mat = new T[numLocalElements_];
    // Memory could not be allocated, end program
    assert(mat != nullptr);
    for (long i = 0; i < numLocalElements_; i++) {
      mat[i] = that.mat[i];
    }
  }
  return *this;
}

template <typename T>
ParallelMatrix<T>::~ParallelMatrix() {
  delete[] mat;
}

template <typename T>
long ParallelMatrix<T>::rows() const {
  return numRows_;
}

template <typename T>
long ParallelMatrix<T>::cols() const {
  return numCols_;
}

template <typename T>
long ParallelMatrix<T>::size() const {
  return cols() * rows();
}

// Get/set element

template <typename T>
T& ParallelMatrix<T>::operator()(const int row, const int col) {
  long localIndex = global2Local(row, col);
  if (localIndex == -1) {
    dummyZero = 0.;
    return dummyZero;
  } else {
    return mat[localIndex];
  }
}

template <typename T>
const T& ParallelMatrix<T>::operator()(const int row, const int col) const {
  long localIndex = global2Local(row, col);
  if (localIndex == -1) {
    return dummyConstZero;
  } else {
    return mat[localIndex];
  }
}

template <typename T>
bool ParallelMatrix<T>::indecesAreLocal(const int& row, const int& col) {
  long localIndex = global2Local(row, col);
  if (localIndex == -1) {
    return false;
  } else {
    return true;
  }
}

template <typename T>
std::tuple<long, long> ParallelMatrix<T>::local2Global(const long& k) {
  // first, we convert this combined local index k
  // into local row / col indeces
  // k = j * numLocalRows_ + i
  int j = k / numLocalRows_;
  int i = k - j * numLocalRows_;

  // now we can convert local row/col indices into global indices

  int l_i = i / blockSizeRows_;  // which block
  int x_i = i % blockSizeRows_;  // where within that block
  // global row
  int I = (l_i * numBlasRows_ + myBlasRow_) * blockSizeRows_ + x_i;

  int l_j = j / blockSizeCols_;  // which block
  int x_j = j % blockSizeCols_;  // where within that block
  // global col
  int J = (l_j * numBlasCols_ + myBlasCol_) * blockSizeCols_ + x_j;
  return {I, J};
}

template <typename T>
long ParallelMatrix<T>::global2Local(const long& row, const long& col) {
  // note: row and col indices use the c++ convention of running from 0 to N-1
  // fortran (infog2l_) wants indices from 1 to N.
  int row_ = row + 1;
  int col_ = col + 1;

  // first, we find the local index
  int iia, jja, iarow, iacol;
  infog2l_(&row_, &col_, &descMat_[0], &numBlasRows_, &numBlasCols_,
           &myBlasRow_, &myBlasCol_, &iia, &jja, &iarow, &iacol);
  if (myBlasRow_ == iarow && myBlasCol_ == iacol) {
    return iia + (jja - 1) * descMat_[8] - 1;
  } else {
    return -1;
  }
}

template <typename T>
std::vector<std::tuple<long, long>> ParallelMatrix<T>::getAllLocalStates() {
  std::vector<std::tuple<long, long>> x;
  for (long k = 0; k < numLocalElements_; k++) {
    std::tuple<long, long> t = local2Global(k);  // bloch indices
    x.push_back(t);
  }
  return x;
}

template <typename T>
std::vector<std::pair<int, int>> ParallelMatrix<T>::getAllLocalWavevectors(
    BaseBandStructure& bandStructure) {

  std::set<std::pair<int,int>> x;
  for (long k = 0; k < numLocalElements_; k++) {
    auto tup = local2Global(k);
    auto is1 = std::get<0>(tup);
    auto is2 = std::get<1>(tup);  // bloch indices
    auto tup1 = bandStructure.getIndex(is1);
    auto ik1 = std::get<0>(tup1);
    auto ib1 = std::get<1>(tup1);
    auto tup2 = bandStructure.getIndex(is2);
    auto ik2 = std::get<0>(tup2);
    auto ib2 = std::get<1>(tup2);
    std::pair<int,int> xx = std::make_pair(ik1.get(), ik2.get());
    x.insert(xx);
  }
  std::vector<std::pair<int, int>> wavevectorPairs(x.size());
  for ( auto t : x ) {
    wavevectorPairs.push_back(t);
  }
  return wavevectorPairs;
}

template <typename T>
ParallelMatrix<T>& ParallelMatrix<T>::operator*=(const T& that) {
  for (long i = 0; i < numLocalElements_; i++) {
    *(mat + i) *= that;
  }
  return *this;
}

template <typename T>
ParallelMatrix<T>& ParallelMatrix<T>::operator/=(const T& that) {
  for (long i = 0; i < numLocalElements_; i++) {
    *(mat + i) /= that;
  }
  return *this;
}

template <typename T>
ParallelMatrix<T>& ParallelMatrix<T>::operator+=(const ParallelMatrix<T>& that) {
  for (long i = 0; i < numLocalElements_; i++) {
    *(mat + i) += *(that.mat + i);
  }
  return *this;
}

template <typename T>
ParallelMatrix<T>& ParallelMatrix<T>::operator-=(const ParallelMatrix<T>& that) {
  for (long i = 0; i < numLocalElements_; i++) {
    *(mat + i) -= *(that.mat + i);
  }
  return *this;
}

template <typename T>
void ParallelMatrix<T>::eye() {
  if (numRows_ != numCols_) {
    Error e("Cannot build an identity matrix with non-square matrix");
  }
  for (int i = 0; i < numLocalElements_; i++) {
    *(mat + i) = 0.;
  }
  for (int i = 0; i < numRows_; i++) {
    operator()(i, i) = 1.;
  }
}

template <typename T>
T ParallelMatrix<T>::squaredNorm() {
  return dot(*this);
}

template <typename T>
T ParallelMatrix<T>::norm() {
  return sqrt(squaredNorm());
}

template <typename T>
T ParallelMatrix<T>::dot(const ParallelMatrix<T>& that) {
  T scalar = dummyConstZero;
  T scalarOut = dummyConstZero;
  for (int i = 0; i < numLocalElements_; i++) {
    scalar += (*(mat + i)) * (*(that.mat + i));
  }
  mpi->allReduceSum(&scalar, &scalarOut);
  return scalarOut;
}

template <typename T>
ParallelMatrix<T> ParallelMatrix<T>::operator-() const {
  ParallelMatrix<T> result = *this;
  for (long i = 0; i < numLocalElements_; i++) {
    *(result.mat + i) = -*(result.mat + i);
  }
  return result;
}

#endif  // include safeguard

#endif  // mpi_avail
