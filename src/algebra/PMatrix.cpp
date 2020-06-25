#ifdef MPI_AVAIL

#include "PMatrix.h"
#include "Blacs.h"
#include "mpiHelper.h"
#include "constants.h"

template <typename T>
Matrix<T>::Matrix(const int& numRows, const int& numCols,
                 const int& numBlocksRows, const int& numBlocksCols) {

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
  if ( numRows_ % numBlocksRows_ != 0 ) blockSizeRows_ += 1;
  blockSizeCols_ = numCols_ / numBlocksCols_;
  if ( numCols_ % numBlocksCols_ != 0 ) blockSizeCols_ += 1;

  // determine the number of local rows and columns
  int iZero = 0; // helper variable
  int numLocalRows_ = numroc_(&numRows_, &blockSizeRows_, &myBlasRow_, &iZero,
          &numBlasRows_);
  int numLocalCols_ = numroc_(&numCols_, &blockSizeCols_, &myBlasCol_, &iZero,
          &numBlasCols_);
  numLocalElements_ = numLocalRows_ * numLocalCols_;

  mat = new T[numLocalElements_];

  // Memory could not be allocated, end program
  assert(mat != nullptr);

  // fill the matrix with zeroes
  for (long i = 0; i < numLocalElements_; ++i) *(mat+i) = 0.; // mat[i] = 0.;

  // Create descriptor
  int descMat_[9]; // descriptor
  int info; // error code
  int lddA = numLocalRows_ > 1 ? numLocalRows_ : 1; // if mpA>1, ldda=mpA, else 1
  int blacsContext = mpi->getBlacsContext();
  descinit_( descMat_,  &numRows_, &numCols_, &blockSizeRows_, &blockSizeCols_,
          &iZero, &iZero, &blacsContext, &lddA, &info);
  if (info != 0) {
    Error e("Something wrong calling descinit", info);
  }
}

template <typename T>
Matrix<T>::Matrix() {
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
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T>& that) {
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
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& that) {
  if ( this != &that) {
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
  }
  return *this;
}

template <typename T>
Matrix<T>::~Matrix() { delete[] mat; }

template <typename T>
long Matrix<T>::rows() const { return numRows_; }

template <typename T>
long Matrix<T>::cols() const { return numCols_; }

template <typename T>
long Matrix<T>::size() const { return cols()*rows(); }

// Get element

template <typename T>
T& Matrix<T>::operator()(const int row, const int col) {

  // first, we find the local index
  int localIndex;
  int iia, jja, iarow, iacol;
  infog2l_(&row, &col, &descMat_[0], &numBlasRows_, &numBlasCols_, &myBlasRow_,
           &myBlasCol_, &iia, &jja, &iarow, &iacol);
  if (myBlasRow_ == iarow && myBlasCol_ == iacol) {
    localIndex = iia + (jja - 1) * descMat_[8] - 1;
  } else {
    localIndex = -1;
  }

  if (localIndex == -1) {
    dummyZero = 0.;
    return dummyZero;
  } else {
    return mat[localIndex];
  }
}

template <typename T>
const T& Matrix<T>::operator()(const int row, const int col) const {

  // first, we find the local index
  int localIndex;
  int iia, jja, iarow, iacol;
  infog2l_(&row, &col, &descMat_[0], &numBlasRows_, &numBlasCols_, &myBlasRow_,
           &myBlasCol_, &iia, &jja, &iarow, &iacol);
  if (myBlasRow_ == iarow && myBlasCol_ == iacol) {
    localIndex = iia + (jja - 1) * descMat_[8] - 1;
  } else {
    localIndex = -1;
  }

  if (localIndex == -1) {
    return dummyConstZero;
  } else {
    return mat[localIndex];
  }
}

template <typename T>
std::tuple<long,long> Matrix<T>::local2Global(const int &k) {
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
std::vector<std::tuple<long, long>> Matrix<T>::getAllLocalWavevectors(
        BaseBandStructure &bandStructure) {
  std::vector<std::tuple<long, long>> wavevectorPairs;
  for (long k = 0; k < numLocalElements_; k++) {
    auto [is1, is2] = local2Global(k);  // bloch indices
    auto [ik1, ib1] = bandStructure.getIndex(is1);
    auto [ik2, ib2] = bandStructure.getIndex(is2);
    // make a pair of these wavevectors
    auto t = std::make_tuple(ik1.get(), ik2.get());
    // add to list if unique
    if (std::find(wavevectorPairs.begin(), wavevectorPairs.end(), t) ==
        wavevectorPairs.end()) {
      wavevectorPairs.push_back(t);
    }
  }
  return wavevectorPairs;
}

template <>
Matrix<double> Matrix<double>::prod(const Matrix<double>& that,
                                    const char& trans1, const char& trans2) {
  Matrix<double> result(numRows_, numCols_, numBlocksRows_, numBlocksCols_);

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
  pdgemm_(&transN, &transN, &m, &n, &k, &alpha, mat, &one, &one, &descMat_[0],
          that.mat, &one, &one, &that.descMat_[0], &beta, result.mat, &one,
          &one, &result.descMat_[0]);
  return result;
}

template <>
Matrix<std::complex<double>> Matrix<std::complex<double>>::prod(
    const Matrix<std::complex<double>>& that, const char& trans1,
    const char& trans2) {
  Matrix<std::complex<double>> result(numRows_, numCols_, numBlocksRows_, numBlocksCols_);

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
  pzgemm_(&transN, &transN, &m, &n, &k, &alpha, mat, &one, &one, &descMat_[0],
          that.mat, &one, &one, &that.descMat_[0], &beta, result.mat, &one,
          &one, &result.descMat_[0]);
  return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator*=(const T& that) {
    Matrix<T> result(numRows_, numCols_, numBlocksRows_, numBlocksCols_);
    for ( long i=0; i<numLocalElements_; i++ ) {
        *(result.mat+i) = *(mat+i) * that;
    }
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator/=(const T& that) {
    Matrix<T> result(numRows_, numCols_, numBlocksRows_, numBlocksCols_);
    for ( long i=0; i<numLocalElements_; i++ ) {
        *(result.mat+i) = *(mat+i) / that;
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator+=(const Matrix<T>& that) {
    Matrix result = *this;
    for ( long i=0; i<numLocalElements_; i++ ) {
        *(result.mat+i) += *(that.mat+i);
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator-=(const Matrix<T>& that) {
    Matrix<T> result = *this;
    for ( long i=0; i<numLocalElements_; i++ ) {
        *(result.mat+i) -= *(that.mat+i);
    }
    return result;
}

template<typename T>
void Matrix<T>::eye() {
    if ( numRows_ != numCols_ ) {
        Error e("Cannot build an identity matrix with non-square matrix");
    }
    for ( int  i=0; i<numLocalElements_; i++ ) {
        *(mat+i) = 0.;
    }
    for ( int  i=0; i<numRows_; i++ ) {
        operator()(i,i) = 1.;
    }
}

template <>
std::tuple<std::vector<double>, Matrix<double>> Matrix<double>::diagonalize() {
  if (numRows_ != numCols_) {
    Error e("Can not diagonalize non-square matrix");
  }
  double* eigenvalues = nullptr;
  eigenvalues = new double[numRows_];

  Matrix<double> eigenvectors(numRows_, numCols_, numBlocksRows_,
                              numBlocksCols_);

  // find the value of lwork. These are internal "scratch" arrays
  int nn = std::max(std::max(numRows_, 2), numBlocksRows_);
  int izero = 0;
  int np = numroc_(&nn, &numBlocksRows_, &izero, &izero, &numBlasRows_);
  int sizesytrd = std::max(numBlocksRows_ * (np + 1), 3 * numBlocksRows_);
  int lwork = 5 * numRows_ + sizesytrd + 1;

  // double work[lwork];
  double* work = nullptr;
  work = new double[lwork];

  char jobz = 'V';  // also eigenvectors
  char uplo = 'U';  // upper triangolar
  int ia = 1;       // row index from which we diagonalize
  int ja = 1;       // row index from which we diagonalize
  int info = 0;
  pdsyev_(&jobz, &uplo, &numRows_, mat, &ia, &ja, &descMat_[0], eigenvalues,
         eigenvectors.mat, &ia, &ja, &eigenvectors.descMat_[0], work, &lwork,
         &info);

  if (info != 0) {
    Error e("PDSYEV failed", info);
  }

  std::vector<double> eigenvalues_(numRows_);
  for ( long i=0; i<numRows_; i++ ) {
    eigenvalues_[i] = *(eigenvalues+i);
  }
  delete[] eigenvalues;
  delete[] work;
  // note that the scattering matrix now has different values

  return {eigenvalues_, eigenvectors};
}

template <>
std::tuple<std::vector<double>, Matrix<std::complex<double>>>
        Matrix<std::complex<double>>::diagonalize() {
  if (numRows_ != numCols_) {
    Error e("Can not diagonalize non-square matrix");
  }
  double* eigenvalues = nullptr;
  eigenvalues = new double[numRows_];

  Matrix<std::complex<double>> eigenvectors(numRows_, numCols_, numBlocksRows_,
                                            numBlocksCols_);

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
    Error e("PDSYEV failed", info);
  }

  std::vector<double> eigenvalues_(numRows_);
  for ( long i=0; i<numRows_; i++ ) {
    eigenvalues_[i] = *(eigenvalues+i);
  }
  delete[] eigenvalues;
  delete[] work;
  delete[] rwork;
  // note that the scattering matrix now has different values

  return {eigenvalues_, eigenvectors};
}

template <typename T>
T Matrix<T>::squaredNorm() {
    return dot(*this);
}

template <typename T>
T Matrix<T>::norm() {
    return sqrt(squaredNorm());
}

template<typename T>
T Matrix<T>::dot(const Matrix<T>& that) {
    T scalar = dummyConstZero;
    T scalarOut = dummyConstZero;
    for ( int i=0; i<numLocalElements_; i++ ) {
        scalar += ( *(mat+i) ) * ( *(that.mat+i) );
    }
    mpi->reduceSum(&scalar, &scalarOut);
    return scalarOut;
}

template<typename T>
Matrix<T> Matrix<T>::operator-() const {
    Matrix<T> result = *this;
    for ( long i=0; i<numLocalElements_; i++ ) {
        *(result.mat+i) = - *(result.mat+i);
    }
    return result;
}

#endif // MPI_AVAIL
