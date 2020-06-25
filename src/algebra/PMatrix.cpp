#include "PMatrix.h"
#include "Blacs.h"
#include "mpiHelper.h"

PMatrix::PMatrix(const int& numRows, const int& numCols,
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

  mat = new double[numLocalElements_];

  // Memory could not be allocated, end program
  assert(mat != nullptr);

  // fill the matrix with zeroes
  for (int i = 0; i < numLocalElements_; ++i) *(mat+i) = 0.; // mat[i] = 0.;

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

PMatrix::PMatrix() {
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

PMatrix::PMatrix(const PMatrix& that) {
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

PMatrix& PMatrix::operator=(const PMatrix& that) {
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

PMatrix::~PMatrix() { delete[] mat; }

long PMatrix::rows() const { return numRows_; }
long PMatrix::cols() const { return numCols_; }
long PMatrix::size() const { return cols()*rows(); }

// Get element

double& PMatrix::operator()(const int row, const int col) {

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

const double& PMatrix::operator()(const int row, const int col) const {

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

//double& PMatrix::operator()(const int row, const int col) const {
//    double value;
//    char bcast = ' '; // only the process containing mat(row,col) updates value
//    char bcastTopology = ' ';
//    pdelget_(&bcast, &bcastTopology, &value, mat, &row, &col, &descMat_[0]);
//    return value;
//}
//
//// Set element
//void PMatrix::operator()(const int row, const int col, const double value){
//    pdelset_(mat, &row, &col, &descMat_, &value);
//}

void PMatrix::eye() {
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

std::tuple<long,long> PMatrix::local2Global(const int &k) {
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

std::vector<std::tuple<long, long>> PMatrix::getAllLocalWavevectors(
        BaseBandStructure &bandStructure) {
  std::vector<std::tuple<long, long>> wavevectorPairs;
  for (int k = 0; k < numLocalElements_; k++) {
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

PMatrix PMatrix::prod(const PMatrix& that, const char & trans1,
        const char & trans2) {
    PMatrix result(numRows_, numCols_, numBlocksRows_, numBlocksCols_);

    int m;
    if ( trans1 == transN ) {
        m = numRows_;
    } else {
        m = numCols_;
    }
    int n;
    if ( trans2 == transN ) {
        n = that.numCols_;
    } else {
        n = that.numRows_;
    }
    int k;
    if ( trans1 == transN ) {
        k = numCols_;
    } else {
        k = numRows_;
    }
    // check on k being consistent
    if ( trans2 == transN ) {
        assert(k == that.numRows_);
    } else {
        assert(k == that.numCols_);
    }
    double alpha = 1.;
    double beta = 1.;
    int one = 1;
    pdgemm_(&transN, &transN, &m, &n, &k, &alpha, mat, &one, &one, &descMat_[0],
            that.mat, &one, &one, &that.descMat_[0], &beta, result.mat, &one,
            &one, &result.descMat_[0]);
    return result;
}

PMatrix PMatrix::operator*=(const double& that) {
    PMatrix result(numRows_, numCols_, numBlocksRows_, numBlocksCols_);
    for ( int  i=0; i<numLocalElements_; i++ ) {
        *(result.mat+i) = *(mat+i) * that;
    }
    return result;
}

PMatrix PMatrix::operator/=(const double& that) {
    PMatrix result(numRows_, numCols_, numBlocksRows_, numBlocksCols_);
    for ( int  i=0; i<numLocalElements_; i++ ) {
        *(result.mat+i) = *(mat+i) / that;
    }
    return result;
}

std::tuple<std::vector<double>, PMatrix> PMatrix::diagonalize() {
  if (numRows_ != numCols_) {
    Error e("Can not diagonalize non-square matrix");
  }
  double* eigenvalues = nullptr;
  eigenvalues = new double[numRows_];

  PMatrix eigenvectors(numRows_, numCols_, numBlocksRows_, numBlocksCols_);

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
  int ia = 0;       // row index from which we diagonalize
  int ja = 0;       // row index from which we diagonalize
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

double PMatrix::squaredNorm() {
    return dot(*this);
}

double PMatrix::norm() {
    return sqrt(squaredNorm());
}

double PMatrix::dot(const PMatrix& that) {
    double scalar = 0.;
    double scalarOut = 0.;
    for ( int i=0; i<numLocalElements_; i++ ) {
        scalar += ( *(mat+i) ) * ( *(that.mat+i) );
    }
    mpi->reduceSum(&scalar, &scalarOut);
    return scalarOut;
}

PMatrix PMatrix::operator-() const {
    PMatrix result = *this;
    for ( int i=0; i<numLocalElements_; i++ ) {
        *(result.mat+i) = - *(result.mat+i);
    }
    return result;
}

PMatrix PMatrix::operator+=(const PMatrix& that) {
    PMatrix result = *this;
    for ( int i=0; i<numLocalElements_; i++ ) {
        *(result.mat+i) += *(that.mat+i);
    }
    return result;
}

PMatrix PMatrix::operator-=(const PMatrix& that) {
    PMatrix result = *this;
    for ( int i=0; i<numLocalElements_; i++ ) {
        *(result.mat+i) -= *(that.mat+i);
    }
    return result;
}
