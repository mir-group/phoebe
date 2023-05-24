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
  char uplo = 'U';  // upper triangular
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

  return std::make_tuple(eigenvalues_, eigenvectors);
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
  char uplo = 'U';  // upper triangular
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

  return std::make_tuple(eigenvalues_, eigenvectors);
}

// function to only compute some eigenvectors/values
template <>
std::tuple<std::vector<double>, ParallelMatrix<double>>
                ParallelMatrix<double>::diagonalize(int numEigenvalues_) {

  int numEigenvalues = numEigenvalues_;

  if (numRows_ != numCols_) {
    Error("Cannot diagonalize non-square matrix");
  }
  if ( numBlasRows_ != numBlasCols_ ) {
    Error("Cannot diagonalize via scalapack with a non-square process grid!");
  }
  if ( numRows_ < numEigenvalues) {
    numEigenvalues = numRows_;
  }

  // eigenvalue return container
  double *eigenvalues;
  allocate(eigenvalues, numEigenvalues);

  // NOTE even though we only need numEigenvalues + numEigenvalues worth of z matrix,
  // scalapack absolutely requires this matrix be the same size as the original.
  // It's a huge waste of memory, and we should check to see if another code
  // can get around this.

  // Make a new PMatrix to receive the output
  ParallelMatrix<double> eigenvectors(numRows_, numCols_,
                                      numBlocksRows_,numBlocksCols_);

  // docs for this scalapack function
  // https://www.ibm.com/docs/en/pessl/5.3.0?topic=easva-pdsyevx-pzheevx-selected-
  //  eigenvalues-optionally-eigenvectors-real-symmetric-complex-hermitian-matrix

  char jobz = 'V';  // also eigenvectors
  char uplo = 'U';  // upper triangular
  int ia = 1;       // row index of start of A
  int ja = 1;       // col index of start of A
  int iz = 1;       // row index of start of Z
  int jz = 1;       // row index of start of Z

  int info = 0;     // error code on return
  int m = 0;        // filled on return with number of eigenvalues found
  int nz = 0;       // filled on return with number of eigenvectors found

  char range = 'I'; // compute a range (from smallest to largest) of the eigenvalues
  int il = 1; //numRows_ - numEigenvalues;  // lower eigenvalue index (indexed from 1)
  //if(il == 0) il = 1;
  int iu = numEigenvalues;       // higher eigenvalue index (indexed from 1)
  double vl = -1000;             // not used unless range = V
  double vu = 0;                 // not used unless range = V

  // VARIABLES related to eigenvalue orthogonalization ===============================

  // use abstol value that gives most orthogonal eigenvectors
  char cmach = 'S';           // used by abstol
  double abstol = 2 * pdlamch_(&blacsContext_, &cmach);
  //This array contains the gap between eigenvalues whose
  //           eigenvectors could not be reorthogonalized.
  double *gap;
  allocate(gap, numBlasRows_*numBlasCols_);
  // specifies which eigenvectors should be reorthogonalized.
  double orfac = 0; // -1; //a default value of 10-3 is used.
  // If jobz = 'V', vector ifail is a vector telling us
  // which of the eigenvectors were not converged
  int *ifail;
  allocate(ifail, numRows_);
  // If jobz = 'V', vector iclustr contains the indices of the eigenvectors
  // corresponding to a cluster of eigenvalues that could not be reorthogonalized
  // due to insufficient workspace.
  int *iclustr;
  allocate(iclustr, 2*numBlasRows_*numBlasCols_);

  // we will let pdseyv determine lwork for us. if we run it with
  // lwork = -1 and work of length 1, it will fill work with
  // an appropriate lwork number
  double *work;
  double *iwork;
  int lwork = -1;
  int liwork = -1; // liwork seems not to be determined this way
  allocate(work, 1);
  allocate(iwork, 1);

  pdsyevx_(&jobz, &range, &uplo,  &numRows_, mat, &ia, &ja, &descMat_[0],
        &vl, &vu, &il, &iu, &abstol, &m, &nz, eigenvalues, &orfac,
        eigenvectors.mat, &iz, &jz, &descMat_[0],
        work, &lwork, iwork, &liwork, ifail, iclustr, gap, &info);

  lwork=int(work[0])*100;
  if(lwork < 0) lwork = 2147483000;
  delete[] work;
  delete[] iwork;
  allocate(work, lwork);

  // for some reason, we must calculate liwork manually. The automatic
  // determination doesn't work. Scalapack docs say:
  // max(isizestein, isizestebz)+2n
  //  isizestein = (n+jz-1) +2n+nprocs+1
  //  isizestebz = max(4n, 14, nprocs)
  int isizestein = (numRows_ + jz - 1) + 2 * numRows_ * mpi->getSize() + 1;
  int isizestebz =  std::max(std::max(4*numRows_, 14), mpi->getSize());
  liwork = std::max(isizestein, isizestebz)+2*numRows_;
  allocate(iwork, liwork);

  // first, call and report any negative eigenvalues -------------------
  // we want to note these for convergence reasons
  //  if(context.getNegativeRelaxons()) {

  if(mpi->mpiHead())
    std::cout << "Checking scattering matrix for negative eigenvalues." << std::endl;
  eigenvectors = *this; // use the space of this matrix first to
                // placehold the actual matrix, as it's changed by pdsyevx
  range = 'V';
  jobz = 'N';
  vl = -1.e3;
  vu = 0;
  if(mpi->mpiHead()) mpi->time();

  pdsyevx_(&jobz, &range, &uplo, &numRows_, eigenvectors.mat, &ia, &ja, &descMat_[0],
        &vl, &vu, &il, &iu, &abstol, &m, &nz, eigenvalues, &orfac,
        eigenvectors.mat, &iz, &jz, &descMat_[0],
        work, &lwork, iwork, &liwork, ifail, iclustr, gap, &info);

  if( m > 3 ) { // more than just the zero eigenmode was found
    if(mpi->mpiHead()) {
      Warning("Relaxons diagonalization found " + std::to_string(m) +
                " in the range -1000 < eigenvalues <= 0."
                "\nA proper relaxons solve will not have any negative eigenvalues,"
                "\nand finding them usually indicates the calculation is unconverged."
                "\nYou should run with more wavevectors.");
      std::cout << "These eigenvalues are:" << std::endl;
      for (int i = 0; i < m; i++) {
        std::cout << i << " " << eigenvalues[i] << std::endl;
      }
      std::cout << std::endl;
    }
  }
  if(mpi->mpiHead()) mpi->time();

  // now we perform the regular call to get the largest ones -------------
  if(mpi->mpiHead()) {
    std::cout << "Now computing first " << numEigenvalues <<
        " eigenvalues and vectors of the scattering matrix." << std::endl;
  }

  range = 'I';
  jobz = 'V';
  eigenvectors.zeros();

  if(mpi->mpiHead()) mpi->time();
  // We could make sure these two matrices have identical blacsContexts.
  // However, as dim(Z) must = dim(A), we can just pass A's desc twice.
  pdsyevx_(&jobz, &range, &uplo, &numRows_, mat, &ia, &ja, &descMat_[0],
        &vl, &vu, &il, &iu, &abstol, &m, &nz, eigenvalues, &orfac,
        eigenvectors.mat, &iz, &jz, &descMat_[0],
        work, &lwork, iwork, &liwork, ifail, iclustr, gap, &info);
  if(mpi->mpiHead()) mpi->time();

  if (info != 0 && mpi->mpiHead()) {
    if(info < 0) {
      std::cout << "Developer Error: "
                "One of the input params to PDSYEVX is wrong!" << std::endl;
    } else {
      if ( (info/2 % 2) != 0 ) {
        std::cout << "PDSYEVX had insufficient memory to "
                        "reorthogonalize the eigenvectors." << std::endl;
      } else if ( (info/4 % 2) != 0 ) {
        std::cout << "PDSYEVX had insufficient memory to " <<
                        "compute all of the eigenvectors." << std::endl;
      } else if ( (info/8 % 2) != 0 ) {
        std::cout << "PDSYEVX had an internal error and " <<
                "failed to compute the eigenvalues." << std::endl;
      }
    }
    Error("PDSYEVX failed.", info);
  }

  // copy to return containers and free the no longer used containers.
  std::vector<double> eigenvalues_(numEigenvalues);
  for (int i = 0; i < numEigenvalues; i++) {
    eigenvalues_[i] = *(eigenvalues + i);
//    std::cout << " eigenvalue " << eigenvalues_[i] << std::endl;
//    std::cout << " eigenvector " << numRows_ << " ";
//    for(int is = 0; is < numRows_; is++) {
//      std::cout << " " << eigenvectors(is,i);
//    }
//    std::cout << std::endl;
  }
  delete[] eigenvalues; delete[] ifail;
  delete[] work; //delete[] iwork;
  delete[] iclustr; delete[] gap;

  // note that the scattering matrix now has different values
  // it's going to be the upper triangle and diagonal of A
  return std::make_tuple(eigenvalues_, eigenvectors);
}

#endif  // MPI_AVAIL
