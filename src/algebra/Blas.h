#ifndef BLAS_H
#define BLAS_H
/** \file   blas.h
 *  \brief  Header file for forward declaration of BLAS and LAPACK functions.
 */

extern "C" {

/* Level 1 BLAS */
double dscal_(const int& n, const double&, double* x, const int&);
double zscal_(const int& n, const std::complex<double>&, std::complex<double>* x, const int&);

void zcopy_(const int& n, const std::complex<double>* x, const int& incx, std::complex<double>* y, const int& incy);
void dcopy_(const int& n, const double*, const int&, double*, const int&);

/* Level 2 BLAS */

void dgemv_(const char& trans, const int& nr, const int& nc, const double& alpha,
             const double* amat, const int& lda, const double* bv,
             const int& incx, const double& beta, double* cv, const int& incy);

void zgemv_(const char& trans, const int& nr, const int& nc,
             const std::complex<double>& alpha, const std::complex<double>* amat, 
             const int& lda, const std::complex<double>* bv,
             const int& incx, const std::complex<double>& beta,
             std::complex<double>* cv, const int& incy);

void daxpy_(const int& n, const double& da, const double* dx, const int& incx, double* dy, const int& incy);

/* Level 3 BLAS */
void dgemm_(const char&, const char&, const int&, const int&, const int&, const double&, const double*,
             const int&, const double*, const int&, const double&, double*, const int&);
void zgemm_(const char&, const char&,
             const int&, const int&, const int&, const std::complex<double>&,
             const std::complex<double>*, const int&, const std::complex<double>*,
             const int&, const std::complex<double>&, std::complex<double>*, const int&);

//Lapack forward declarations
void dgeev_(char* JOBVL, char* JOBVR, int* N, complex* A, int* LDA,
        complex* W, complex* VL, int* LDVL, complex* VR, int* LDVR,
        complex* WORK, int* LWORK, double* RWORK, int* INFO);
void zgeev_(char* JOBVL, char* JOBVR, int* N, complex* A, int* LDA,
        complex* W, complex* VL, int* LDVL, complex* VR, int* LDVR,
        complex* WORK, int* LWORK, double* RWORK, int* INFO);
void dgetrf_(int* M, int* N, complex* A, int* LDA, int* IPIV, int* INFO);
void zgetrf_(int* M, int* N, complex* A, int* LDA, int* IPIV, int* INFO);

}
#endif /* BLAS_H */


