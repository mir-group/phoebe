#ifndef BLACS_H
#define BLACS_H
/** \file   blacs.h
 *  \brief  Header file for forward declaration of BLACS and SCALAPACK
 *  functions.
 */

#include <complex>

extern "C" {

void blacs_get_(int *, int *, int *);
void blacs_pinfo_(int *, int *);
void blacs_gridinit_(int *, char *, int *, int *);
void blacs_gridinfo_(int *, int *, int *, int *, int *);
void descinit_(int *, int *, int *, int *, int *, int *, int *, int *, int *,
               int *);
void blacs_gridexit_(const int *);
int numroc_(int *, int *, int *, int *, int *);

//void pdelset_(double *, int *, int *, int *, double *);
//void pdelget_(char *, char *, double *, double *, const int *, const int *,
//              const int *);
void infog2l_(const int *, const int *, const int *, const int *, const int *,
              const int *, const int *, int *, int *, int *, int *);

int indxg2p_(int *, int *, int *, int *, int *);

int indxl2g_(int *, int *, int *, int *, int *);

void pdgemm_(const char *, const char *, int *, int *, const int *, double *,
             double *, int *, int *, const int *, double *, int *, int *,
             const int *, double *, double *, int *, int *, int *);
void pdsyev_(char *, char *, int *, double *, int *, int *, int *, double *,
             double *, int *, int *, int *, double *, int *, int *);
//void pzelset_(std::complex<double> *, int *, int *, int *,
//              std::complex<double> *);
//void pzelget_(char *, char *, std::complex<double> *, std::complex<double> *,
//              int *, int *, int *);
void pzgemm_(const char *, const char *, int *, int *, const int *,
             std::complex<double> *, std::complex<double> *, int *, int *,
             const int *, std::complex<double> *, int *, int *, const int *,
             std::complex<double> *, std::complex<double> *, int *, int *,
             int *);
void pzheev_(char *, char *, int *, std::complex<double> *, int *, int *, int *,
             double *, std::complex<double> *, int *, int *, int *,
             std::complex<double> *, int *, std::complex<double> *, int *,
             int *);
};

#endif /* BLACS_H */
