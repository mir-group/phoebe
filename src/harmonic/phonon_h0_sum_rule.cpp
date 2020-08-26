#include "phonon_h0.h"

#include <math.h>
#include <complex>
#include <iostream>

#include "constants.h"
#include "eigen.h"
#include "exceptions.h"

void PhononH0::setAcousticSumRule(const std::string &sumRule) {
  double norm2;
  //  VectorXi u_less(6*3*numAtoms)
  //  indices of the vectors u that are not independent to the preceding ones
  //  n_less = number of such vectors
  //
  //  Tensor<int> ind_v(:,:,:)
  //  Tensor<double> v(:,:)
  //  These are the "vectors" associated with symmetry conditions, coded by
  //  indicating the positions (i.e. the seven indices) of the non-zero
  //  elements (there should be only 2 of them) and the value of that element
  //  We do so in order to limit the amount of memory used.
  //
  //  Tensor<double> zeu_u(6*3,3,3,numAtoms)
  //  These are vectors associated with the sum rules on effective charges
  //
  //  Tensor<int> zeu_less(6*3)
  //  indices of zeu_u vectors that are not independent to the preceding ones
  //  ! nzeu_less = number of such vectors

  std::string sr = sumRule;
  std::transform(sr.begin(), sr.end(), sr.begin(), ::tolower);

  if (sr == "") {
    return;
  }

  if ((sr != "simple") && (sr != "crystal")) {
    Error e("invalid Acoustic Sum Rule", 1);
  }

  if (mpi->mpiHead()) {
    std::cout << "Start imposing " << sumRule << " acoustic sum rule."
	      << std::endl;
  }

  if (sr == "simple") {
    // Simple Acoustic Sum Rule on effective charges

    double sum;

    for (long i = 0; i < 3; i++) {
      for (long j = 0; j < 3; j++) {
        sum = 0.;
        for (long na = 0; na < numAtoms; na++) {
          sum += bornCharges(na, i, j);
        }
        for (long na = 0; na < numAtoms; na++) {
          bornCharges(na, i, j) -= sum / numAtoms;
        }
      }
    }

    // Simple Acoustic Sum Rule on force constants in real space

    for (long i = 0; i < 3; i++) {
      for (long j = 0; j < 3; j++) {
        for (long na = 0; na < numAtoms; na++) {
          sum = 0.;
          for (long nb = 0; nb < numAtoms; nb++) {
            for (long n1 = 0; n1 < qCoarseGrid(0); n1++) {
              for (long n2 = 0; n2 < qCoarseGrid(1); n2++) {
                for (long n3 = 0; n3 < qCoarseGrid(2); n3++) {
                  sum += forceConstants(i, j, n1, n2, n3, na, nb);
                }
              }
            }
          }
          forceConstants(i, j, 0, 0, 0, na, na) -= sum;
        }
      }
    }
  } else {
    // Acoustic Sum Rule on effective charges

    // generating the vectors of the orthogonal of the subspace to project
    // the effective charges matrix on

    Eigen::Tensor<double, 4> zeu_u(6 * 3, 3, 3, numAtoms);
    zeu_u.setZero();
    Eigen::Tensor<double, 3> zeu_new(3, 3, numAtoms);
    zeu_new.setZero();

    for (long i = 0; i < 3; i++) {
      for (long j = 0; j < 3; j++) {
        for (long iat = 0; iat < numAtoms; iat++) {
          zeu_new(i, j, iat) = bornCharges(iat, i, j);
        }
      }
    }

    long p = 0;
    for (long i = 0; i < 3; i++) {
      for (long j = 0; j < 3; j++) {
        for (long iat = 0; iat < numAtoms; iat++) {
          // These are the 3*3 vectors associated with the
          // translational acoustic sum rules
          zeu_u(p, i, j, iat) = 1.;
        }
        p += 1;
      }
    }

    // Gram-Schmidt orthonormalization of the set of vectors created.

    // temporary vectors
    Eigen::Tensor<double, 3> zeu_w(3, 3, numAtoms), zeu_x(3, 3, numAtoms);
    Eigen::Tensor<double, 3> tempZeu(3, 3, numAtoms);
    // note: it's important to initialize these tensors
    zeu_w.setZero();
    zeu_x.setZero();
    tempZeu.setZero();
    Eigen::VectorXi zeu_less(6 * 3);
    zeu_less.setZero();
    double scal;
    long nzeu_less = 0;
    long r;

    for (long k = 0; k < p; k++) {
      for (long i = 0; i < 3; i++) {
        for (long j = 0; j < 3; j++) {
          for (long iat = 0; iat < numAtoms; iat++) {
            zeu_w(i, j, iat) = zeu_u(k, i, j, iat);
            zeu_x(i, j, iat) = zeu_u(k, i, j, iat);
          }
        }
      }

      for (long q = 0; q < k - 1; q++) {
        r = 1;
        for (long izeu_less = 0; izeu_less < nzeu_less; izeu_less++) {
          if (zeu_less(izeu_less) == q) {
            r = 0;
          };
        }
        if (r != 0) {
          for (long i = 0; i < 3; i++) {
            for (long j = 0; j < 3; j++) {
              for (long iat = 0; iat < numAtoms; iat++) {
                tempZeu(i, j, iat) = zeu_u(q, i, j, iat);
              }
            }
          }
          // i.e. zeu_u(q,:,:,:)
          sp_zeu(zeu_x, tempZeu, scal);
          zeu_w -= scal * tempZeu;
        }
      }
      sp_zeu(zeu_w, zeu_w, norm2);

      if (norm2 > 1.0e-16) {
        for (long i = 0; i < 3; i++) {
          for (long j = 0; j < 3; j++) {
            for (long iat = 0; iat < numAtoms; iat++) {
              zeu_u(k, i, j, iat) = zeu_w(i, j, iat) / sqrt(norm2);
            }
          }
        }
      } else {
        zeu_less(nzeu_less) = k;
        nzeu_less += 1;
      }
    }

    // Projection of the effective charge "vector" on the orthogonal of the
    // subspace of the vectors verifying the sum rules

    zeu_w.setZero();
    for (long k = 0; k < p; k++) {
      r = 1;
      for (long izeu_less = 0; izeu_less < nzeu_less; izeu_less++) {
        if (zeu_less(izeu_less) == k) {
          r = 0;
        };
      }
      if (r != 0) {
        // copy vector
        for (long i = 0; i < 3; i++) {
          for (long j = 0; j < 3; j++) {
            for (long iat = 0; iat < numAtoms; iat++) {
              zeu_x(i, j, iat) = zeu_u(k, i, j, iat);
            }
          }
        }
        // get rescaling factor
        sp_zeu(zeu_x, zeu_new, scal);
        // rescale vector
        for (long i = 0; i < 3; i++) {
          for (long j = 0; j < 3; j++) {
            for (long iat = 0; iat < numAtoms; iat++) {
              zeu_w(i, j, iat) += scal * zeu_u(k, i, j, iat);
            }
          }
        }
      }
    }

    // Final substraction of the former projection to the initial zeu, to
    // get the new "projected" zeu

    zeu_new -= zeu_w;
    sp_zeu(zeu_w, zeu_w, norm2);
    if ( mpi->mpiHead() ) {
      std::cout << "Norm of the difference between old and new effective "
                   "charges: " << sqrt(norm2) << std::endl;
    }

    for (long i = 0; i < 3; i++) {
      for (long j = 0; j < 3; j++) {
        for (long iat = 0; iat < numAtoms; iat++) {
          bornCharges(iat, i, j) = zeu_new(i, j, iat);
        }
      }
    }

    ///////////////////////////////////////////////////////////////////////////

    // Acoustic Sum Rule on force constants

    // generating the vectors of the orthogonal of the subspace to project
    // the force-constants matrix on

    int nr1 = qCoarseGrid(0);
    int nr2 = qCoarseGrid(1);
    int nr3 = qCoarseGrid(2);

    Eigen::Tensor<double,8> uvec(18*numAtoms,nr1,nr2,nr3,3,3,numAtoms,numAtoms);
    uvec.setZero();

    Eigen::Tensor<double,7> frc_new(nr1,nr2,nr3,3,3,numAtoms,numAtoms);
    for (long nb = 0; nb < numAtoms; nb++) {
      for (long na = 0; na < numAtoms; na++) {
        for (long j = 0; j < 3; j++) {
          for (long i = 0; i < 3; i++) {
            for (long n3 = 0; n3 < nr3; n3++) {
              for (long n2 = 0; n2 < nr2; n2++) {
                for (long n1 = 0; n1 < nr1; n1++) {
                  frc_new(n1,n2,n3,i,j,na,nb) = forceConstants(i,j,n1,n2,n3,na,nb);
                }
              }
            }
          }
        }
      }
    }

    p = 0;
    for ( int na=0; na<numAtoms; na++ ) {
      for ( int j=0; j<3; j++ ) {
        for ( int i=0; i<3; i++ ) {
          // These are the 3*3*nat vectors associated with the
          // translational acoustic sum rules
          for ( int nb=0; nb<numAtoms; nb++ ) {
            for ( int n3=0; n3<nr3; n3++ ) {
              for ( int n2=0; n2<nr2; n2++ ) {
                for ( int n1=0; n1<nr1; n1++ ) {
                  uvec(p,n1,n2,n3,i,j,na,nb) = 1.;
                }
              }
            }
          }
          p += 1;
        }
      }
    }

    Eigen::Tensor<int,3> ind_v(9*numAtoms*numAtoms*nr1*nr2*nr3,2,7);
    Eigen::Tensor<double,2> v(9*numAtoms*numAtoms*nr1*nr2*nr3,2);
    v.setZero();
    ind_v.setZero();

    int m = 0;
    int q, l;

    for ( int nb=1; nb<=numAtoms; nb++ ) {
      for ( int na=1; na<=numAtoms; na++ ) {
        for ( int j=1; j<=3; j++ ) {
          for ( int i=1; i<=3; i++ ) {
            for ( int n3=1; n3<=nr3; n3++ ) {
              for ( int n2=1; n2<=nr2; n2++ ) {
                for ( int n1=1; n1<=nr1; n1++ ) {
                  // These are the vectors associated with the symmetry constraints
                  q = 1;
                  l = 1;
                  while ( ( l <= m ) && ( q != 0 ) ) {
                    if ( (ind_v(l-1,1-1,1-1)==n1) && (ind_v(l-1,1-1,2-1)==n2) &&
                        (ind_v(l-1,1-1,3-1)==n3) && (ind_v(l-1,1-1,4-1)==i) &&
                        (ind_v(l-1,1-1,5-1)==j)  && (ind_v(l-1,1-1,6-1)==na) &&
                        (ind_v(l-1,1-1,7-1)==nb)) { q=0; }
                    if ( (ind_v(l-1,2-1,1-1)==n1) && (ind_v(l-1,2-1,2-1)==n2) &&
                        (ind_v(l-1,2-1,3-1)==n3) && (ind_v(l-1,2-1,4-1)==i) &&
                        (ind_v(l-1,2-1,5-1)==j)  && (ind_v(l-1,2-1,6-1)==na) &&
                        (ind_v(l-1,2-1,7-1)==nb)) {q = 0;}
                    l += 1;
                  }
                  if ( ( n1==((nr1+1-n1)%nr1)+1 ) &&
                      ( n2==((nr2+1-n2)%nr2)+1 ) &&
                      ( n3==((nr3+1-n3)%nr3)+1 ) &&
                      ( i==j ) && ( na==nb ) ) { q=0; }
                  if ( q != 0 ) {
                    m += 1;
                    ind_v(m-1,1-1,1-1) = n1;
                    ind_v(m-1,1-1,2-1) = n2;
                    ind_v(m-1,1-1,3-1) = n3;
                    ind_v(m-1,1-1,4-1) = i;
                    ind_v(m-1,1-1,5-1) = j;
                    ind_v(m-1,1-1,6-1) = na;
                    ind_v(m-1,1-1,7-1) = nb;
                    v(m-1,1-1) = 1. / sqrt(2.);
                    ind_v(m-1,2-1,1-1) = ((nr1+1-n1) % nr1) + 1;
                    ind_v(m-1,2-1,2-1) = ((nr2+1-n2) % nr2) + 1;
                    ind_v(m-1,2-1,3-1) = ((nr3+1-n3) % nr3) + 1;
                    ind_v(m-1,2-1,4-1) = j;
                    ind_v(m-1,2-1,5-1) = i;
                    ind_v(m-1,2-1,6-1) = nb;
                    ind_v(m-1,2-1,7-1) = na;
                    v(m-1,2-1) = - 1. / sqrt(2.);
                  }
                }
              }
            }
          }
        }
      }
    }

    // Gram-Schmidt orthonormalization of the set of vectors created.
    // Note that the vectors corresponding to symmetry constraints are already
    // orthonormalized by construction.

    int n_less = 0;
    Eigen::Tensor<double,7> w(nr1,nr2,nr3,3,3,numAtoms,numAtoms);
    Eigen::Tensor<double,7> x(nr1,nr2,nr3,3,3,numAtoms,numAtoms);
    w.setZero();
    x.setZero();

    Eigen::VectorXi u_less(6*3*numAtoms);
    u_less.setZero();

    int n1, n2, n3, i, j, na, nb, na1, i1, j1;

    for ( int k=1; k<=p; k++ ) {
      //        w(:,:,:,:,:,:,:) = uvec(k-1,:,:,:,:,:,:,:);
      //          x(:,:,:,:,:,:,:) = uvec(k-1,:,:,:,:,:,:,:);
      for ( int nb=0; nb<numAtoms; nb++ ) {
        for ( int na=0; na<numAtoms; na++ ) {
          for ( int j=0; j<3; j++ ) {
            for ( int i=0; i<3; i++ ) {
              for ( int n3=0; n3<nr3; n3++ ) {
                for ( int n2=0; n2<nr2; n2++ ) {
                  for ( int n1=0; n1<nr1; n1++ ) {
                    w(n1,n2,n3,i,j,na,nb) = uvec(k-1,n1,n2,n3,i,j,na,nb);
                    x(n1,n2,n3,i,j,na,nb) = uvec(k-1,n1,n2,n3,i,j,na,nb);
                  }
                }
              }
            }
          }
        }
      }

      for ( int l=1; l<=m; l++ ) {
//        call sp2(x,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal);
        scal = 0.;
        for ( int r=1; r<=2; r++ ) {
          n1 = ind_v(l-1,r-1,1-1);
          n2 = ind_v(l-1,r-1,2-1);
          n3 = ind_v(l-1,r-1,3-1);
          i = ind_v(l-1,r-1,4-1);
          j = ind_v(l-1,r-1,5-1);
          na = ind_v(l-1,r-1,6-1);
          nb = ind_v(l-1,r-1,7-1);
          scal += w(n1-1,n2-1,n3-1,i-1,j-1,na-1,nb-1) * v(l-1,r-1);
        }

        for ( int r=1; r<=2; r++ ) {
          n1 = ind_v(l-1,r-1,1-1);
          n2 = ind_v(l-1,r-1,2-1);
          n3 = ind_v(l-1,r-1,3-1);
          i = ind_v(l-1,r-1,4-1);
          j = ind_v(l-1,r-1,5-1);
          na = ind_v(l-1,r-1,6-1);
          nb = ind_v(l-1,r-1,7-1);
          w(n1-1,n2-1,n3-1,i-1,j-1,na-1,nb-1) -= scal * v(l-1,r-1);
        }
      }

      if ( k <= ( 9*numAtoms) ) {
        na1 = k % numAtoms;
        if ( na1==0 ) { na1 = numAtoms; };
        j1 = (((k-na1)/numAtoms) % 3) + 1;
        i1 = (((((k-na1)/numAtoms)-j1+1)/3) % 3) + 1;
      } else {
        q = k - 9 * numAtoms;
        //            if ( n == 4 ) {
          //               na1 = q % nat;
        //               if ( na1 == 0 ) { na1 = nat; }
        //               i1 = MOD((q-na1)/nat,3)+1
        //            } else {
        na1 = q % numAtoms;
        if ( na1 == 0 ) na1 = numAtoms;
        j1 = (((q-na1)/numAtoms) % 3) + 1;
        i1 = (((((q-na1)/numAtoms)-j1+1)/3) % 3) + 1;
        //            }
      }
      for ( int q=1; q<=k-1; q++ ) {
        r = 1;
        for ( int i_less=1; i_less<=n_less; i_less++ ) {
          if ( u_less(i_less-1) == q ) { r = 0; }
        }
        if ( r != 0 ) {
//          call sp3(x,uvec(q-1,:,:,:,:,:,:,:),i1,na1,nr1,nr2,nr3,nat,scal)
          scal = 0.;
          for ( int nb=0; nb<numAtoms; nb++ ) {
            for ( int j=0; j<3; j++ ) {
              for ( int n3=0; n3<nr3; n3++ ) {
                for ( int n2=0; n2<nr2; n2++ ) {
                  for ( int n1=0; n1<nr1; n1++ ) {
                    scal += x(n1,n2,n3,i1-1,j,na1-1,nb)
                        * uvec(q-1,n1,n2,n3,i1-1,j,na1-1,nb);
                  }
                }
              }
            }
          }

          //               w(:,:,:,:,:,:,:) -= scal * uvec(q-1,:,:,:,:,:,:,:)
          for ( int nb=0; nb<numAtoms; nb++ ) {
            for ( int na=0; na<numAtoms; na++ ) {
              for ( int j=0; j<3; j++ ) {
                for ( int i=0; i<3; i++ ) {
                  for ( int n3=0; n3<nr3; n3++ ) {
                    for ( int n2=0; n2<nr2; n2++ ) {
                      for ( int n1=0; n1<nr1; n1++ ) {
                        w(n1,n2,n3,i,j,na,nb) -= scal * uvec(q-1,n1,n2,n3,i,j,na,nb);
                      }
                    }
                  }
                }
              }
            }
          }

        }
      }

//      call sp1(w,w,nr1,nr2,nr3,nat,norm2)
      norm2 = 0.;
      for ( int nb=0; nb<numAtoms; nb++ ) {
        for ( int na=0; na<numAtoms; na++ ) {
          for ( int j=0; j<3; j++ ) {
            for ( int i=0; i<3; i++ ) {
              for ( int n3=0; n3<nr3; n3++ ) {
                for ( int n2=0; n2<nr2; n2++ ) {
                  for ( int n1=0; n1<nr1; n1++ ) {
                    norm2 += w(n1,n2,n3,i,j,na,nb)
                        * w(n1,n2,n3,i,j,na,nb);
                  }
                }
              }
            }
          }
        }
      }

      if ( norm2 > 1.0e-16 ) {
        // uvec(k-1,:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) / sqrt(norm2);
        for ( int nb=0; nb<numAtoms; nb++ ) {
          for ( int na=0; na<numAtoms; na++ ) {
            for ( int j=0; j<3; j++ ) {
              for ( int i=0; i<3; i++ ) {
                for ( int n3=0; n3<nr3; n3++ ) {
                  for ( int n2=0; n2<nr2; n2++ ) {
                    for ( int n1=0; n1<nr1; n1++ ) {
                      uvec(k-1,n1,n2,n3,i,j,na,nb) =
                          w(n1,n2,n3,i,j,na,nb) / sqrt(norm2);
                    }
                  }
                }
              }
            }
          }
        }
      } else {
        n_less += 1;
        u_less(n_less-1) = k;
      }
    }

    // Projection of the force-constants "vector" on the orthogonal of the
    // subspace of the vectors verifying the sum rules and symmetry contraints

    w.setZero();
    for ( int l=1; l<=m; l++ ) {

//      call sp2(frc_new,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
      scal = 0.;
      for ( int i=1; i<=2; i++ ) {
        scal += frc_new(ind_v(l-1,i-1,1-1), ind_v(l-1,i-1,2-1), ind_v(l-1,i-1,3-1),
            ind_v(l-1,i-1,4-1), ind_v(l-1,i-1,5-1), ind_v(l-1,i-1,6-1),
            ind_v(l-1,i-1,7-1)) * v(l-1,i-1);
      }

      for ( int r=1; r<=2; r++ ) {
        n1 = ind_v(l-1,r-1,1-1);
        n2 = ind_v(l-1,r-1,2-1);
        n3 = ind_v(l-1,r-1,3-1);
        i = ind_v(l-1,r-1,4-1);
        j = ind_v(l-1,r-1,5-1);
        na = ind_v(l-1,r-1,6-1);
        nb = ind_v(l-1,r-1,7-1);
        w(n1-1,n2-1,n3-1,i-1,j-1,na-1,nb-1) += scal * v(l-1,r-1);
      }
    }
    for ( int k=1; k<=p; k++ ) {
      r = 1;
      for ( int i_less=1; i_less<=n_less; i_less++ ) {
        if ( u_less(i_less-1) == k ) { r = 0; }
      }
      if ( r != 0 ) {

        scal = 0.;
        for ( int nb=0; nb<numAtoms; nb++ ) {
          for ( int na=0; na<numAtoms; na++ ) {
            for ( int j=0; j<3; j++ ) {
              for ( int i=0; i<3; i++ ) {
                for ( int n3=0; n3<nr3; n3++ ) {
                  for ( int n2=0; n2<nr2; n2++ ) {
                    for ( int n1=0; n1<nr1; n1++ ) {
                      scal += uvec(k-1,n1,n2,n3,i,j,na,nb)
                          * frc_new(n1,n2,n3,i,j,na,nb);
                    }
                  }
                }
              }
            }
          }
        }

        for ( int nb=0; nb<numAtoms; nb++ ) {
          for ( int na=0; na<numAtoms; na++ ) {
            for ( int j=0; j<3; j++ ) {
              for ( int i=0; i<3; i++ ) {
                for ( int n3=0; n3<nr3; n3++ ) {
                  for ( int n2=0; n2<nr2; n2++ ) {
                    for ( int n1=0; n1<nr1; n1++ ) {
                      w(n1,n2,n3,i,j,na,nb) +=
                          scal * uvec(k-1,n1,n2,n3,i,j,na,nb);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    // Final substraction of the former projection to the initial frc,
    // to get the new "projected" frc

    frc_new -= w;
    scal = 0.;
    for ( int nb=0; nb<numAtoms; nb++ ) {
      for ( int na=0; na<numAtoms; na++ ) {
        for ( int j=0; j<3; j++ ) {
          for ( int i=0; i<3; i++ ) {
            for ( int n3=0; n3<nr3; n3++ ) {
              for ( int n2=0; n2<nr2; n2++ ) {
                for ( int n1=0; n1<nr1; n1++ ) {
                  scal += w(n1,n2,n3,i,j,na,nb) * w(n1,n2,n3,i,j,na,nb);
                }
              }
            }
          }
        }
      }
    }

    if (mpi->mpiHead()) {
      std::cout << "Norm of the difference between old and new "
                   "force-constants: " << sqrt(scal) << std::endl;
    }

//    forceConstants = frc_new;
    for (long nb = 0; nb < numAtoms; nb++) {
      for (long na = 0; na < numAtoms; na++) {
        for (long i = 0; i < 3; i++) {
          for (long j = 0; j < 3; j++) {
            for (long n3 = 0; n3 < nr3; n3++) {
              for (long n2 = 0; n2 < nr2; n2++) {
                for (long n1 = 0; n1 < nr1; n1++) {
                  forceConstants(i,j,n1,n2,n3,na,nb) = frc_new(n1,n2,n3,i,j,na,nb);
                }
              }
            }
          }
        }
      }
    }

  }
  if (mpi->mpiHead()) {
    std::cout << "Finished imposing " << sumRule << " acoustic sum rule."
	      << std::endl;
  }
}

void PhononH0::sp_zeu(Eigen::Tensor<double, 3> &zeu_u,
                      Eigen::Tensor<double, 3> &zeu_v, double &scal) {
  // get the scalar product of two effective charges matrices zeu_u and zeu_v
  // (considered as vectors in the R^(3*3*nat) space)

  scal = 0.;
  for (long i = 0; i < 3; i++) {
    for (long j = 0; j < 3; j++) {
      for (long na = 0; na < numAtoms; na++) {
        scal += zeu_u(i, j, na) * zeu_v(i, j, na);
      }
    }
  }
}
