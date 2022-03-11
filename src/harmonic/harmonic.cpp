#include "harmonic.h"

std::tuple<DoubleView2D, ComplexView3D, ComplexView4D>
HarmonicHamiltonian::kokkosBatchedDiagonalizeWithVelocities(
    const DoubleView2D &cartesianCoordinates) {

  int numK = cartesianCoordinates.extent(0);

  double delta = 1.0e-8;
  double threshold = 0.000001 / energyRyToEv;// = 1 micro-eV

  // prepare all the wavevectors at which we need the hamiltonian
  DoubleView2D allVectors("wavevectors_k", numK * 7, 3);
  Kokkos::parallel_for(
      "elHamiltonianShifted_d", numK, KOKKOS_LAMBDA(int iK) {
        for (int i = 0; i < 3; ++i) {
          allVectors(iK * 7, i) = cartesianCoordinates(iK, i);
        }
        for (int iDir = 0; iDir < 3; ++iDir) {
          for (int i = 0; i < 3; ++i) {
            allVectors(iK * 7 + iDir * 2 + 1, i) = cartesianCoordinates(iK, i);
            allVectors(iK * 7 + iDir * 2 + 2, i) = cartesianCoordinates(iK, i);
            if (iDir == i) {
              allVectors(iK * 7 + iDir * 2 + 1, i) += delta;
              allVectors(iK * 7 + iDir * 2 + 2, i) -= delta;
            }
          }
        }
      });

  // compute the electronic properties at all wavevectors
  ComplexView3D allHamiltonians = kokkosBatchedBuildBlochHamiltonian(allVectors);

  int numBands = allHamiltonians.extent(1);

  // save energies and eigenvectors to results
  DoubleView2D resultEnergies("energies", numK, numBands);
  ComplexView3D resultEigenvectors("eigenvectors", numK, numBands, numBands);
  ComplexView4D resultVelocities("velocities", numK, numBands, numBands, 3);

  // put the Hamiltonian matrix in resultEigenvectors
  Kokkos::parallel_for(
      "eigenvectors", Range3D({0, 0, 0}, {numK, numBands, numBands}),
      KOKKOS_LAMBDA(int iK, int m, int n) {
        resultEigenvectors(iK, m, n) = allHamiltonians(iK*7, m, n);
      });
  // now, diagonalize the H matrix in place
  kokkosZHEEV(resultEigenvectors, resultEnergies);

  // these are temporary "scratch" memory spaces
  ComplexView3D der("der", numK, numBands, numBands);
  ComplexView3D tmpV("tmpV", numK, numBands, numBands);

  for (int i = 0; i < 3; ++i) {

    // To build the velocity operator, we first compute the matrix
    // der = dH/dk = ( H(k+dk)-H(k-dk) ) / 2dk
    // we also reinforce the Hermiticity (better numerical properties)
    Kokkos::parallel_for(
        "der", Range3D({0, 0, 0}, {numK, numBands, numBands}),
        KOKKOS_LAMBDA(int iK, int m, int n) {
          ComplexView2D HPlus = Kokkos::subview(allHamiltonians, iK * 7 + i * 2 + 1, Kokkos::ALL, Kokkos::ALL);
          ComplexView2D HMins = Kokkos::subview(allHamiltonians, iK * 7 + i * 2 + 2, Kokkos::ALL, Kokkos::ALL);
          der(iK, m, n) = 0.25 / delta * ((HPlus(m, n) - HMins(m, n))
                                          + Kokkos::conj(HPlus(n, m) - HMins(n, m)));
        });

    // Now we complete the Hellman-Feynman theorem
    // and compute the velocity as v = U(k)^* der * U(k)
    Kokkos::parallel_for(
        "tmpV", Range3D({0, 0, 0}, {numK, numBands, numBands}),
        KOKKOS_LAMBDA(int iK, int m, int n) {
          auto L = Kokkos::subview(resultEigenvectors, iK, Kokkos::ALL, Kokkos::ALL);
          auto R = Kokkos::subview(der, iK, Kokkos::ALL, Kokkos::ALL);
          auto A = Kokkos::subview(tmpV, iK, Kokkos::ALL, Kokkos::ALL);
          Kokkos::complex<double> tmp(0.,0.);
          for (int l = 0; l < numBands; ++l) {
            tmp += Kokkos::conj(L(l,m)) * R(l, n);
          }
          A(m, n) = tmp;
        });

    Kokkos::parallel_for(
        "vel", Range3D({0, 0, 0}, {numK, numBands, numBands}),
        KOKKOS_LAMBDA(int iK, int m, int n) {
          double norm = 0.;
          for (int i=0; i<3; ++i) {
            norm += cartesianCoordinates(iK,i)*cartesianCoordinates(iK,i);
          }
          Kokkos::complex<double> tmp(0.,0.);
          if ( norm > 1.0e-6 ) {// skip the gamma point
            for (int l = 0; l < numBands; ++l) {
              tmp += tmpV(iK, m, l) * resultEigenvectors(iK, l, n);
            }
          }
          resultVelocities(iK, m, n, i) = tmp;
        });
  }

  // deallocate the scratch
  Kokkos::resize(der, 0, 0, 0);
  Kokkos::resize(tmpV, 0, 0, 0);

  // the degeneracy diagonalization is a bit tricky.
  // in the way it was implemented for the CPU case can't work for GPU
  // that's because it's hard to do several diagonalizations on GPU
  // (due to the fact that the GPU can't launch kernels by itself)

  // ALGO:
  // - find all the degrees of degeneracy that have to be run.
  // - loop over the degree of degeneracy
  //   - set up all sub-matrices of velocities of rank (ndegree)
  //   - diagonalize them
  //   - rotate velocity elements

  // degCount
  IntView2D degCount("tmpV", numK, numBands);
  // maxDegeneracy is the largest degenerate block of energies
  int maxDegeneracy = 1;
  // here I do two things: find the maxDegeneracy for the current system
  // and set degCount such that, degCount(ik,ib) = 0 if it is a state that is
  // not the starting point of a degenerate block, or degCount(ik,ib) = iDeg
  // if ik,ib is at the start of a block iDegxiDeg that needs to be
  // diagonalized, with iDeg >= 2
  Kokkos::parallel_reduce(
      "Loop1", numK,
      KOKKOS_LAMBDA(const int &iK, int &iDeg) {
        for (int ib = 0; ib < numBands; ib++) {
          // first, we check if the band is degenerate, and the size of the
          // degenerate subspace
          int sizeSubspace = 1;
          for (int ib2 = ib + 1; ib2 < numBands; ib2++) {
            // I consider bands degenerate if their frequencies are the same
            // within 0.0001 cm^-1
            if (abs(resultEnergies(iK, ib) - resultEnergies(iK, ib2)) > threshold) {
              break;
            }
            ++sizeSubspace;
          }

          if (sizeSubspace > iDeg) {
            iDeg = sizeSubspace; // iDeg = std::max(iDeg, sizeSubspace);
          }
          if (sizeSubspace == 1) {
            degCount(iK, ib) = 0;
          } else {
            degCount(iK, ib) = sizeSubspace;
            for (int i=1; i<sizeSubspace; ++i) {
              degCount(iK, ib+i) = 0;
            }
          }

          // we skip the bands in the subspace, since we corrected them already
          ib += sizeSubspace - 1;
        }
      },
      Kokkos::Max<int>(maxDegeneracy));

  // we now do a diagonalization of all degenerate velocity sub-blocks
  // since cuda can only launch the diagonalization of N matrices at fixed
  // matrix size, we now do a loop over the degree of degeneracy
  for (int iDeg = 2; iDeg<=maxDegeneracy; ++iDeg) {

    // count how many degenerate points we have at this degree of degeneracy
    int numMatrices = 0;
    Kokkos::parallel_reduce("Loop1", numK,
        KOKKOS_LAMBDA(const int &iK, int &iMat) {
          for (int ib = 0; ib < numBands; ib++) {
            if (degCount(iK, ib) == iDeg) {
              ++iMat;
            }
          }
        }, numMatrices);

    /**
     * The next part is a bit convoluted because it's tricky to do in parallel
     * We want to build N matrices of size iDegxiDeg containing the velocity
     * degenerate blocks.
     * To this aim, we need two functions funcK, funcB. Such that funcK(iMat),
     * funcB(iMat) return the iK and ib index locating the start of the block
     * that will be diagonalized by the iMat-th matrix
     *
     * We take advantage of an exclusive scan. Say we have a vector
     * 0 0 0 1 0 0 0 0 1;
     * which spans the bloch states, and 1 locates the beginning of a block
     * of size iDeg that we want to diagonalize. The indices locating the "1",
     * i.e. [3, 8], are given by an exclusive scan:
     * 0 0 0 0 1 1 1 1 1 2
     * So, scan(ik,ib) = iMat. And from this, I can do the inverse mapping.
     */

    IntView1D scan("scan", numK*numBands);
    Kokkos::parallel_for(
        Range2D({0,0},{numK, numBands}), KOKKOS_LAMBDA(int ik, int ib) {
          if ( iDeg == degCount(ik,ib) ) {
            scan(ik*numBands+ib) = 1;
          } else {
            scan(ik*numBands+ib) = 0;
          }
        });

    Kokkos::parallel_scan("scan", numK*numBands, KOKKOS_LAMBDA(const int i,
                                                 int& update, const bool final) {
          // Load old value in case we update it before accumulating
          const int val_i = scan(i);
          if (final) {
            scan(i) = update; // only update array on final pass
          }
          // For exclusive scan, change the update value after
          // updating array, like we do here. For inclusive scan,
          // change the update value before updating array.
          update += val_i;
        });

    IntView1D funcK("funcK", numMatrices);
    IntView1D funcB("funcB", numMatrices);
    Kokkos::parallel_for(
        "funcKB", Range2D({0,0},{numK, numBands}), KOKKOS_LAMBDA(int ik, int ib) {
          int ikb = ik*numBands + ib;
          int iMat = scan(ikb);
          if (degCount(ik, ib) == iDeg) { // if we actually had a 1
            // the if is needed to avoid race conditions in funcK/funcB
            funcK(iMat) = ik;
            funcB(iMat) = ib;
          }
        });
    Kokkos::resize(scan, 0);

    // now, I copy the degenerate blocks of the velocity operator in a scratch
    // memory space "subVelocity".
    ComplexView3D subVelocity("subVel", 3*numMatrices, iDeg, iDeg);
    Kokkos::parallel_for(
        "vel", Range3D({0, 0, 0}, {3*numMatrices, iDeg, iDeg}),
        KOKKOS_LAMBDA(int iMatCart, int i, int j) {
          int iMat = iMatCart / 3;
          int iCart = iMatCart % 3;
          int iK = funcK(iMat);
          int ib = funcB(iMat);
          subVelocity(iMatCart, i, j) = 0.5
              * (resultVelocities(iK, ib + i, ib + j, iCart)
                 + Kokkos::conj(resultVelocities(iK, ib + j, ib + i, iCart)));
        });

    // now diagonalize the blocks
    ComplexView3D newEigenvectors("newEig", 3*numMatrices, iDeg, iDeg);
    Kokkos::deep_copy(newEigenvectors, subVelocity);
    // Kokkos::parallel_for(
    //     "vel", Range3D({0, 0, 0}, {3*numMatrices, iDeg, iDeg}),
    //     KOKKOS_LAMBDA(int iMatCart, int i, int j) {
    //       newEigenvectors(iMatCart, i, j) = subVelocity(iMatCart, i, j);
    //     });
    // I don't need the eigenvalues, but I need the allocation to run zheev
    DoubleView2D newEigenvalues("newEigv", 3*numMatrices, iDeg);
    kokkosZHEEV(newEigenvectors, newEigenvalues);

    // rotate the original matrix in the new basis
    // that diagonalizes the subspace.
    // subMat = newEigenvectors.adjoint() * subMat * newEigenvectors;
    ComplexView3D tmpSubMat("tmpSubMat", 3*numMatrices, iDeg, iDeg);
    Kokkos::parallel_for(
        "vel", Range3D({0, 0, 0}, {3*numMatrices, iDeg, iDeg}),
        KOKKOS_LAMBDA(int iMatCart, int i, int j) {
          Kokkos::complex<double> tmp(0., 0.);
          for (int k = 0; k < iDeg; k++) {
            tmp += Kokkos::conj(newEigenvectors(iMatCart, k, i)) * subVelocity(iMatCart, k, j);
          }
          tmpSubMat(iMatCart, i, j) = tmp;
        });

    // finish the rotation and substitute back in the final results
    Kokkos::parallel_for(
        "vel", Range3D({0, 0, 0}, {3*numMatrices, iDeg, iDeg}),
        KOKKOS_LAMBDA(int iMatCart, int i, int j) {
          int iMat = iMatCart / 3;
          int iCart = iMatCart % 3;
          int iK = funcK(iMat);
          int ib = funcB(iMat);
          Kokkos::complex<double> tmp = zero;
          for (int k = 0; k < iDeg; ++k) {
            tmp += tmpSubMat(iMatCart, i, k) * newEigenvectors(iMatCart, k, j);
          }
          resultVelocities(iK, ib + i, ib + j, iCart) = tmp;
        });
  }

  return {resultEnergies, resultEigenvectors, resultVelocities};
}

