#include "harmonic.h"

void HarmonicHamiltonian::kokkosBatchedTreatDegenerateVelocities(
    const DoubleView2D& cartesianCoordinates,
    const DoubleView2D& resultEnergies, ComplexView4D& resultVelocities,
    const double& threshold) {

  int numK = cartesianCoordinates.extent(0);
  int numBands = resultEnergies.extent(1);

//  double threshold = 0.000001 / energyRyToEv;// = 1 micro-eV

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
    // I don't need the eigenvalues, but I need the allocation to run zheev
    DoubleView2D newEigenvalues("newEigv", 3*numMatrices, iDeg);
    kokkosZHEEV(newEigenvectors, newEigenvalues);
    Kokkos::realloc(newEigenvalues, 0, 0);

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
          Kokkos::complex<double> tmp(0.0,0.0);
          for (int k = 0; k < iDeg; ++k) {
            tmp += tmpSubMat(iMatCart, i, k) * newEigenvectors(iMatCart, k, j);
          }
          resultVelocities(iK, ib + i, ib + j, iCart) = tmp;
        });
  }

//  return {resultEnergies, resultEigenvectors, resultVelocities};
}

