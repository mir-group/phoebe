#include "phonon_h0.h"

StridedComplexView3D PhononH0::kokkosBatchedBuildBlochHamiltonian(
    const DoubleView2D &cartesianCoordinates) {

  int numK = cartesianCoordinates.extent(0);
  int numBands = this->numBands;
  int numBravaisVectors = this->numBravaisVectors;
  Kokkos::complex<double> complexI(0.0, 1.0);

  auto bravaisVectors_d = this->bravaisVectors_d;
  auto weights_d = this->weights_d;
  auto mat2R_d = this->mat2R_d;
  auto numAtoms = this->numAtoms;

  ComplexView3D dynamicalMatrices("dynMat", numK, numBands, numBands);

  ComplexView2D phases_d("phPhases_d", numK, numBravaisVectors);
  Kokkos::parallel_for(
      "el_phases", Range2D({0, 0}, {numK, numBravaisVectors}),
      KOKKOS_LAMBDA(int iK, int iR) {
        double arg = 0.0;
        for (int i = 0; i < 3; i++) {
          arg += cartesianCoordinates(iK, i) * bravaisVectors_d(iR, i);
        }
        phases_d(iK, iR) = exp(-complexI * arg) * weights_d(iR);
      });
  Kokkos::fence();
  //print2DComplex("new = ", phases_d);
  //print2D("new = ", bravaisVectors_d);
  //print3D("new = ", mat2R_d);

  Kokkos::parallel_for(
      "el_hamilton", Range3D({0, 0, 0}, {numK, numBands, numBands}),
      KOKKOS_LAMBDA(int iK, int m, int n) {
        Kokkos::complex<double> tmp(0.0);
        for (int iR = 0; iR < numBravaisVectors; iR++) {
          tmp += phases_d(iK, iR) * mat2R_d(m, n, iR);
        }
        dynamicalMatrices(iK, m, n) = tmp;
      });
  Kokkos::realloc(phases_d, 0, 0);
  //print3DComplex("new = ", dynamicalMatrices);

  if (hasDielectric) {
    auto longRangeCorrection1_d = this->longRangeCorrection1_d;

    Kokkos::parallel_for(
        "el_hamilton", Range3D({0, 0, 0}, {numK, 3, 3}),
        KOKKOS_LAMBDA(int iK, int i, int j) {
          auto D = Kokkos::subview(dynamicalMatrices, iK, Kokkos::ALL, Kokkos::ALL);
          for (int iAt=0; iAt<numAtoms; ++iAt) {
            int m = iAt * 3 + i;
            int n = iAt * 3 + j;
            D(m, n) += longRangeCorrection1_d(i, j, iAt);
          }
        });
    Kokkos::fence();

    double norm = e2 * fourPi / volumeUnitCell;
    int numG = gVectors_d.extent(0);

    auto gVectors_d = this->gVectors_d;
    auto dielectricMatrix_d = this->dielectricMatrix_d;
    auto bornCharges_d = this->bornCharges_d;
    auto atomicPositions_d = this->atomicPositions_d;
    auto gMax = this->gMax;

    int byte_size = Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace::scratch_memory_space>::shmem_size(numBands);
    Kokkos::parallel_for("long-range-ph-H0",
        Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>(numG, Kokkos::AUTO()).set_scratch_size(
          0, Kokkos::PerThread(byte_size)),
        KOKKOS_LAMBDA (const Kokkos::TeamPolicy<>::member_type &team) {

          int iG = team.league_rank();
          auto G = Kokkos::subview(gVectors_d, iG, Kokkos::ALL);

          Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace::scratch_memory_space> GQZ(team.thread_scratch(0), 3*numAtoms);

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, numK), [&] (const int iK){

            auto Q = Kokkos::subview(cartesianCoordinates, iK, Kokkos::ALL);
            auto D = Kokkos::subview(dynamicalMatrices, iK, Kokkos::ALL, Kokkos::ALL);

            double GQ[3];
            for (int i=0; i<3; ++i) {
              GQ[i] = G(i) + Q(i);
            }

            double geg = 0.;
            Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, 9), [&] (const int ij, double &geg){
                int i = ij / 3;
                int j = ij % 3;
                geg += GQ[i] * dielectricMatrix_d(i, j) * GQ[j];
            }, geg);

            if (geg > 0. && geg < 4. * gMax) {
              double normG = norm * exp(-geg * 0.25) / geg;

              Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, 3*numAtoms), [&] (const int inb){
                  int nb = inb / 3;
                  int i = inb % 3;
                  double tmp = 0.0;
                  for (int j : {0, 1, 2}) {
                     tmp += GQ[j] * bornCharges_d(i, j, nb);
                  }
                  GQZ[inb] = tmp;
              });

              Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, numAtoms*numAtoms), [&] (const int nanb){
                  int na = nanb / numAtoms;
                  int nb = nanb % numAtoms;
                  double arg = 0.;
                  for (int i = 0; i < 3; i++) {
                    double xi = atomicPositions_d(na,i) - atomicPositions_d(nb,i);
                    arg += xi * GQ[i];
                  }
                  Kokkos::complex<double> phase = normG * exp(complexI * arg);
                  for (int j : {0, 1, 2}) {
                    for (int i : {0, 1, 2}) {
                      Kokkos::complex<double> x = phase * GQZ[na*3+i] * GQZ[nb*3+j];
                      Kokkos::atomic_add(&D(na*3+i, nb*3+j), x);
                      // there's a reduce on iG!
                    }
                  }
              });
            }
          });
        });
  //print3DComplex("new = ", dynamicalMatrices);
  }


  // Col-major matrices stored contiguously
  Kokkos::LayoutStride DDlayout(
      numK, numBands*numBands,
      numBands, 1,
      numBands, numBands
  );

  StridedComplexView3D DD("dynMat", DDlayout);
  Kokkos::parallel_for(
      "el_hamilton", Range3D({0, 0, 0}, {numK, numBands, numBands}),
      KOKKOS_LAMBDA(int iK, int m, int n) {
        auto D = Kokkos::subview(dynamicalMatrices, iK, Kokkos::ALL, Kokkos::ALL);
        DD(iK, m, n) = 0.5 * ( D(m,n) + Kokkos::conj(D(n,m)) );
      });
  Kokkos::realloc(dynamicalMatrices, 0, 0, 0);

  auto atomicMasses_d = this->atomicMasses_d;

  Kokkos::parallel_for(
      "el_hamilton", Range3D({0, 0, 0}, {numK, numBands, numBands}),
      KOKKOS_LAMBDA(int iK, int m, int n) {
        DD(iK, m, n) /= sqrt(atomicMasses_d(m) * atomicMasses_d(n));
      });

  //print3DComplex("new = ", DD);
  return DD;
}

std::tuple<DoubleView2D, StridedComplexView3D> PhononH0::kokkosBatchedDiagonalizeFromCoordinates(
    const DoubleView2D &cartesianCoordinates, const bool withMassScaling) {

  StridedComplexView3D dynamicalMatrices =
      kokkosBatchedBuildBlochHamiltonian(cartesianCoordinates);

  // now the diagonalization.

  int numK = dynamicalMatrices.extent(0);
  int numBands = this->numBands;

  //print3DComplex("new = ", dynamicalMatrices);
  DoubleView2D frequencies("energies", numK, numBands);
  kokkosZHEEV(dynamicalMatrices, frequencies);

  Kokkos::parallel_for(
      "el_hamilton", Range2D({0, 0}, {numK, numBands}),
      KOKKOS_LAMBDA(int iK, int m) {
        if (frequencies(iK, m) > 0.) {
          frequencies(iK, m) = sqrt(frequencies(iK, m));
        } else {
          frequencies(iK, m) = -sqrt(-frequencies(iK, m));
        }
      });
  //print2D("new = ", frequencies);
  //print3DComplex("new = ", dynamicalMatrices);
  if(withMassScaling) kokkosBatchedScaleEigenvectors(dynamicalMatrices);
  return std::make_tuple(frequencies, dynamicalMatrices);
}

void PhononH0::kokkosBatchedScaleEigenvectors(StridedComplexView3D& eigenvectors) {
  int numK = eigenvectors.extent(0);
  int numBands = this->numBands;
  auto atomicMasses_d = this->atomicMasses_d;

  Kokkos::parallel_for(
      "mass_rescale", Range3D({0, 0, 0}, {numK, numBands, numBands}),
      KOKKOS_LAMBDA(int iK, int m, int n) {
        eigenvectors(iK, m, n) /= sqrt(atomicMasses_d(m));
      });
}

double PhononH0::getDeviceMemoryUsage() {
  double memory = 16 * longRangeCorrection1_d.size()
      + 8 * (atomicMasses_d.size() + gVectors_d.size()
           + dielectricMatrix_d.size() + bornCharges_d.size()
           + atomicPositions_d.size() + bravaisVectors_d.size()
           + weights_d.size() + mat2R_d.size() );
  return memory;
}

int PhononH0::estimateBatchSize(const bool& withVelocity) {
  double memoryAvailable = kokkosDeviceMemory->getAvailableMemory();

  int matrixSize = numBands * numBands;
  int transformSize = matrixSize + numBravaisVectors;
  transformSize = std::max(transformSize, 2*matrixSize + numBands);

  double memoryPerPoint = 0.;
  if (withVelocity) {
    memoryPerPoint += 16. * (transformSize*7 + 2*matrixSize + matrixSize*3);
  } else {
    memoryPerPoint += 16. * transformSize;
  }
  //printf("memory per point = %g, available = %g\n", memoryPerPoint, memoryAvailable);

  // we try to use 90% of the available memory (leave some buffer)
  int numBatches = int(memoryAvailable / memoryPerPoint * 0.95);
  if (numBatches < 1) {
    Error("Not enough memory available on device, try reduce memory usage");
  }

  return numBatches;
}

FullBandStructure PhononH0::kokkosPopulate(Points &fullPoints,
                                           const bool &withVelocities,
                                           const bool &withEigenvectors,
                                           const bool isDistributed) {

  int numBands = getNumBands();
  auto particle = getParticle();

  FullBandStructure fullBandStructure(numBands, particle, withVelocities,
                                      withEigenvectors, fullPoints,
                                      isDistributed);

  std::vector<std::vector<int>> ikBatches;
  {
    std::vector<int> ikIterator = fullBandStructure.getWavevectorIndices();
    int batchSize = estimateBatchSize(withVelocities);
    ikBatches = kokkosDeviceMemory->splitToBatches(ikIterator, batchSize);
  }

  for (auto ikBatch : ikBatches) {

    int numK = ikBatch.size();

    DoubleView2D cartesianWavevectors_d("el_cartWav_d", numK, 3);
    {
      auto cartesianWavevectors_h = Kokkos::create_mirror_view(cartesianWavevectors_d);
#pragma omp parallel for
      for (int iik = 0; iik < numK; ++iik) {
        int ik = ikBatch[iik];
        WavevectorIndex ikIdx(ik);
        Eigen::Vector3d k = fullBandStructure.getWavevector(ikIdx);
        for (int i = 0; i < 3; ++i) {
          cartesianWavevectors_h(iik, i) = k(i);
        }
      }
      Kokkos::deep_copy(cartesianWavevectors_d, cartesianWavevectors_h);
    }

    {
      DoubleView2D allEnergies_d;
      StridedComplexView3D allEigenvectors_d;
      ComplexView4D allVelocities_d;
      if (withVelocities) {
        auto t = kokkosBatchedDiagonalizeWithVelocities(cartesianWavevectors_d);
        allEnergies_d = std::get<0>(t);
        allEigenvectors_d = std::get<1>(t);
        allVelocities_d = std::get<2>(t);
      } else {
        auto t = kokkosBatchedDiagonalizeFromCoordinates(cartesianWavevectors_d);
        allEnergies_d = std::get<0>(t);
        allEigenvectors_d = std::get<1>(t);
      }
      Kokkos::realloc(cartesianWavevectors_d, 0, 0);

      // TODO: this is the only difference from the ElectronH0 kokkosPopulate()
      // can we unify the two?
      kokkosBatchedScaleEigenvectors(allEigenvectors_d);

      auto allEnergies_h = Kokkos::create_mirror_view(allEnergies_d);
      auto allEigenvectors_h = Kokkos::create_mirror_view(allEigenvectors_d);
      Kokkos::deep_copy(allEnergies_h, allEnergies_d);
      Kokkos::deep_copy(allEigenvectors_h, allEigenvectors_d);

#pragma omp parallel for
      for (int iik = 0; iik < numK; iik++) {
        int ik = ikBatch[iik];
        Point point = fullBandStructure.getPoint(ik);
        Eigen::VectorXd ens(numBands);
        for (int ib = 0; ib < numBands; ++ib) {
          ens(ib) = allEnergies_h(iik, ib);
        }
        fullBandStructure.setEnergies(point, ens);
        if (withEigenvectors) {
          Eigen::MatrixXcd eigVec(numBands, numBands);
          for (int ib1 = 0; ib1 < numBands; ++ib1) {
            for (int ib2 = 0; ib2 < numBands; ++ib2) {
              eigVec(ib1, ib2) = allEigenvectors_h(iik, ib1, ib2);
            }
          }
          fullBandStructure.setEigenvectors(point, eigVec);
        }
      }

      if (withVelocities) {
        auto allVelocities_h = Kokkos::create_mirror_view(allVelocities_d);
        Kokkos::deep_copy(allVelocities_h, allVelocities_d);
        for (int iik = 0; iik < numK; iik++) {
          int ik = ikBatch[iik];
          Point point = fullBandStructure.getPoint(ik);
          Eigen::Tensor<std::complex<double>,3> v(numBands, numBands, 3);
          for (int ib1 = 0; ib1 < numBands; ++ib1) {
            for (int ib2 = 0; ib2 < numBands; ++ib2) {
              for (int i = 0; i < 3; ++i) {
                v(ib1, ib2, i) = allVelocities_h(iik, ib1, ib2, i);
              }
            }
          }
          fullBandStructure.setVelocities(point, v);
        }
      }
      Kokkos::realloc(allEnergies_d, 0, 0);
      allEigenvectors_d = decltype(allEigenvectors_d)();
      Kokkos::realloc(allVelocities_d, 0, 0, 0, 0);
    }
  }

  return fullBandStructure;
}


std::tuple<DoubleView2D, StridedComplexView3D, ComplexView4D>
PhononH0::kokkosBatchedDiagonalizeWithVelocities(
    const DoubleView2D &cartesianCoordinates) {

  // Note: this is slightly different than electronH0Wannier
  // here, we need to compute the derivative of sqrt(DynamicalMatrix)
  // while for electrons we derive the BlochHamiltonian directly

  int numK = cartesianCoordinates.extent(0);

  double delta = 1.0e-8;

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

  //print2D("new = ", allVectors);

  // compute the electronic properties at all wavevectors
  auto t = kokkosBatchedDiagonalizeFromCoordinates(allVectors, false);
  DoubleView2D allEnergies = std::get<0>(t);
  StridedComplexView3D allEigenvectors = std::get<1>(t);

  int numBands = allEnergies.extent(1);

  //print2D("new = ", allEnergies);
  //print3DComplex("new = ", allEigenvectors);

  Kokkos::LayoutStride eigenvectorLayout(
      numK, numBands*numBands,
      numBands, 1,
      numBands, numBands);

  // save energies and eigenvectors to results
  DoubleView2D resultEnergies("energies", numK, numBands);
  StridedComplexView3D resultEigenvectors("eigenvectors", eigenvectorLayout);
  ComplexView4D resultVelocities("velocities", numK, numBands, numBands, 3);

  Kokkos::parallel_for(
      "der", Range2D({0, 0}, {numK, numBands}), KOKKOS_LAMBDA(int iK, int m) {
        auto E = Kokkos::subview(allEnergies, iK * 7, Kokkos::ALL);
        auto X = Kokkos::subview(allEigenvectors, iK * 7, Kokkos::ALL, Kokkos::ALL);
        resultEnergies(iK, m) = E(m);
        for (int n=0; n<numBands; ++n) {
          resultEigenvectors(iK, m, n) = X(m, n);
        }
      });

  // these are temporary "scratch" memory spaces
  ComplexView3D der("der", numK, numBands, numBands);
  ComplexView3D tmpV("tmpV", numK, numBands, numBands);
  //ComplexView3D plus("plus", numK, numBands, numBands);
  //ComplexView3D minus("plus", numK, numBands, numBands);

  for (int i = 0; i < 3; ++i) {

    // To build the velocity operator, we first compute the matrix
    // der = dH/dk = ( H(k+dk)-H(k-dk) ) / 2dk
    // we also reinforce the Hermiticity (better numerical properties)
    Kokkos::parallel_for(
        "der", Range3D({0, 0, 0}, {numK, numBands, numBands}),
        KOKKOS_LAMBDA(int iK, int m, int n) {
          auto XPlus = Kokkos::subview(allEigenvectors, iK * 7 + i * 2 + 1, Kokkos::ALL, Kokkos::ALL);
          auto XMins = Kokkos::subview(allEigenvectors, iK * 7 + i * 2 + 2, Kokkos::ALL, Kokkos::ALL);
          DoubleView1D EPlus = Kokkos::subview(allEnergies, iK * 7 + i * 2 + 1, Kokkos::ALL);
          DoubleView1D EMins = Kokkos::subview(allEnergies, iK * 7 + i * 2 + 2, Kokkos::ALL);
          Kokkos::complex<double> x(0.,0.);
          //Kokkos::complex<double> p(0.,0.);
          //Kokkos::complex<double> pm(0.,0.);
          for (int l=0; l<numBands; ++l) {
            x += XPlus(m,l) * EPlus(l) * Kokkos::conj(XPlus(n,l))
                - XMins(m,l) * EMins(l) * Kokkos::conj(XMins(n,l));
            //p += XPlus(m,l) * EPlus(l) * Kokkos::conj(XPlus(n,l));
            //pm += XMins(m,l) * EMins(l) * Kokkos::conj(XMins(n,l));
          }
          der(iK, m, n) = x/(2*delta);
          //plus(iK,m,n) = p;
          //minus(iK,m,n) = pm;
        });
    Kokkos::fence();
    //print3DComplex("new = ", der);
    //print3DComplex("new = ", plus);
    //print3DComplex("new = ", minus);

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
            tmp += Kokkos::conj(L(l,m)) * 0.5 * (R(l, n) + Kokkos::conj(R(n, l)));
          }
          A(m, n) = tmp;
        });
    Kokkos::fence();

    Kokkos::parallel_for(
        "vel", Range3D({0, 0, 0}, {numK, numBands, numBands}),
        KOKKOS_LAMBDA(int iK, int m, int n) {
          double norm = 0.;
          for (int i=0; i<3; ++i) {
            norm += cartesianCoordinates(iK,i) * cartesianCoordinates(iK,i);
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

  //print4DComplex("new = ", resultVelocities);

  kokkosBatchedTreatDegenerateVelocities(cartesianCoordinates, resultEnergies,
                                         resultVelocities, 0.0001 / ryToCmm1);

  //print4DComplex("new = ", resultVelocities);

//  // the degeneracy diagonalization is a bit tricky.
//  // in the way it was implemented for the CPU case can't work for GPU
//  // that's because it's hard to do several diagonalizations on GPU
//  // (due to the fact that the GPU can't launch kernels by itself)
//
//  // ALGO:
//  // - find all the degrees of degeneracy that have to be run.
//  // - loop over the degree of degeneracy
//  //   - set up all sub-matrices of velocities of rank (ndegree)
//  //   - diagonalize them
//  //   - rotate velocity elements
//
//  // degCount
//  IntView2D degCount("tmpV", numK, numBands);
//  // maxDegeneracy is the largest degenerate block of energies
//  int maxDegeneracy = 1;
//  // here I do two things: find the maxDegeneracy for the current system
//  // and set degCount such that, degCount(ik,ib) = 0 if it is a state that is
//  // not the starting point of a degenerate block, or degCount(ik,ib) = iDeg
//  // if ik,ib is at the start of a block iDegxiDeg that needs to be
//  // diagonalized, with iDeg >= 2
//  Kokkos::parallel_reduce(
//      "Loop1", numK,
//      KOKKOS_LAMBDA(const int &iK, int &iDeg) {
//        for (int ib = 0; ib < numBands; ib++) {
//          // first, we check if the band is degenerate, and the size of the
//          // degenerate subspace
//          int sizeSubspace = 1;
//          for (int ib2 = ib + 1; ib2 < numBands; ib2++) {
//            // I consider bands degenerate if their frequencies are the same
//            // within 0.0001 cm^-1
//            if (abs(resultEnergies(iK, ib) - resultEnergies(iK, ib2)) > threshold) {
//              break;
//            }
//            ++sizeSubspace;
//          }
//
//          if (sizeSubspace > iDeg) {
//            iDeg = sizeSubspace; // iDeg = std::max(iDeg, sizeSubspace);
//          }
//          if (sizeSubspace == 1) {
//            degCount(iK, ib) = 0;
//          } else {
//            degCount(iK, ib) = sizeSubspace;
//            for (int i=1; i<sizeSubspace; ++i) {
//              degCount(iK, ib+i) = 0;
//            }
//          }
//
//          // we skip the bands in the subspace, since we corrected them already
//          ib += sizeSubspace - 1;
//        }
//      },
//      Kokkos::Max<int>(maxDegeneracy));
//
//  // we now do a diagonalization of all degenerate velocity sub-blocks
//  // since cuda can only launch the diagonalization of N matrices at fixed
//  // matrix size, we now do a loop over the degree of degeneracy
//  for (int iDeg = 2; iDeg<=maxDegeneracy; ++iDeg) {
//
//    // count how many degenerate points we have at this degree of degeneracy
//    int numMatrices = 0;
//    Kokkos::parallel_reduce("Loop1", numK,
//        KOKKOS_LAMBDA(const int &iK, int &iMat) {
//          for (int ib = 0; ib < numBands; ib++) {
//            if (degCount(iK, ib) == iDeg) {
//              ++iMat;
//            }
//          }
//        }, numMatrices);
//
//    /**
//     * The next part is a bit convoluted because it's tricky to do in parallel
//     * We want to build N matrices of size iDegxiDeg containing the velocity
//     * degenerate blocks.
//     * To this aim, we need two functions funcK, funcB. Such that funcK(iMat),
//     * funcB(iMat) return the iK and ib index locating the start of the block
//     * that will be diagonalized by the iMat-th matrix
//     *
//     * We take advantage of an exclusive scan. Say we have a vector
//     * 0 0 0 1 0 0 0 0 1;
//     * which spans the bloch states, and 1 locates the beginning of a block
//     * of size iDeg that we want to diagonalize. The indices locating the "1",
//     * i.e. [3, 8], are given by an exclusive scan:
//     * 0 0 0 0 1 1 1 1 1 2
//     * So, scan(ik,ib) = iMat. And from this, I can do the inverse mapping.
//     */
//
//    IntView1D scan("scan", numK*numBands);
//    Kokkos::parallel_for(
//        Range2D({0,0},{numK, numBands}), KOKKOS_LAMBDA(int ik, int ib) {
//          if ( iDeg == degCount(ik,ib) ) {
//            scan(ik*numBands+ib) = 1;
//          } else {
//            scan(ik*numBands+ib) = 0;
//          }
//        });
//
//    Kokkos::parallel_scan("scan", numK*numBands, KOKKOS_LAMBDA(const int i,
//                                               int& update, const bool final) {
//          // Load old value in case we update it before accumulating
//          const int val_i = scan(i);
//          if (final) {
//            scan(i) = update; // only update array on final pass
//          }
//          // For exclusive scan, change the update value after
//          // updating array, like we do here. For inclusive scan,
//          // change the update value before updating array.
//          update += val_i;
//        });
//
//    IntView1D funcK("funcK", numMatrices);
//    IntView1D funcB("funcB", numMatrices);
//    Kokkos::parallel_for(
//        "funcKB", Range2D({0,0},{numK, numBands}), KOKKOS_LAMBDA(int ik, int ib) {
//          int ikb = ik*numBands + ib;
//          int iMat = scan(ikb);
//          if (degCount(ik, ib) == iDeg) { // if we actually had a 1
//            // the if is needed to avoid race conditions in funcK/funcB
//            funcK(iMat) = ik;
//            funcB(iMat) = ib;
//          }
//        });
//    Kokkos::resize(scan, 0);
//
//    // now, I copy the degenerate blocks of the velocity operator in a scratch
//    // memory space "subVelocity".
//    ComplexView3D subVelocity("subVel", 3*numMatrices, iDeg, iDeg);
//    Kokkos::parallel_for(
//        "vel", Range3D({0, 0, 0}, {3*numMatrices, iDeg, iDeg}),
//        KOKKOS_LAMBDA(int iMatCart, int i, int j) {
//          int iMat = iMatCart / 3;
//          int iCart = iMatCart % 3;
//          int iK = funcK(iMat);
//          int ib = funcB(iMat);
//          subVelocity(iMatCart, i, j) = 0.5
//              * (resultVelocities(iK, ib + i, ib + j, iCart)
//                 + Kokkos::conj(resultVelocities(iK, ib + j, ib + i, iCart)));
//        });
//
//    // now diagonalize the blocks
//    ComplexView3D newEigenvectors("newEig", 3*numMatrices, iDeg, iDeg);
//    Kokkos::deep_copy(newEigenvectors, subVelocity);
//    // I don't need the eigenvalues, but I need the allocation to run zheev
//    DoubleView2D newEigenvalues("newEigv", 3*numMatrices, iDeg);
//    kokkosZHEEV(newEigenvectors, newEigenvalues);
//    Kokkos::realloc(newEigenvalues, 0, 0);
//
//    // rotate the original matrix in the new basis
//    // that diagonalizes the subspace.
//    // subMat = newEigenvectors.adjoint() * subMat * newEigenvectors;
//    ComplexView3D tmpSubMat("tmpSubMat", 3*numMatrices, iDeg, iDeg);
//    Kokkos::parallel_for(
//        "vel", Range3D({0, 0, 0}, {3*numMatrices, iDeg, iDeg}),
//        KOKKOS_LAMBDA(int iMatCart, int i, int j) {
//          Kokkos::complex<double> tmp(0., 0.);
//          for (int k = 0; k < iDeg; k++) {
//            tmp += Kokkos::conj(newEigenvectors(iMatCart, k, i)) * subVelocity(iMatCart, k, j);
//          }
//          tmpSubMat(iMatCart, i, j) = tmp;
//        });
//
//    // finish the rotation and substitute back in the final results
//    Kokkos::parallel_for(
//        "vel", Range3D({0, 0, 0}, {3*numMatrices, iDeg, iDeg}),
//        KOKKOS_LAMBDA(int iMatCart, int i, int j) {
//          int iMat = iMatCart / 3;
//          int iCart = iMatCart % 3;
//          int iK = funcK(iMat);
//          int ib = funcB(iMat);
//          Kokkos::complex<double> tmp = zero;
//          for (int k = 0; k < iDeg; ++k) {
//            tmp += tmpSubMat(iMatCart, i, k) * newEigenvectors(iMatCart, k, j);
//          }
//          resultVelocities(iK, ib + i, ib + j, iCart) = tmp;
//        });
//  }

  return std::make_tuple(resultEnergies, resultEigenvectors, resultVelocities);
}

