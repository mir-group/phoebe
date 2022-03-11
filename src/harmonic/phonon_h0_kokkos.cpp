#include "phonon_h0.h"

ComplexView3D PhononH0::kokkosBatchedBuildBlochHamiltonian(
    const DoubleView2D &cartesianCoordinates) {

  int numK = cartesianCoordinates.extent(0);
  int numBands = this->numBands;
  int numBravaisVectors = this->numBravaisVectors;
  Kokkos::complex<double> complexI(0.0, 1.0);

  ComplexView3D dynamicalMatrices("dynMat", numK, numBands, numBands);

  ComplexView2D phases_d("elPhases_d", numK, numBravaisVectors);
  Kokkos::parallel_for(
      "el_phases", Range2D({0, 0}, {numK, numBravaisVectors}),
      KOKKOS_LAMBDA(int iK, int iR) {
        double arg = 0.0;
        for (int i = 0; i < 3; i++) {
          arg += cartesianCoordinates(iK, i) * bravaisVectors_d(iR, i);
        }
        phases_d(iK, iR) = exp(-complexI * arg) * weights_d(iR);
      });

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

  if (hasDielectric) {

    Kokkos::parallel_for(
        "el_hamilton", Range3D({0, 0, 0}, {numK, 3, 3}),
        KOKKOS_LAMBDA(int iK, int i, int j) {
          auto D = Kokkos::subview(dynamicalMatrices, iK, Kokkos::ALL, Kokkos::ALL);
          for (int iAt=0; iAt<numAtoms; ++iAt) {
            int m = iAt * 3 + i;
            int n = iAt * 3 + j;
            D(m, n) += longRangeCorrection1_d(m, n, iAt);
          }
        });

    double norm = e2 * fourPi / volumeUnitCell;
    int numG = gVectors_d.extent(0);

    Kokkos::parallel_for(
        "el_hamilton", Range2D({0,0},{numK, numG}), KOKKOS_LAMBDA(int iK, int iG) {
          auto G = Kokkos::subview(gVectors_d, iG, Kokkos::ALL);
          auto Q = Kokkos::subview(cartesianCoordinates, iK, Kokkos::ALL);
          auto D = Kokkos::subview(dynamicalMatrices, iK, Kokkos::ALL, Kokkos::ALL);

          double GQ[3];
          for (int i=0; i<3; ++i) {
            GQ[i] = G(i) + Q(i);
          }

          double geg = 0.;
          for (int i=0; i<3; ++i) {
            for (int j = 0; j < 3; ++j) {
              geg += GQ[i] * dielectricMatrix_d(i, j) * GQ[j];
            }
          }

          if (geg > 0. && geg < 4. * gMax) {
            double normG = norm * exp(-geg * 0.25) / geg;

            //TODO: this instruction may not work on a GPU. Fix this!
            std::vector<double> GQZ(3*numAtoms, 0.);
            for (int i : {0, 1, 2}) {
              for (int j : {0, 1, 2}) {
                for (int nb = 0; nb < numAtoms; nb++) {
                  GQZ[nb * 3 + i] = GQ[j] * bornCharges_d(i, j, nb);
                }
              }
            }

            for (int nb = 0; nb < numAtoms; nb++) {
              for (int na = 0; na < numAtoms; na++) {
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
                  }
                }
              }
            }
          }

        });

  }

  ComplexView3D DD("dynMat", numK, numBands, numBands);
  Kokkos::parallel_for(
      "el_hamilton", Range3D({0, 0, 0}, {numK, numBands, numBands}),
      KOKKOS_LAMBDA(int iK, int m, int n) {
        auto D = Kokkos::subview(dynamicalMatrices, iK, Kokkos::ALL, Kokkos::ALL);
        DD(iK, m, n) = 0.5 * ( D(m,n) + Kokkos::conj(D(n,m)) );
      });
  Kokkos::realloc(dynamicalMatrices, 0, 0, 0);

  Kokkos::parallel_for(
      "el_hamilton", Range3D({0, 0, 0}, {numK, numBands, numBands}),
      KOKKOS_LAMBDA(int iK, int m, int n) {
        auto D = Kokkos::subview(DD, iK, Kokkos::ALL, Kokkos::ALL);
        D(m, n) /= sqrt(atomicMasses_d(m) * atomicMasses_d(n));
      });

  return DD;
}

std::tuple<DoubleView2D, ComplexView3D> PhononH0::kokkosBatchedDiagonalizeFromCoordinates(
    const DoubleView2D &cartesianCoordinates) {

  ComplexView3D dynamicalMatrices =
      kokkosBatchedBuildBlochHamiltonian(cartesianCoordinates);

  // now the diagonalization.

  int numK = dynamicalMatrices.extent(0);
  int numBands = this->numBands;

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
  return {frequencies, dynamicalMatrices};
}

void PhononH0::kokkosBatchedScaleEigenvectors(ComplexView3D& eigenvectors) {
  int numK = eigenvectors.extent(0);
  int numBands = this->numBands;

  Kokkos::parallel_for(
      "mass_rescale", Range3D({0, 0, 0}, {numK, numBands, numBands}),
      KOKKOS_LAMBDA(int iK, int m, int n) {
        auto D = Kokkos::subview(eigenvectors, iK, Kokkos::ALL, Kokkos::ALL);
        D(m, n) /= sqrt(atomicMasses_d(m));
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

  // we try to use 90% of the available memory (leave some buffer)
  int numBatches = int(memoryAvailable / memoryPerPoint * 0.95);
  if (numBatches < 1) {
    Error("Not enough memory available on device, try reduce memory usage");
  }

  return numBatches;
}

std::tuple<DoubleView2D, ComplexView3D, ComplexView4D>
PhononH0::kokkosBatchedDiagonalizeWithVelocities(
    const DoubleView2D &cartesianCoordinates) {
  auto t = HarmonicHamiltonian::kokkosBatchedDiagonalizeWithVelocities(cartesianCoordinates);
  auto w = std::get<0>(t);
  auto e = std::get<1>(t);
  auto v = std::get<2>(t);
  kokkosBatchedScaleEigenvectors(e);
  return {w,e,v};
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
      ComplexView3D allEigenvectors_d;
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
      Kokkos::realloc(allEigenvectors_d, 0, 0, 0);
      Kokkos::realloc(allVelocities_d, 0, 0, 0, 0);
    }
  }

  return fullBandStructure;
}
