#include "electron_h0_wannier.h"
#include "constants.h"
#include "exceptions.h"

ElectronH0Wannier::ElectronH0Wannier(
    const Eigen::Matrix3d &directUnitCell_,
    const Eigen::Matrix<double, 3, Eigen::Dynamic> &bravaisVectors_,
    const Eigen::VectorXd &vectorsDegeneracies_,
    const Eigen::Tensor<std::complex<double>, 3> &h0R_,
    const Eigen::Tensor<std::complex<double>, 4> &rMatrix_)
    : particle(Particle::electron) {

  h0R = h0R_;
  rMatrix = rMatrix_;
  directUnitCell = directUnitCell_;
  bravaisVectors = bravaisVectors_;
  vectorsDegeneracies = vectorsDegeneracies_;

  if (h0R.dimension(1) != h0R.dimension(2)) {
    Error("WannierH0(): h0R should have dimensions (R,bands,bands)");
  }
  if (h0R.dimension(0) != bravaisVectors.cols()) {
    Error("WannierH0(): h0R and bravaisVectors not aligned");
  }
  if (vectorsDegeneracies.size() != bravaisVectors.cols()) {
    Error("WannierH0(): degeneracies not aligned with vectors");
  }

  if ((rMatrix.dimension(1) != h0R.dimension(0)) ||
      (rMatrix.dimension(2) != h0R.dimension(1)) ||
      (rMatrix.dimension(3) != h0R.dimension(2))) {
    Error("WannierH0(): h0R and rMatrix should be aligned");
  }

  if (rMatrix.dimension(0) != 3) {
    Error("WannierH0(): rMatrix should be a vector");
  }

  numWannier = int(h0R.dimension(1));
  numVectors = int(vectorsDegeneracies.size());

  // copy data to the GPU
  {
    Kokkos::resize(h0R_d, numWannier, numWannier, numVectors);
    Kokkos::resize(vectorsDegeneracies_d, numVectors);
    Kokkos::resize(bravaisVectors_d, numVectors, 3);

    auto h0R_h = create_mirror_view(h0R_d);
    auto vectorsDegeneracies_h = create_mirror_view(vectorsDegeneracies_d);
    auto bravaisVectors_h = create_mirror_view(bravaisVectors_d);

    for (int iR = 0; iR < numVectors; iR++) {
      for (int i = 0; i < numWannier; i++) {
        for (int j = 0; j < numWannier; j++) {
          h0R_h(i, j, iR) = h0R(iR, i, j);
        }
      }
      vectorsDegeneracies_h(iR) = vectorsDegeneracies(iR);
      for (int i = 0; i < 3; ++i) {
        bravaisVectors_h(iR, i) = bravaisVectors(i, iR);
      }
    }
    Kokkos::deep_copy(h0R_d, h0R_h);
    Kokkos::deep_copy(vectorsDegeneracies_d, vectorsDegeneracies_h);
    Kokkos::deep_copy(bravaisVectors_d, bravaisVectors_h);
  }
}

// copy constructor
ElectronH0Wannier::ElectronH0Wannier(const ElectronH0Wannier &that)
    : particle(Particle::electron) {
  h0R = that.h0R;
  rMatrix = that.rMatrix;
  directUnitCell = that.directUnitCell;
  numWannier = that.numWannier;
  bravaisVectors = that.bravaisVectors;
  numVectors = that.numVectors;
  vectorsDegeneracies = that.vectorsDegeneracies;
  hasShiftedVectors = that.hasShiftedVectors;
  vectorsShifts = that.vectorsShifts;
  degeneracyShifts = that.degeneracyShifts;
  h0R_d = that.h0R_d;
  degeneracyShifts_d = that.degeneracyShifts_d;
  vectorsShifts_d = that.vectorsShifts_d;
  vectorsDegeneracies_d = that.vectorsDegeneracies_d;
  bravaisVectors_d = that.bravaisVectors_d;
}

// copy assignment
ElectronH0Wannier &ElectronH0Wannier::operator=(const ElectronH0Wannier &that) {
  if (this != &that) {
    particle = that.particle;
    numVectors = that.numVectors;
    numWannier = that.numWannier;
    bravaisVectors = that.bravaisVectors;
    vectorsDegeneracies = that.vectorsDegeneracies;
    directUnitCell = that.directUnitCell;
    h0R = that.h0R;
    rMatrix = that.rMatrix;
    hasShiftedVectors = that.hasShiftedVectors;
    vectorsShifts = that.vectorsShifts;
    degeneracyShifts = that.degeneracyShifts;
    h0R_d = that.h0R_d;
    degeneracyShifts_d = that.degeneracyShifts_d;
    vectorsShifts_d = that.vectorsShifts_d;
    vectorsDegeneracies_d = that.vectorsDegeneracies_d;
    bravaisVectors_d = that.bravaisVectors_d;
  }
  return *this;
}

ElectronH0Wannier::~ElectronH0Wannier() {
  // Deallocate stuff from GPU
  // Eigen class attributes should be deallocated automatically
  Kokkos::realloc(h0R_d, 0, 0, 0);
  Kokkos::realloc(degeneracyShifts_d, 0, 0, 0);
  Kokkos::realloc(vectorsShifts_d, 0, 0, 0, 0, 0);
  Kokkos::realloc(vectorsDegeneracies_d, 0);
  Kokkos::realloc(bravaisVectors_d, 0, 0);
}

Particle ElectronH0Wannier::getParticle() { return particle; }

int ElectronH0Wannier::getNumBands() { return numWannier; }

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
ElectronH0Wannier::diagonalize(Point &point) {
  Eigen::Vector3d k = point.getCoordinates(Points::cartesianCoordinates);

  auto tup = diagonalizeFromCoordinates(k);
  auto energies = std::get<0>(tup);
  auto eigenvectors = std::get<1>(tup);

  // note: the eigenvector matrix is the unitary transformation matrix U
  // from the Bloch to the Wannier gauge.

  return {energies, eigenvectors};
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
ElectronH0Wannier::diagonalizeFromCoordinates(Eigen::Vector3d &k) {

  std::vector<Eigen::Vector3d> cartesianWavevectors;
  cartesianWavevectors.push_back(k);
  auto t = batchedDiagonalizeFromCoordinates(cartesianWavevectors);
  auto energies = std::get<0>(t)[0];
  auto eigenVectors = std::get<1>(t)[0];
  return {energies, eigenVectors};
}

std::vector<Eigen::MatrixXcd> ElectronH0Wannier::batchedBuildHamiltonians(
    std::vector<Eigen::Vector3d>& cartesianWavevectors) {

  int numK = cartesianWavevectors.size();
  std::vector<Eigen::MatrixXcd> Hs(numK, Eigen::MatrixXcd::Zero(numWannier, numWannier));

  if (!hasShiftedVectors) {
    Eigen::MatrixXcd phases(numVectors, numK);
#pragma omp parallel for collapse(2)
    for (int iK = 0; iK < numK; ++iK) {
      for (int iR = 0; iR < numVectors; iR++) {
        double phase = cartesianWavevectors[iK].dot(bravaisVectors.col(iR));
        std::complex<double> phaseFactor = {cos(phase), sin(phase)};
        phases(iR, iK) = phaseFactor / vectorsDegeneracies(iR);
      }
    }
#pragma omp parallel for collapse(3)
    for (int iK = 0; iK < numK; ++iK) {
      for (int n = 0; n < numWannier; n++) {
        for (int m = 0; m < numWannier; m++) {
          for (int iR = 0; iR < numVectors; iR++) {
            Hs[iK](m, n) += phases(iR, iK) * h0R(iR, m, n);
          }
        }
      }
    }

  } else {

    for (int iK = 0; iK < numK; ++iK) {
      for (int iR = 0; iR < numVectors; iR++) {
        double phaseArg;
        std::complex<double> phase;
        for (int iw2 = 0; iw2 < numWannier; iw2++) {
          for (int iw1 = 0; iw1 < numWannier; iw1++) {
            for (int iDeg = 0; iDeg < degeneracyShifts(iw1, iw2, iR); ++iDeg) {
              phaseArg = 0.;
              for (int i : {0, 1, 2}) {
                phaseArg += cartesianWavevectors[iK](i) * vectorsShifts(i, iDeg, iw1, iw2, iR);
              }
              phase = {cos(phaseArg), sin(phaseArg)};
              phase /= vectorsDegeneracies(iR) * degeneracyShifts(iw1, iw2, iR);
              Hs[iK](iw1, iw2) += phase * h0R(iR, iw1, iw2);
            }
          }
        }
      }
    }
  }
  return Hs;
}

std::tuple<std::vector<Eigen::VectorXd>, std::vector<Eigen::MatrixXcd>>
ElectronH0Wannier::batchedDiagonalizeFromCoordinates(std::vector<Eigen::Vector3d>& cartesianWavevectors) {

  auto Hs = batchedBuildHamiltonians(cartesianWavevectors);

  int numK = cartesianWavevectors.size();

  std::vector<Eigen::VectorXd> allEnergies(numK);
  std::vector<Eigen::MatrixXcd> allEigenvectors(numK);
#pragma omp parallel for
  for (int iK = 0; iK < numK; ++iK) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigenSolver(Hs[iK]);
    allEnergies[iK] = eigenSolver.eigenvalues();
    allEigenvectors[iK] = eigenSolver.eigenvectors();
  }
  return {allEnergies, allEigenvectors};
}

Eigen::Tensor<std::complex<double>, 3>
ElectronH0Wannier::diagonalizeVelocity(Point &point) {
  Eigen::Vector3d coordinates = point.getCoordinates(Points::cartesianCoordinates);
  return diagonalizeVelocityFromCoordinates(coordinates);
}

Eigen::Tensor<std::complex<double>, 3>
ElectronH0Wannier::diagonalizeVelocityFromCoordinates(
    Eigen::Vector3d &coordinates) {
  std::vector<Eigen::Vector3d> cartesianWavevectors;
  cartesianWavevectors.push_back(coordinates);
  auto t = batchedDiagonalizeWithVelocities(cartesianWavevectors);
  return std::get<2>(t)[0];
}

FullBandStructure ElectronH0Wannier::populate(Points &fullPoints,
                                              bool &withVelocities,
                                              bool &withEigenvectors,
                                              bool isDistributed) {

  FullBandStructure fullBandStructure(numWannier, particle, withVelocities,
                                      withEigenvectors, fullPoints,
                                      isDistributed);
    std::vector<int> iks = fullBandStructure.getWavevectorIndices();
    int numK = iks.size();

    std::vector<Eigen::Vector3d> cartesianWavevectors(numK);

  #pragma omp parallel for
    for (int iik = 0; iik < numK; iik++) {
      int ik = iks[iik];
      WavevectorIndex ikIdx(ik);
      cartesianWavevectors[iik] = fullBandStructure.getWavevector(ikIdx);
    }

    std::vector<Eigen::VectorXd> allEnergies;
    std::vector<Eigen::MatrixXcd> allEigenvectors;
    std::vector<Eigen::Tensor<std::complex<double>,3>> allVelocities;

    if (withVelocities) {
      auto t = batchedDiagonalizeWithVelocities(cartesianWavevectors);
      allEnergies = std::get<0>(t);
      allEigenvectors = std::get<1>(t);
      allVelocities = std::get<2>(t);
    } else {
      auto t = batchedDiagonalizeFromCoordinates(cartesianWavevectors);
      allEnergies = std::get<0>(t);
      allEigenvectors = std::get<1>(t);
    }

#pragma omp parallel for
    for (int iik = 0; iik < numK; iik++) {
      int ik = iks[iik];
      Point point = fullBandStructure.getPoint(ik);
      auto tup = diagonalize(point);
      auto ens = allEnergies[iik];
      fullBandStructure.setEnergies(point, ens);
      if (withVelocities) {
        auto velocities = diagonalizeVelocity(point);
        fullBandStructure.setVelocities(point, velocities);
      }
      if (withEigenvectors) {
        auto eigenVectors = allEigenvectors[iik];
        fullBandStructure.setEigenvectors(point, eigenVectors);
      }
    }
    return fullBandStructure;
}

std::vector<Eigen::MatrixXcd>
ElectronH0Wannier::getBerryConnection(Point &point) {
  Eigen::Vector3d k = point.getCoordinates(Points::cartesianCoordinates);

  // first we diagonalize the hamiltonian
  auto tup = diagonalize(point);
  // auto ens = std::get<0>(tup);
  auto eigenvectors = std::get<1>(tup);

  // note: the eigenvector matrix is the unitary transformation matrix U
  // from the Bloch to the Wannier gauge.

  std::vector<Eigen::MatrixXcd> berryConnection;

  std::vector<std::complex<double>> phases(bravaisVectors.cols());
  for (int iR = 0; iR < bravaisVectors.cols(); iR++) {
    Eigen::Vector3d R = bravaisVectors.col(iR);
    double phase = k.dot(R);
    std::complex<double> phaseFactor = {cos(phase), sin(phase)};
    phases[iR] = phaseFactor / vectorsDegeneracies(iR);
  }

  for (int i : {0, 1, 2}) {
    // now construct the berryConnection in reciprocal space and Wannier gauge
    Eigen::MatrixXcd berryConnectionW = Eigen::MatrixXcd::Zero(numWannier, numWannier);
    for (int n = 0; n < numWannier; n++) {
      for (int m = 0; m < numWannier; m++) {
          for (int iR = 0; iR < bravaisVectors.cols(); iR++) {
          berryConnectionW(m, n) +=
              phases[iR] * rMatrix(iR, m, n, i);
        }
      }
    }
    Eigen::MatrixXcd thisBerryConnection(numWannier, numWannier);
    thisBerryConnection = eigenvectors.adjoint() * berryConnectionW * eigenvectors;
    berryConnection.push_back(thisBerryConnection);
  }
  return berryConnection;
}

void ElectronH0Wannier::addShiftedVectors(Eigen::Tensor<double,3> degeneracyShifts_,
                       Eigen::Tensor<double,5> vectorsShifts_) {
  // some validation
  if (degeneracyShifts_.dimension(0) != numWannier ||
      degeneracyShifts_.dimension(1) != numWannier ||
      degeneracyShifts_.dimension(2) != numVectors ) {
    Error("Inconsistent dimensions on degeneracyShifts");
  }
  if (vectorsShifts_.dimension(0) != 3 ||
      vectorsShifts_.dimension(1) != 8 ||
      vectorsShifts_.dimension(2) != numWannier ||
      vectorsShifts_.dimension(3) != numWannier ||
      vectorsShifts_.dimension(4) != numVectors) {
    Error("Inconsistent dimensions on vectorsShifts");
  }
  hasShiftedVectors = true;
  vectorsShifts = vectorsShifts_;
  degeneracyShifts = degeneracyShifts_;

  // note: for better performance, I shift all vectors by R,
  // so that in the Fourier transform we only work with one Eigen object
  for (int iR = 0; iR < numVectors; iR++) {
    for (int iw1 = 0; iw1 < numWannier; iw1++) {
      for (int iw2 = 0; iw2 < numWannier; iw2++) {
        for (int iDeg = 0; iDeg < degeneracyShifts(iw1, iw2, iR); ++iDeg) {
          Eigen::Vector3d r = bravaisVectors.col(iR);
          for (int i : {0, 1, 2}) {
            vectorsShifts(i, iDeg, iw1, iw2, iR) += r(i);
          }
        }
      }
    }
  }


  { // copy to GPU
    Kokkos::resize(degeneracyShifts_d, numWannier, numWannier, numVectors);
    Kokkos::resize(vectorsShifts_d, 3, 8, numWannier, numWannier, numVectors);
    auto degeneracyShifts_h = create_mirror_view(degeneracyShifts_d);
    auto vectorsShifts_h = create_mirror_view(vectorsShifts_d);
    for (int iR = 0; iR < numVectors; iR++) {
      for (int i = 0; i < numWannier; i++) {
        for (int j = 0; j < numWannier; j++) {
          degeneracyShifts_h(i, j, iR) = degeneracyShifts(i, j, iR);
          for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 8; ++b) {
              vectorsShifts_h(a, b, i, j, iR) = vectorsShifts(a, b, i, j, iR);
            }
          }
        }
      }
    }
    Kokkos::deep_copy(degeneracyShifts_d, degeneracyShifts_h);
    Kokkos::deep_copy(vectorsShifts_d, vectorsShifts_h);
  }
}


std::tuple<std::vector<Eigen::VectorXd>,
           std::vector<Eigen::MatrixXcd>,
           std::vector<Eigen::Tensor<std::complex<double>, 3>>>
ElectronH0Wannier::batchedDiagonalizeWithVelocities(
    std::vector<Eigen::Vector3d> cartesianCoordinates) {

  int numK = cartesianCoordinates.size();

  double delta = 1.0e-8;
  double threshold = 0.000001 / energyRyToEv; // = 1 micro-eV

  std::vector<Eigen::Vector3d> allVectors(numK * 7);
#pragma omp parallel for
  for (int iK = 0; iK < numK; ++iK) {
    Eigen::Vector3d coordinates = cartesianCoordinates[iK];
    allVectors[iK * 7] = coordinates;
    for (int i : {0, 1, 2}) {
      // define q+ and q- from finite differences.
      Eigen::Vector3d qPlus = coordinates;
      Eigen::Vector3d qMins = coordinates;
      qPlus(i) += delta;
      qMins(i) -= delta;
      allVectors[iK * 7 + i * 2 + 1] = qPlus; // 1, 3, 5, 8
      allVectors[iK * 7 + i * 2 + 2] = qMins;// 2, 4, 6, 9
    }
  }

  auto Hs = batchedBuildHamiltonians(allVectors);

  std::vector<Eigen::VectorXd> resultsEnergies(numK);
  std::vector<Eigen::MatrixXcd> resultsEigenvectors(numK);
  std::vector<Eigen::Tensor<std::complex<double>, 3>> resultsVelocities(numK);

  // find energies and eigenvectors at the desired wavevectors

#pragma omp parallel for
  for (int iK=0; iK<numK; ++iK) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigenSolver(Hs[iK*7]);
    resultsEnergies[iK] = eigenSolver.eigenvalues();
    resultsEigenvectors[iK] = eigenSolver.eigenvectors();
  }

  // now build velocity operator

#pragma omp parallel for
  for (int iK=0; iK<numK; ++iK) {

    Eigen::Tensor<std::complex<double>, 3> velocity(numWannier, numWannier, 3);
    velocity.setZero();

    // if we are working at gamma, we set all velocities to zero.
    if (cartesianCoordinates[iK].norm() < 1.0e-6) {
      resultsVelocities[iK] = velocity;
      continue;
    }

    // now we compute the velocity operator, diagonalizing the expectation
    // value of the derivative of the dynamical matrix.
    // This works better than doing finite differences on the frequencies.
    for (int i : {0, 1, 2}) {
      // build the dynamical matrix at the two wavevectors
      // since we diagonalized it before, A = M.U.M*
      Eigen::MatrixXcd sqrtDPlus = Hs[iK * 7 + i * 2 + 1];
      Eigen::MatrixXcd sqrtDMins = Hs[iK * 7 + i * 2 + 2];

      // now we can build the velocity operator
      Eigen::MatrixXcd der = (sqrtDPlus - sqrtDMins) / (2. * delta);

      // and to be safe, we reimpose hermiticity
      der = 0.5 * (der + der.adjoint());

      // now we rotate in the basis of the eigenvectors at q.
      der = resultsEigenvectors[iK].adjoint() * der * resultsEigenvectors[iK];

      for (int ib1 = 0; ib1 < numWannier; ib1++) {
        for (int ib2 = 0; ib2 < numWannier; ib2++) {
          velocity(ib1, ib2, i) = der(ib1, ib2);
        }
      }
    }

    // turns out that the above algorithm has problems with degenerate bands
    // so, we diagonalize the velocity operator in the degenerate subspace,

    for (int ib = 0; ib < numWannier; ib++) {

      // first, we check if the band is degenerate, and the size of the
      // degenerate subspace
      int sizeSubspace = 1;
      for (int ib2 = ib + 1; ib2 < numWannier; ib2++) {
        // I consider bands degenerate if their frequencies are the same
        // within 0.0001 cm^-1
        if (abs(resultsEnergies[iK](ib) - resultsEnergies[iK](ib2)) > threshold) {
          break;
        }
        sizeSubspace += 1;
      }

      if (sizeSubspace > 1) {
        Eigen::MatrixXcd subMat(sizeSubspace, sizeSubspace);
        // we have to repeat for every direction
        for (int iCart : {0, 1, 2}) {

          // take the velocity matrix of the degenerate subspace
          for (int i = 0; i < sizeSubspace; i++) {
            for (int j = 0; j < sizeSubspace; j++) {
              subMat(i, j) = velocity(ib + i, ib + j, iCart);
            }
          }

          // reinforce hermiticity
          subMat = 0.5 * (subMat + subMat.adjoint());

          // diagonalize the subMatrix
          Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigenSolver(subMat);
          const Eigen::MatrixXcd &newEigenvectors = eigenSolver.eigenvectors();

          // rotate the original matrix in the new basis
          // that diagonalizes the subspace.
          subMat = newEigenvectors.adjoint() * subMat * newEigenvectors;

          // reinforce hermiticity
          subMat = 0.5 * (subMat + subMat.adjoint());

          // substitute back
          for (int i = 0; i < sizeSubspace; i++) {
            for (int j = 0; j < sizeSubspace; j++) {
              velocity(ib + i, ib + j, iCart) = subMat(i, j);
            }
          }
        }
      }

      // we skip the bands in the subspace, since we corrected them already
      ib += sizeSubspace - 1;
    }
    resultsVelocities[iK] = velocity;
  }

  return {resultsEnergies, resultsEigenvectors, resultsVelocities};
}

ComplexView3D ElectronH0Wannier::kokkosBatchedBuildBlochHamiltonian(
    const DoubleView2D &cartesianCoordinates) {

  int numK = cartesianCoordinates.extent(0);

  ComplexView3D hamiltonians("hamiltonians", numK, numWannier, numWannier);

  if (!hasShiftedVectors) {
    ComplexView2D elPhases_d("elPhases_d", numK, numVectors);
    Kokkos::parallel_for(
        "el_phases", Range2D({0, 0}, {numK, numVectors}),
        KOKKOS_LAMBDA(int iK, int iR) {
          double arg = 0.0;
          for (int i = 0; i < 3; i++) {
            arg += cartesianCoordinates(iK, i) * bravaisVectors_d(iR, i);
          }
          elPhases_d(iK, iR) = exp(complexI * arg) / vectorsDegeneracies_d(iR);
        });

    Kokkos::parallel_for(
        "el_hamilton", Range3D({0, 0, 0}, {numK, numWannier, numWannier}),
        KOKKOS_LAMBDA(int iK, int m, int n) {
          Kokkos::complex<double> tmp(0.0);
          for (int iR = 0; iR < numVectors; iR++) {
            tmp += elPhases_d(iK, iR) * h0R_d(m, n, iR);
          }
          hamiltonians(iK, m, n) = tmp;
        });
    Kokkos::realloc(elPhases_d, 0, 0);

  } else {

    Kokkos::parallel_for(
        "elHamiltonianShifted_d",
        Range3D({0, 0, 0}, {numK, numWannier, numWannier}),
        KOKKOS_LAMBDA(int iK, int iw1, int iw2) {
          Kokkos::complex<double> tmp(0.0);
          for (int iR = 0; iR < numVectors; iR++) {
            for (int iDeg = 0; iDeg < degeneracyShifts(iw1, iw2, iR); ++iDeg) {
              double arg = 0.;
              for (int i = 0; i < 3; ++i) {
                arg += cartesianCoordinates(iK, i) * vectorsShifts(i, iDeg, iw1, iw2, iR);
              }
              Kokkos::complex<double> phase = exp(complexI * arg)
                  / vectorsDegeneracies_d(iR) * degeneracyShifts_d(iw1, iw2, iR);
              tmp += h0R_d(iR, iw1, iw2) * phase;
            }
          }
          hamiltonians(iK, iw1, iw2) = tmp;
        });
  }

  return hamiltonians;
}

std::tuple<DoubleView2D, ComplexView3D> ElectronH0Wannier::kokkosBatchedDiagonalizeFromCoordinates(
    const DoubleView2D &cartesianCoordinates) {

  ComplexView3D blochHamiltonians =
      kokkosBatchedBuildBlochHamiltonian(cartesianCoordinates);

  // now the diagonalization.

  int numK = blochHamiltonians.extent(0);
  DoubleView2D allEnergies("energies_d", numK, numWannier);

  kokkosZHEEV(blochHamiltonians, allEnergies);
  // blochHamiltonians now contains eigenvectors

  return {allEnergies, blochHamiltonians};
}

std::tuple<DoubleView2D, ComplexView3D, ComplexView4D>
ElectronH0Wannier::kokkosBatchedDiagonalizeWithVelocities(
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

  // save energies and eigenvectors to results
  DoubleView2D resultEnergies("energies", numK, numWannier);
  ComplexView3D resultEigenvectors("eigenvectors", numK, numWannier, numWannier);
  ComplexView4D resultVelocities("velocities", numK, numWannier, numWannier, 3);

  // put the Hamiltonian matrix in resultEigenvectors
  Kokkos::parallel_for(
      "eigenvectors", Range3D({0, 0, 0}, {numK, numWannier, numWannier}),
      KOKKOS_LAMBDA(int iK, int m, int n) {
        resultEigenvectors(iK, m, n) = allHamiltonians(iK*7, m, n);
      });
  // now, diagonalize the H matrix in place
  kokkosZHEEV(resultEigenvectors, resultEnergies);

  // these are temporary "scratch" memory spaces
  ComplexView3D der("der", numK, numWannier, numWannier);
  ComplexView3D tmpV("tmpV", numK, numWannier, numWannier);

  for (int i = 0; i < 3; ++i) {

    // To build the velocity operator, we first compute the matrix
    // der = dH/dk = ( H(k+dk)-H(k-dk) ) / 2dk
    // we also reinforce the Hermiticity (better numerical properties)
    Kokkos::parallel_for(
        "der", Range3D({0, 0, 0}, {numK, numWannier, numWannier}),
        KOKKOS_LAMBDA(int iK, int m, int n) {
          ComplexView2D HPlus = Kokkos::subview(allHamiltonians, iK * 7 + i * 2 + 1, Kokkos::ALL, Kokkos::ALL);
          ComplexView2D HMins = Kokkos::subview(allHamiltonians, iK * 7 + i * 2 + 2, Kokkos::ALL, Kokkos::ALL);
          der(iK, m, n) = 0.25 / delta * ((HPlus(m, n) - HMins(m, n))
                 + Kokkos::conj(HPlus(n, m) - HMins(n, m)));
        });

    // Now we complete the Hellman-Feynman theorem
    // and compute the velocity as v = U(k)^* der * U(k)
    Kokkos::parallel_for(
        "tmpV", Range3D({0, 0, 0}, {numK, numWannier, numWannier}),
        KOKKOS_LAMBDA(int iK, int m, int n) {
          auto L = Kokkos::subview(resultEigenvectors, iK, Kokkos::ALL, Kokkos::ALL);
          auto R = Kokkos::subview(der, iK, Kokkos::ALL, Kokkos::ALL);
          auto A = Kokkos::subview(tmpV, iK, Kokkos::ALL, Kokkos::ALL);
          Kokkos::complex<double> tmp(0.,0.);
          for (int l = 0; l < numWannier; ++l) {
            tmp += Kokkos::conj(L(l,m)) * R(l, n);
          }
          A(m, n) = tmp;
        });

    Kokkos::parallel_for(
        "vel", Range3D({0, 0, 0}, {numK, numWannier, numWannier}),
        KOKKOS_LAMBDA(int iK, int m, int n) {
          double norm = 0.;
          for (int i=0; i<3; ++i) {
            norm += cartesianCoordinates(iK,i)*cartesianCoordinates(iK,i);
          }
          Kokkos::complex<double> tmp(0.,0.);
          if ( norm > 1.0e-6 ) {// skip the gamma point
            for (int l = 0; l < numWannier; ++l) {
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
  IntView2D degCount("tmpV", numK, numWannier);
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
        for (int ib = 0; ib < numWannier; ib++) {
          // first, we check if the band is degenerate, and the size of the
          // degenerate subspace
          int sizeSubspace = 1;
          for (int ib2 = ib + 1; ib2 < numWannier; ib2++) {
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
          for (int ib = 0; ib < numWannier; ib++) {
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

    IntView1D scan("scan", numK*numWannier);
    Kokkos::parallel_for(
        Range2D({0,0},{numK, numWannier}), KOKKOS_LAMBDA(int ik, int ib) {
          if ( iDeg == degCount(ik,ib) ) {
            scan(ik*numWannier+ib) = 1;
          } else {
            scan(ik*numWannier+ib) = 0;
          }
        });

    Kokkos::parallel_scan("scan", numK*numWannier, KOKKOS_LAMBDA(const int i,
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
        "funcKB", Range2D({0,0},{numK, numWannier}), KOKKOS_LAMBDA(int ik, int ib) {
          int ikb = ik*numWannier + ib;
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

