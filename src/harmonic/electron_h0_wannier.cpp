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

  if ((rMatrix.dimension(1) != h0R.dimension(0)) || (rMatrix.dimension(2) != h0R.dimension(1)) || (rMatrix.dimension(3) != h0R.dimension(2))) {
    Error("WannierH0(): h0R and rMatrix should be aligned");
  }

  if (rMatrix.dimension(0) != 3) {
    Error("WannierH0(): rMatrix should be a vector");
  }

  numWannier = int(h0R.dimension(1));
  numVectors = int(vectorsDegeneracies.size());

  // copy data to the GPU
  {
    HostComplexView3D h0R_h((Kokkos::complex<double> *) h0R.data(),
                            numVectors, numWannier, numWannier);
    HostDoubleView1D vectorsDegeneracies_h(vectorsDegeneracies.data(), numVectors);
    HostDoubleView2D bravaisVectors_h(bravaisVectors.data(), 3, numVectors);
    for (int iR = 0; iR < numVectors; iR++) {
      for (int i = 0; i < numWannier; i++) {
        for (int j = 0; j < numWannier; j++) {
          h0R_h(iR, i, j) = h0R(iR, i, j);
        }
      }
      vectorsDegeneracies_h(iR) = vectorsDegeneracies(iR);
      for (int i = 0; i < 3; ++i) {
        bravaisVectors_h(i, iR) = bravaisVectors(i, iR);
      }
    }
    h0R_d = Kokkos::create_mirror_view(Kokkos::DefaultExecutionSpace(), h0R_h);
    vectorsDegeneracies_d = Kokkos::create_mirror_view(Kokkos::DefaultExecutionSpace(), vectorsDegeneracies_h);
    bravaisVectors_d = Kokkos::create_mirror_view(Kokkos::DefaultExecutionSpace(), bravaisVectors_h);
    Kokkos::deep_copy(h0R_d, h0R_h);
    Kokkos::deep_copy(vectorsDegeneracies_d, vectorsDegeneracies_h);
    Kokkos::deep_copy(bravaisVectors_d, bravaisVectors_h);

    if (hasShiftedVectors) {
      HostDoubleView3D degeneracyShifts_h(degeneracyShifts.data(),
                                          numWannier, numWannier, numVectors);
      HostDoubleView5D vectorsShifts_h(vectorsShifts.data(),
                                       3, 8, numWannier, numWannier, numVectors);
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
      degeneracyShifts_d = Kokkos::create_mirror_view(Kokkos::DefaultExecutionSpace(), degeneracyShifts_h);
      vectorsShifts_d = Kokkos::create_mirror_view(Kokkos::DefaultExecutionSpace(), vectorsShifts_h);
      Kokkos::deep_copy(degeneracyShifts_d, degeneracyShifts_h);
      Kokkos::deep_copy(vectorsShifts_d, vectorsShifts_h);
    }
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
  // internalPopulate acts on a std::vector of wavevectors
  std::vector<Eigen::Vector3d> kVecs;
  kVecs.push_back(k);
  auto t = internalPopulate(kVecs);
  std::vector<Eigen::VectorXd> eigvals = std::get<0>(t);
  std::vector<Eigen::MatrixXcd> eigvecs = std::get<1>(t);
  return {eigvals[0], eigvecs[0]};
}

Eigen::Tensor<std::complex<double>, 3>
ElectronH0Wannier::diagonalizeVelocity(Point &point) {
  Eigen::Vector3d coordinates = point.getCoordinates(Points::cartesianCoordinates);
  return diagonalizeVelocityFromCoordinates(coordinates);
}

Eigen::Tensor<std::complex<double>, 3>
ElectronH0Wannier::diagonalizeVelocityFromCoordinates(
    Eigen::Vector3d &coordinates) {
  std::vector<Eigen::Vector3d> wavevectors;
  wavevectors.push_back(coordinates);
  auto t = populate(wavevectors, true);
  std::vector<Eigen::Tensor<std::complex<double>, 3>> velocities = std::get<2>(t);
  return velocities[0];
}

FullBandStructure ElectronH0Wannier::populate(Points &fullPoints,
                                              bool &withVelocities,
                                              bool &withEigenvectors,
                                              bool isDistributed) {

  FullBandStructure fullBandStructure(numWannier, particle, withVelocities,
                                      withEigenvectors, fullPoints,
                                      isDistributed);

  std::vector<int> iks = fullBandStructure.getWavevectorIndices();
  int niks = iks.size();

  // first prepare the list of wavevectors
  std::vector<Eigen::Vector3d> allWavevectors(niks);
#pragma omp parallel for
  for (int iik = 0; iik < niks; iik++) {
    int ik = iks[iik];
    Eigen::Vector3d k = fullPoints.getPointCoordinates(ik, Points::cartesianCoordinates);
    allWavevectors[iik] = k;
  }

  // then call the function to diagonalize these wavevectors
  auto t = populate(allWavevectors, withVelocities);

  // retrieve results
  std::vector<Eigen::VectorXd> allEnergies = std::get<0>(t);
  std::vector<Eigen::MatrixXcd> allEigenvectors = std::get<1>(t);
  std::vector<Eigen::Tensor<std::complex<double>, 3>> allVelocities = std::get<2>(t);

  for (int iik = 0; iik < niks; iik++) {
    int ik = iks[iik];
    auto point = fullBandStructure.getPoint(ik);
    fullBandStructure.setEnergies(point, allEnergies[iik]);
    if (withVelocities) {
      fullBandStructure.setVelocities(point, allVelocities[iik]);
    }
    if (withEigenvectors) {
      fullBandStructure.setEigenvectors(point, allEigenvectors[iik]);
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

void ElectronH0Wannier::addShiftedVectors(Eigen::Tensor<double, 3> degeneracyShifts_,
                                          Eigen::Tensor<double, 5> vectorsShifts_) {
  // some validation
  if (degeneracyShifts_.dimension(0) != numWannier || degeneracyShifts_.dimension(1) != numWannier || degeneracyShifts_.dimension(2) != numVectors) {
    Error("Inconsistent dimensions on degeneracyShifts");
  }
  if (vectorsShifts_.dimension(0) != 3 || vectorsShifts_.dimension(1) != 8 || vectorsShifts_.dimension(2) != numWannier || vectorsShifts_.dimension(3) != numWannier || vectorsShifts_.dimension(4) != numVectors) {
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
}

std::tuple<std::vector<Eigen::VectorXd>,
           std::vector<Eigen::MatrixXcd>>
ElectronH0Wannier::internalPopulate(
    const std::vector<Eigen::Vector3d> &cartesianCoordinates) {

  // note: this is like a copy of the above function populate()
  // we should rely only on this one, as it is more general.

  int numK = cartesianCoordinates.size();
  std::vector<Eigen::MatrixXcd> blochHamiltonians(numK, Eigen::MatrixXcd::Zero(numWannier, numWannier));

  if (!hasShiftedVectors) {
    Eigen::MatrixXcd phases(numVectors, numK);
    for (int iK = 0; iK < numK; ++iK) {
      Eigen::Vector3d k = cartesianCoordinates[iK];
      for (int iR = 0; iR < numVectors; iR++) {
        double phase = k.dot(bravaisVectors.col(iR));
        std::complex<double> phaseFactor = {cos(phase), sin(phase)};
        phases(iR, iK) = phaseFactor / vectorsDegeneracies(iR);
      }
    }
    for (int iK = 0; iK < numK; ++iK) {
      for (int n = 0; n < numWannier; n++) {
        for (int m = 0; m < numWannier; m++) {
          for (int iR = 0; iR < numVectors; iR++) {
            blochHamiltonians[iK](m, n) += phases(iR, iK) * h0R(iR, m, n);
          }
        }
      }
    }

  } else {

    for (int iK = 0; iK < numK; ++iK) {
      Eigen::Vector3d k = cartesianCoordinates[iK];
      for (int iR = 0; iR < numVectors; iR++) {
        double phaseArg;
        std::complex<double> phase;
        for (int iw2 = 0; iw2 < numWannier; iw2++) {
          for (int iw1 = 0; iw1 < numWannier; iw1++) {
            for (int iDeg = 0; iDeg < degeneracyShifts(iw1, iw2, iR); ++iDeg) {
              phaseArg = 0.;
              for (int i : {0, 1, 2}) {
                phaseArg += k(i) * vectorsShifts(i, iDeg, iw1, iw2, iR);
              }
              phase = {cos(phaseArg), sin(phaseArg)};
              phase /= vectorsDegeneracies(iR) * degeneracyShifts(iw1, iw2, iR);
              blochHamiltonians[iK](iw1, iw2) += phase * h0R(iR, iw1, iw2);
            }
          }
        }
      }
    }
  }

  std::vector<Eigen::VectorXd> allEnergies(numK);
  std::vector<Eigen::MatrixXcd> allEigenvectors(numK);

  for (int iK = 0; iK < numK; ++iK) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigenSolver(blochHamiltonians[iK]);
    allEnergies[iK] = eigenSolver.eigenvalues();
    allEigenvectors[iK] = eigenSolver.eigenvectors();
  }
  return {allEnergies, allEigenvectors};
}

ComplexView3D ElectronH0Wannier::kokkosBuildBlochHamiltonian(
    const DoubleView2D &cartesianCoordinates) {

  int numK = cartesianCoordinates.extent(0);

  ComplexView3D blochHamiltonians_d("hamiltonians", numK, numWannier, numWannier);

  if (!hasShiftedVectors) {

    ComplexView2D elPhases_d("elPhases", numK, numVectors);
    Kokkos::parallel_for(
        "elPhases_d", Range2D({0, 0}, {numK, numVectors}),
        KOKKOS_LAMBDA(int iK, int iR) {
          double arg = 0.0;
          for (int i = 0; i < 3; i++) {
            arg += cartesianCoordinates(iK, i) * bravaisVectors_d(i, iR);
          }
          elPhases_d(iK, iR) = exp(complexI * arg) / vectorsDegeneracies_d(iR);
        });

    Kokkos::parallel_for(
        "elHamiltonian_d",
        Range3D({0, 0, 0}, {numK, numWannier, numWannier}),
        KOKKOS_LAMBDA(int iK, int n, int m) {
          Kokkos::complex<double> tmp(0.0);
          for (int iR = 0; iR < numVectors; iR++) {
            tmp += elPhases_d(iK, iR) * h0R_d(iR, m, n);
          }
          blochHamiltonians_d(iK, m, n) = tmp;
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
          blochHamiltonians_d(iK, iw1, iw2) = tmp;
        });
  }

  return blochHamiltonians_d;
}

std::tuple<DoubleView2D, ComplexView3D> ElectronH0Wannier::kokkosInternalPopulate(
    const DoubleView2D &cartesianCoordinates) {

  ComplexView3D blochHamiltonians = kokkosBuildBlochHamiltonian(cartesianCoordinates);

  // now the diagonalization.

  int numK = blochHamiltonians.extent(0);
  DoubleView2D allEnergies("energies_d", numK, numWannier);

  kokkosZHEEV(blochHamiltonians, allEnergies);
  // blochHamiltonians now contains eigenvectors

  return {allEnergies, blochHamiltonians};
}

std::tuple<std::vector<Eigen::VectorXd>, std::vector<Eigen::MatrixXcd>,
           std::vector<Eigen::Tensor<std::complex<double>, 3>>>
ElectronH0Wannier::populate(
    const std::vector<Eigen::Vector3d> &cartesianCoordinates,
    const bool &withVelocities) {

  if (!withVelocities) {
    auto t = internalPopulate(cartesianCoordinates);
    std::vector<Eigen::VectorXd> allEnergies = std::get<0>(t);
    std::vector<Eigen::MatrixXcd> allEigenvectors = std::get<1>(t);
    std::vector<Eigen::Tensor<std::complex<double>, 3>> allVelocities;
    return {allEnergies, allEigenvectors, allVelocities};

  } else {

    int numK = cartesianCoordinates.size();

    double delta = 1.0e-8;
    double threshold = 0.000001 / energyRyToEv;// = 1 micro-eV

    std::vector<Eigen::Vector3d> allVectors(numK * 7);
#pragma omp parallel for
    for (int iK = 0; iK < numK; ++iK) {
      Eigen::Vector3d coordinates = cartesianCoordinates[iK];
      allVectors[iK * 7] = coordinates;
      for (int i : {0, 1, 2}) {
        // define q+ and q- from finite differences.
        Eigen::Vector3d qPlus = coordinates;
        Eigen::Vector3d qMinus = coordinates;
        qPlus(i) += delta;
        qMinus(i) -= delta;
        allVectors[iK * 7 + i * 2 + 1] = qPlus; // 1, 3, 5, 8
        allVectors[iK * 7 + i * 2 + 2] = qMinus;// 2, 4, 6, 9
      }
    }

    auto t = internalPopulate(allVectors);
    auto allEnergies = std::get<0>(t);
    auto allEigenvectors = std::get<1>(t);

    std::vector<Eigen::VectorXd> resultEnergies(numK);
    std::vector<Eigen::MatrixXcd> resultEigenvectors(numK);
    std::vector<Eigen::Tensor<std::complex<double>, 3>> resultVelocities(numK);

    Eigen::Tensor<std::complex<double>, 3> velocity(numWannier, numWannier, 3);
    for (int iK = 0; iK < numK; ++iK) {
      Eigen::Vector3d k = allVectors[iK * 7];
      Eigen::VectorXd energies = allEnergies[iK * 7];
      Eigen::MatrixXcd eigenvectors = allEigenvectors[iK * 7];

      resultEnergies[iK] = energies;
      resultEigenvectors[iK] = eigenvectors;
      velocity.setZero();

      if (k.norm() < 1.0e-6) {// if not zero, velocity is non-analytical
        resultVelocities[iK] = velocity;
        continue;
      }

      for (int i : {0, 1, 2}) {
        Eigen::VectorXd enPlus = allEnergies[iK * 7 + i * 2 + 1];
        Eigen::VectorXd enMinus = allEnergies[iK * 7 + i * 2 + 2];
        Eigen::MatrixXcd eigPlus = allEigenvectors[iK * 7 + i * 2 + 1];
        Eigen::MatrixXcd eigMinus = allEigenvectors[iK * 7 + i * 2 + 2];

        // build diagonal matrices with frequencies
        Eigen::MatrixXd enPlusMat(numWannier, numWannier);
        Eigen::MatrixXd enMinusMat(numWannier, numWannier);
        enPlusMat.setZero();
        enMinusMat.setZero();
        enPlusMat.diagonal() << enPlus;
        enMinusMat.diagonal() << enMinus;

        // build the dynamical matrix at the two wavevectors
        // since we diagonalized it before, A = M.U.M*
        // Eigen::MatrixXcd sqrtDPlus(numWannier, numWannier);
        Eigen::MatrixXcd sqrtDPlus = eigPlus * enPlusMat * eigPlus.adjoint();
        // Eigen::MatrixXcd sqrtDMinus(numWannier, numWannier);
        Eigen::MatrixXcd sqrtDMinus = eigMinus * enMinusMat * eigMinus.adjoint();

        // now we can build the velocity operator
        // Eigen::MatrixXcd der(numWannier, numWannier);
        Eigen::MatrixXcd der = (sqrtDPlus - sqrtDMinus) / (2. * delta);

        // and to be safe, we reimpose hermiticity
        der = 0.5 * (der + der.adjoint());

        // now we rotate in the basis of the eigenvectors at q.
        der = eigenvectors.adjoint() * der * eigenvectors;

        for (int ib1 = 0; ib1 < numWannier; ib1++) {
          for (int ib2 = 0; ib2 < numWannier; ib2++) {
            velocity(ib1, ib2, i) = der(ib1, ib2);
          }
        }
      }

      {
        // turns out that the above algorithm has problems with degenerate bands
        // so, we diagonalize the velocity operator in the degenerate subspace,
        for (int ib = 0; ib < numWannier; ib++) {

          // first, we check if the band is degenerate, and the size of the
          // degenerate subspace
          int sizeSubspace = 1;
          for (int ib2 = ib + 1; ib2 < numWannier; ib2++) {
            // I consider bands degenerate if their frequencies are the same
            // within 0.0001 cm^-1
            if (abs(energies(ib) - energies(ib2)) > threshold) {
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
      }
      resultVelocities[iK] = velocity;
    }
    return {resultEnergies, resultEigenvectors, resultVelocities};
  }
}

std::tuple<DoubleView2D, ComplexView3D, ComplexView4D>
ElectronH0Wannier::kokkosPopulate(const DoubleView2D &cartesianCoordinates,
                                  const bool &withVelocities) {

  if (!withVelocities) {
    auto t = kokkosInternalPopulate(cartesianCoordinates);
    DoubleView2D allEnergies = std::get<0>(t);
    ComplexView3D allEigenvectors = std::get<1>(t);
    ComplexView4D allVelocities;
    return {allEnergies, allEigenvectors, allVelocities};

  } else {

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
    ComplexView3D allHamiltonians = kokkosBuildBlochHamiltonian(allVectors);

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
      ComplexView3D der("der", numK, numWannier, numWannier);
      Kokkos::parallel_for(
          "der", Range3D({0, 0, 0}, {numK, numWannier, numWannier}),
          KOKKOS_LAMBDA(int iK, int m, int n) {
            der(iK, m, n) = 0.25 / delta
                * ((allHamiltonians(iK * 7 + i * 2 + 1, m, n)
                    - allHamiltonians(iK * 7 + i * 2 + 2, m, n))
                   + Kokkos::conj(allHamiltonians(iK * 7 + i * 2 + 1, n, m)
                                  - allHamiltonians(iK * 7 + i * 2 + 2, n, m)));
          });

      // Now we complete the Hellman-Feynman theorem
      // and compute the velocity as v = U(k)^* der * U(k)
      Kokkos::parallel_for(
          "tmpV", Range3D({0, 0, 0}, {numK, numWannier, numWannier}),
          KOKKOS_LAMBDA(int iK, int m, int n) {
            Kokkos::complex<double> tmp = zero;
            for (int l = 0; l < numWannier; ++l) {
              tmp += Kokkos::conj(resultEigenvectors(iK, l, m)) * der(iK, l, n);
            }
            tmpV(iK, m, n) = tmp;
          });

      Kokkos::parallel_for(
          "vel", Range3D({0, 0, 0}, {numK, numWannier, numWannier}),
          KOKKOS_LAMBDA(int iK, int m, int n) {
            Kokkos::complex<double> tmp = zero;
            for (int l = 0; l < numWannier; ++l) {
              tmp += tmpV(iK, m, l) * resultEigenvectors(iK, l, n);
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
              tmp += tmpSubMat(iMatCart, i, k) * newEigenvectors(iMatCart, k, i);
            }
            resultVelocities(iK, ib + i, ib + j, iCart) = tmp;
          });
    }

    return {resultEnergies, resultEigenvectors, resultVelocities};
  }
}
