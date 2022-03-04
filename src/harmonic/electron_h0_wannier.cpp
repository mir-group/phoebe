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
  double delta = 1.0e-8;
  double threshold = 0.000001 / energyRyToEv; // = 1 micro-eV

  Eigen::Tensor<std::complex<double>, 3> velocity(numWannier, numWannier, 3);
  velocity.setZero();

  // if we are working at gamma, we set all velocities to zero.
  if (coordinates.norm() < 1.0e-6) {
    return velocity;
  }

  // get the eigenvectors and the energies of the q-point
  auto tup = diagonalizeFromCoordinates(coordinates);
  auto energies = std::get<0>(tup);
  auto eigenvectors = std::get<1>(tup);

  // now we compute the velocity operator, diagonalizing the expectation
  // value of the derivative of the dynamical matrix.
  // This works better than doing finite differences on the frequencies.
  for (int i : {0,1,2}) {
    // define q+ and q- from finite differences.
    Eigen::Vector3d qPlus = coordinates;
    Eigen::Vector3d qMinus = coordinates;
    qPlus(i) += delta;
    qMinus(i) -= delta;

    // diagonalize the dynamical matrix at q+ and q-
    auto tup2 = diagonalizeFromCoordinates(qPlus);
    auto enPlus = std::get<0>(tup2);
    auto eigPlus = std::get<1>(tup2);
    auto tup1 = diagonalizeFromCoordinates(qMinus);
    auto enMinus = std::get<0>(tup1);
    auto eigMinus = std::get<1>(tup1);

    // build diagonal matrices with frequencies
    Eigen::MatrixXd enPlusMat(numWannier, numWannier);
    Eigen::MatrixXd enMinusMat(numWannier, numWannier);
    enPlusMat.setZero();
    enMinusMat.setZero();
    enPlusMat.diagonal() << enPlus;
    enMinusMat.diagonal() << enMinus;

    // build the dynamical matrix at the two wavevectors
    // since we diagonalized it before, A = M.U.M*
    Eigen::MatrixXcd sqrtDPlus(numWannier, numWannier);
    sqrtDPlus = eigPlus * enPlusMat * eigPlus.adjoint();
    Eigen::MatrixXcd sqrtDMinus(numWannier, numWannier);
    sqrtDMinus = eigMinus * enMinusMat * eigMinus.adjoint();

    // now we can build the velocity operator
    Eigen::MatrixXcd der(numWannier, numWannier);
    der = (sqrtDPlus - sqrtDMinus) / (2. * delta);

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
      for (int iCart : {0,1,2}) {

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
        const Eigen::MatrixXcd&newEigenvectors = eigenSolver.eigenvectors();

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
  return velocity;
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


//  std::vector<int> iks = fullBandStructure.getWavevectorIndices();
//  int numK = iks.size();
//#pragma omp parallel for
//  for (int iik = 0; iik < numK; iik++){
//    int ik = iks[iik];
//    Point point = fullBandStructure.getPoint(ik);
//    auto tup = diagonalize(point);
//    auto ens = std::get<0>(tup);
//    auto eigenVectors = std::get<1>(tup);
//    fullBandStructure.setEnergies(point, ens);
//    if (withVelocities) {
//      auto velocities = diagonalizeVelocity(point);
//      fullBandStructure.setVelocities(point, velocities);
//    }
//    if (withEigenvectors) {
//      fullBandStructure.setEigenvectors(point, eigenVectors);
//    }
//  }
//  return fullBandStructure;
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
      Eigen::Vector3d qMinus = coordinates;
      qPlus(i) += delta;
      qMinus(i) -= delta;
      allVectors[iK * 7 + i * 2 + 1] = qPlus; // 1, 3, 5, 8
      allVectors[iK * 7 + i * 2 + 2] = qMinus;// 2, 4, 6, 9
    }
  }

  auto Hs = batchedBuildHamiltonians(allVectors);

  std::vector<Eigen::VectorXd> resultsEnergies(numK);
  std::vector<Eigen::MatrixXcd> resultsEigenvectors(numK);
  std::vector<Eigen::Tensor<std::complex<double>, 3>> resultsVelocities(numK);

  // find energies and eigenvectors at the desired wavevectors

  for (int iK=0; iK<numK; ++iK) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigenSolver(Hs[iK*7]);
    resultsEnergies[iK] = eigenSolver.eigenvalues();
    resultsEigenvectors[iK] = eigenSolver.eigenvectors();
  }

  // now build velocity operator

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
      Eigen::MatrixXcd sqrtDMinus = Hs[iK * 7 + i * 2 + 2];

      // now we can build the velocity operator
      Eigen::MatrixXcd der = (sqrtDPlus - sqrtDMinus) / (2. * delta);

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
