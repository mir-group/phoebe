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
    Error e("WannierH0(): h0R should have dimensions (R,bands,bands)");
  }
  if (h0R.dimension(0) != bravaisVectors.cols()) {
    Error e("WannierH0(): h0R and bravaisVectors not aligned");
  }
  if (vectorsDegeneracies.size() != bravaisVectors.cols()) {
    Error e("WannierH0(): degeneracies not aligned with vectors");
  }

  if ((rMatrix.dimension(1) != h0R.dimension(0)) ||
      (rMatrix.dimension(2) != h0R.dimension(1)) ||
      (rMatrix.dimension(3) != h0R.dimension(2))) {
    Error e("WannierH0(): h0R and rMatrix should be aligned");
  }

  if (rMatrix.dimension(0) != 3) {
    Error e("WannierH0(): rMatrix should be a vector");
  }

  numBands = h0R.dimension(1);
  numVectors = vectorsDegeneracies.size();
}

// copy constructor
ElectronH0Wannier::ElectronH0Wannier(const ElectronH0Wannier &that)
    : particle(Particle::electron) {
  h0R = that.h0R;
  rMatrix = that.rMatrix;
  directUnitCell = that.directUnitCell;
  numBands = that.numBands;
  bravaisVectors = that.bravaisVectors;
  numVectors = that.numVectors;
  vectorsDegeneracies = that.vectorsDegeneracies;
}

// copy assignment
ElectronH0Wannier &ElectronH0Wannier::operator=(const ElectronH0Wannier &that) {
  if (this != &that) {
    bravaisVectors.resize(0, 0);
    vectorsDegeneracies.resize(0);
    h0R.resize(0, 0, 0);
    rMatrix.resize(0, 0, 0, 0);
    particle = that.particle;
    numVectors = that.numVectors;
    numBands = that.numBands;
    bravaisVectors = that.bravaisVectors;
    vectorsDegeneracies = that.vectorsDegeneracies;
    directUnitCell = that.directUnitCell;
    h0R = that.h0R;
    rMatrix = that.rMatrix;
  }
  return *this;
}

Particle ElectronH0Wannier::getParticle() { return particle; }

long ElectronH0Wannier::getNumBands() { return numBands; }

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
ElectronH0Wannier::diagonalize(Point &point) {
  Eigen::Vector3d k = point.getCoords(Points::cartesianCoords);

  auto tup = diagonalizeFromCoords(k);
  auto energies = std::get<0>(tup);
  auto eigenvectors = std::get<1>(tup);

  // note: the eigenvector matrix is the unitary transformation matrix U
  // from the Bloch to the Wannier gauge.

  return {energies, eigenvectors};
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
ElectronH0Wannier::diagonalizeFromCoords(Eigen::Vector3d &k) {

  Eigen::MatrixXcd h0K(numBands, numBands);
  h0K.setZero();

  for (long iR = 0; iR < bravaisVectors.cols(); iR++) {
    Eigen::Vector3d R = bravaisVectors.col(iR);
    double phase = k.dot(R);
    std::complex<double> phaseFactor = {cos(phase), sin(phase)};
    for (long m = 0; m < numBands; m++) {
      for (long n = 0; n < numBands; n++) {
        h0K(m, n) += phaseFactor * h0R(iR, m, n) / vectorsDegeneracies(iR);
      }
    }
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(h0K);
  Eigen::VectorXd energies = eigensolver.eigenvalues();
  Eigen::MatrixXcd eigenvectors = eigensolver.eigenvectors();

  return {energies, eigenvectors};
}

Eigen::Tensor<std::complex<double>, 3>
ElectronH0Wannier::diagonalizeVelocity(Point &point) {
  Eigen::Vector3d coords = point.getCoords(Points::cartesianCoords);
  return diagonalizeVelocityFromCoords(coords);
}

Eigen::Tensor<std::complex<double>, 3>
ElectronH0Wannier::diagonalizeVelocityFromCoords(Eigen::Vector3d &coords) {
  double delta = 1.0e-8;
  double threshold = 0.000001 / energyRyToEv; // = 1 micro-eV

  Eigen::Tensor<std::complex<double>, 3> velocity(numBands, numBands, 3);
  velocity.setZero();

  // if we are working at gamma, we set all velocities to zero.
  if (coords.norm() < 1.0e-6) {
    return velocity;
  }

  // get the eigenvectors and the energies of the q-point
  auto tup = diagonalizeFromCoords(coords);
  auto energies = std::get<0>(tup);
  auto eigenvectors = std::get<1>(tup);

  // now we compute the velocity operator, diagonalizing the expectation
  // value of the derivative of the dynamical matrix.
  // This works better than doing finite differences on the frequencies.
  for (long i = 0; i < 3; i++) {
    // define q+ and q- from finite differences.
    Eigen::Vector3d qPlus = coords;
    Eigen::Vector3d qMins = coords;
    qPlus(i) += delta;
    qMins(i) -= delta;

    // diagonalize the dynamical matrix at q+ and q-
    auto tup = diagonalizeFromCoords(qPlus);
    auto enPlus = std::get<0>(tup);
    auto eigPlus = std::get<1>(tup);
    auto tup1 = diagonalizeFromCoords(qMins);
    auto enMins = std::get<0>(tup1);
    auto eigMins = std::get<1>(tup1);

    // build diagonal matrices with frequencies
    Eigen::MatrixXd enPlusMat(numBands, numBands);
    Eigen::MatrixXd enMinsMat(numBands, numBands);
    enPlusMat.setZero();
    enMinsMat.setZero();
    enPlusMat.diagonal() << enPlus;
    enMinsMat.diagonal() << enMins;

    // build the dynamical matrix at the two wavevectors
    // since we diagonalized it before, A = M.U.M*
    Eigen::MatrixXcd sqrtDPlus(numBands, numBands);
    sqrtDPlus = eigPlus * enPlusMat * eigPlus.adjoint();
    Eigen::MatrixXcd sqrtDMins(numBands, numBands);
    sqrtDMins = eigMins * enMinsMat * eigMins.adjoint();

    // now we can build the velocity operator
    Eigen::MatrixXcd der(numBands, numBands);
    der = (sqrtDPlus - sqrtDMins) / (2. * delta);

    // and to be safe, we reimpose hermiticity
    der = 0.5 * (der + der.adjoint());

    // now we rotate in the basis of the eigenvectors at q.
    der = eigenvectors.adjoint() * der * eigenvectors;

    for (long ib1 = 0; ib1 < numBands; ib1++) {
      for (long ib2 = 0; ib2 < numBands; ib2++) {
        velocity(ib1, ib2, i) = der(ib1, ib2);
      }
    }
  }

  // turns out that the above algorithm has problems with degenerate bands
  // so, we diagonalize the velocity operator in the degenerate subspace,

  for (long ib = 0; ib < numBands; ib++) {

    // first, we check if the band is degenerate, and the size of the
    // degenerate subspace
    long sizeSubspace = 1;
    for (long ib2 = ib + 1; ib2 < numBands; ib2++) {
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
      for (long iCart = 0; iCart < 3; iCart++) {

        // take the velocity matrix of the degenerate subspace
        for (long i = 0; i < sizeSubspace; i++) {
          for (long j = 0; j < sizeSubspace; j++) {
            subMat(i, j) = velocity(ib + i, ib + j, iCart);
          }
        }

        // reinforce hermiticity
        subMat = 0.5 * (subMat + subMat.adjoint());

        // diagonalize the subMatrix
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(subMat);
        Eigen::MatrixXcd newEigvecs = eigensolver.eigenvectors();

        // rotate the original matrix in the new basis
        // that diagonalizes the subspace.
        subMat = newEigvecs.adjoint() * subMat * newEigvecs;

        // reinforce hermiticity
        subMat = 0.5 * (subMat + subMat.adjoint());

        // substitute back
        for (long i = 0; i < sizeSubspace; i++) {
          for (long j = 0; j < sizeSubspace; j++) {
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
                                              bool &withEigenvectors) {

  FullBandStructure fullBandStructure(numBands, particle, withVelocities,
                                      withEigenvectors, fullPoints);

  for (long ik = 0; ik < fullBandStructure.getNumPoints(); ik++) {
    Point point = fullBandStructure.getPoint(ik);
    auto tup = diagonalize(point);
    auto ens = std::get<0>(tup);
    auto eigvecs = std::get<1>(tup);
    fullBandStructure.setEnergies(point, ens);
    if (withVelocities) {
      auto vels = diagonalizeVelocity(point);
      fullBandStructure.setVelocities(point, vels);
    }
    if (withEigenvectors) {
      fullBandStructure.setEigenvectors(point, eigvecs);
    }
  }
  return fullBandStructure;
}

std::vector<Eigen::MatrixXcd>
ElectronH0Wannier::getBerryConnection(Point &point) {
  Eigen::Vector3d k = point.getCoords(Points::cartesianCoords);

  // first we diagonalize the hamiltonian
  auto tup = diagonalize(point);
  auto ens = std::get<0>(tup);
  auto eigvecs = std::get<1>(tup);

  // note: the eigenvector matrix is the unitary transformation matrix U
  // from the Bloch to the Wannier gauge.

  std::vector<Eigen::MatrixXcd> bc;

  for (long i = 0; i < 3; i++) {

    // now construct the berryConnection in reciprocal space and Wannier gauge
    Eigen::MatrixXcd berryConnectionW(numBands, numBands);
    berryConnectionW.setZero();

    for (long iR = 0; iR < bravaisVectors.cols(); iR++) {
      Eigen::Vector3d R = bravaisVectors.col(iR);
      double phase = k.dot(R);
      std::complex<double> phaseFactor = {cos(phase), sin(phase)};
      for (long m = 0; m < numBands; m++) {
        for (long n = 0; n < numBands; n++) {
          berryConnectionW(m, n) +=
              phaseFactor * rMatrix(i, iR, m, n) / vectorsDegeneracies(iR);
        }
      }
    }

    Eigen::MatrixXcd berryConnection(numBands, numBands);
    berryConnection = eigvecs.adjoint() * berryConnectionW * eigvecs;
    bc.push_back(berryConnection);
  }
  return bc;
}
