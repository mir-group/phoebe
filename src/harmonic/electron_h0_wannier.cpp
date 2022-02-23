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
  }
  return *this;
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
  std::vector<Eigen::Tensor<std::complex<double>,3>> velocities = std::get<2>(t);
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
  for(int iik = 0; iik < niks; iik++) {
    int ik = iks[iik];
    Eigen::Vector3d k = fullPoints.getPointCoordinates(ik, Points::cartesianCoordinates);
    allWavevectors[iik] = k;
  }

  // then call the function to diagonalize these wavevectors
  auto t = populate(allWavevectors, withVelocities);

  // retrieve results
  std::vector<Eigen::VectorXd> allEnergies = std::get<0>(t);
  std::vector<Eigen::MatrixXcd> allEigenvectors = std::get<1>(t);
  std::vector<Eigen::Tensor<std::complex<double>,3>> allVelocities = std::get<2>(t);

  for(int iik = 0; iik < niks; iik++) {
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
           std::vector<Eigen::MatrixXcd>> ElectronH0Wannier::internalPopulate(
    const std::vector<Eigen::Vector3d>& cartesianCoordinates) {

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

  for (int iK=0; iK<numK; ++iK) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigenSolver(blochHamiltonians[iK]);
    allEnergies[iK] = eigenSolver.eigenvalues();
    allEigenvectors[iK] = eigenSolver.eigenvectors();
  }
  return {allEnergies, allEigenvectors};
}

std::tuple<std::vector<Eigen::VectorXd>, std::vector<Eigen::MatrixXcd>,
           std::vector<Eigen::Tensor<std::complex<double>,3>>
           > ElectronH0Wannier::populate(
    const std::vector<Eigen::Vector3d>& cartesianCoordinates,
    const bool& withVelocities) {

  int numK = cartesianCoordinates.size();

  if (!withVelocities) {
    auto t = internalPopulate(cartesianCoordinates);
    std::vector<Eigen::VectorXd> allEnergies = std::get<0>(t);
    std::vector<Eigen::MatrixXcd> allEigenvectors = std::get<1>(t);
    std::vector<Eigen::Tensor<std::complex<double>,3>> allVelocities;
    return {allEnergies, allEigenvectors, allVelocities};

  } else {

    double delta = 1.0e-8;
    double threshold = 0.000001 / energyRyToEv;// = 1 micro-eV

    std::vector<Eigen::Vector3d> allVectors(numK * 7);

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

      if (k.norm() < 1.0e-6) { // if not zero, velocity is non-analytical
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
