#include "phonon_h0.h"

#include <cmath>
#include <complex>
#include <iostream>

#include "constants.h"
#include "eigen.h"
#include "exceptions.h"
#include "points.h"
#include "utilities.h"

PhononH0::PhononH0(Crystal &crystal, const Eigen::Matrix3d &dielectricMatrix_,
                   const Eigen::Tensor<double, 3> &bornCharges_,
                   const Eigen::Tensor<double, 7> &forceConstants_,
                   const std::string &sumRule)
    : particle(Particle::phonon), bornCharges(bornCharges_),
      forceConstants(forceConstants_) {
  // in this section, we save as class properties a few variables
  // that are needed for the diagonalization of phonon frequencies

  directUnitCell = crystal.getDirectUnitCell();
  reciprocalUnitCell = crystal.getReciprocalUnitCell();
  volumeUnitCell = crystal.getVolumeUnitCell();
  atomicSpecies = crystal.getAtomicSpecies();
  speciesMasses = crystal.getSpeciesMasses();
  atomicPositions = crystal.getAtomicPositions();
  dielectricMatrix = dielectricMatrix_;

  if (dielectricMatrix.sum() > 0.) {
    hasDielectric = true; // defaulted to false in *.h file
  }

  Eigen::VectorXi qCoarseGrid_(3);
  qCoarseGrid_(0) = int(forceConstants.dimension(2));
  qCoarseGrid_(1) = int(forceConstants.dimension(3));
  qCoarseGrid_(2) = int(forceConstants.dimension(4));
  qCoarseGrid = qCoarseGrid_;

  numAtoms = crystal.getNumAtoms();
  numBands = numAtoms * 3;

  // now, I initialize an auxiliary set of vectors that are needed
  // for the diagonalization, which are precomputed once and for all.

  Eigen::MatrixXd directUnitCellSup(3, 3);
  directUnitCellSup.col(0) = directUnitCell.col(0) * qCoarseGrid(0);
  directUnitCellSup.col(1) = directUnitCell.col(1) * qCoarseGrid(1);
  directUnitCellSup.col(2) = directUnitCell.col(2) * qCoarseGrid(2);

  nr1Big = 2 * qCoarseGrid(0);
  nr2Big = 2 * qCoarseGrid(1);
  nr3Big = 2 * qCoarseGrid(2);

  PhononH0::wsInit(directUnitCellSup);

  setAcousticSumRule(sumRule);

  reorderDynamicalMatrix();

  if (hasDielectric) { // prebuild terms useful for long range corrections
    double cutoff = gMax * 4.;

    // Estimate of nr1x,nr2x,nr3x generating all vectors up to G^2 < geg
    // Only for dimensions where periodicity is present, e.g. if nr1=1
    // and nr2=1, then the G-vectors run along nr3 only.
    // (useful if system is in vacuum, e.g. 1D or 2D)

    int nr1x = 0;
    if (qCoarseGrid(0) > 1) {
      nr1x = (int)(sqrt(cutoff) / reciprocalUnitCell.col(0).norm()) + 1;
    }
    int nr2x = 0;
    if (qCoarseGrid(1) > 1) {
      nr2x = (int)(sqrt(cutoff) / reciprocalUnitCell.col(1).norm()) + 1;
    }
    int nr3x = 0;
    if (qCoarseGrid(2) > 1) {
      nr3x = (int)(sqrt(cutoff) / reciprocalUnitCell.col(2).norm()) + 1;
    }

    double norm = e2 * fourPi / volumeUnitCell;

    int numG = (2*nr1x+1) * (2*nr2x+1) * (2*nr3x+1);
    gVectors.resize(3,numG);
    {
      int ig = 0;
      for (int m1 = -nr1x; m1 <= nr1x; m1++) {
        for (int m2 = -nr2x; m2 <= nr2x; m2++) {
          for (int m3 = -nr3x; m3 <= nr3x; m3++) {
            gVectors(0, ig) = double(m1) * reciprocalUnitCell(0, 0) + double(m2) * reciprocalUnitCell(0, 1) + double(m3) * reciprocalUnitCell(0, 2);
            gVectors(1, ig) = double(m1) * reciprocalUnitCell(1, 0) + double(m2) * reciprocalUnitCell(1, 1) + double(m3) * reciprocalUnitCell(1, 2);
            gVectors(2, ig) = double(m1) * reciprocalUnitCell(2, 0) + double(m2) * reciprocalUnitCell(2, 1) + double(m3) * reciprocalUnitCell(2, 2);
            ++ig;
          }
        }
      }
    }

    longRangeCorrection1.resize(3,3,numAtoms);
    longRangeCorrection1.setZero();
    for (int ig=0; ig<numG; ++ig) {
      Eigen::Vector3d g = gVectors.col(ig);

      double geg = (g.transpose() * dielectricMatrix * g).value();

      if (0. < geg && geg < 4. * gMax) {
        double normG = norm * exp(-geg * 0.25) / geg;

        Eigen::MatrixXd gZ(3, numAtoms);
        for (int na = 0; na < numAtoms; na++) {
          for (int i : {0, 1, 2}) {
            gZ(i, na) = g(0) * bornCharges(na, 0, i) + g(1) * bornCharges(na, 1, i) + g(2) * bornCharges(na, 2, i);
          }
        }

        Eigen::MatrixXd fnAt = Eigen::MatrixXd::Zero(numAtoms, 3);
        for (int na = 0; na < numAtoms; na++) {
          for (int i : {0, 1, 2}) {
            for (int nb = 0; nb < numAtoms; nb++) {
              double arg =
                  (atomicPositions.row(na) - atomicPositions.row(nb)).dot(g);
              fnAt(na, i) += gZ(i, nb) * cos(arg);
            }
          }
        }
        for (int na = 0; na < numAtoms; na++) {
          for (int j : {0, 1, 2}) {
            for (int i : {0, 1, 2}) {
              longRangeCorrection1(i, j, na) += -normG * gZ(i, na) * fnAt(na, j);
            }
          }
        }
      }
    }
  }
}

// copy constructor
PhononH0::PhononH0(const PhononH0 &that)
    : particle(that.particle), hasDielectric(that.hasDielectric),
      numAtoms(that.numAtoms), numBands(that.numBands),
      directUnitCell(that.directUnitCell),
      reciprocalUnitCell(that.reciprocalUnitCell),
      volumeUnitCell(that.volumeUnitCell), atomicSpecies(that.atomicSpecies),
      speciesMasses(that.speciesMasses), atomicPositions(that.atomicPositions),
      dielectricMatrix(that.dielectricMatrix), bornCharges(that.bornCharges),
      qCoarseGrid(that.qCoarseGrid), forceConstants(that.forceConstants),
      wsCache(that.wsCache), nr1Big(that.nr1Big), nr2Big(that.nr2Big),
      nr3Big(that.nr3Big), numBravaisVectors(that.numBravaisVectors),
      bravaisVectors(that.bravaisVectors), weights(that.weights),
      mat2R(that.mat2R), gVectors(that.gVectors),
      longRangeCorrection1(that.longRangeCorrection1) {}

// copy assignment
PhononH0 &PhononH0::operator=(const PhononH0 &that) {
  if (this != &that) {
    particle = that.particle;
    hasDielectric = that.hasDielectric;
    numAtoms = that.numAtoms;
    numBands = that.numBands;
    directUnitCell = that.directUnitCell;
    reciprocalUnitCell = that.reciprocalUnitCell;
    volumeUnitCell = that.volumeUnitCell;
    atomicSpecies = that.atomicSpecies;
    speciesMasses = that.speciesMasses;
    atomicPositions = that.atomicPositions;
    dielectricMatrix = that.dielectricMatrix;
    bornCharges = that.bornCharges;
    qCoarseGrid = that.qCoarseGrid;
    forceConstants = that.forceConstants;
    wsCache = that.wsCache;
    nr1Big = that.nr1Big;
    nr2Big = that.nr2Big;
    nr3Big = that.nr3Big;
    numBravaisVectors = that.numBravaisVectors;
    bravaisVectors = that.bravaisVectors;
    weights = that.weights;
    mat2R = that.mat2R;

    gVectors = that.gVectors;
    longRangeCorrection1 = that.longRangeCorrection1;
  }
  return *this;
}

int PhononH0::getNumBands() { return numBands; }

Particle PhononH0::getParticle() { return particle; }

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
PhononH0::diagonalize(Point &point) {
  Eigen::Vector3d q = point.getCoordinates(Points::cartesianCoordinates);
  bool withMassScaling = true;
  auto tup = diagonalizeFromCoordinates(q, withMassScaling);
  auto energies = std::get<0>(tup);
  auto eigenvectors = std::get<1>(tup);
  return {energies, eigenvectors};
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
PhononH0::diagonalizeFromCoordinates(Eigen::Vector3d &q) {
  bool withMassScaling = true;
  return diagonalizeFromCoordinates(q, withMassScaling);
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
PhononH0::diagonalizeFromCoordinates(Eigen::Vector3d &q,
                                     const bool &withMassScaling) {
  // to be executed at every q-point to get phonon frequencies and wavevectors
  std::vector<Eigen::Vector3d> qVecs;
  qVecs.push_back(q);
  auto t = internalPopulate(qVecs, withMassScaling);
  std::vector<Eigen::VectorXd> eigvals = std::get<0>(t);
  std::vector<Eigen::MatrixXcd> eigvecs = std::get<1>(t);
  return {eigvals[0], eigvecs[0]};
}

FullBandStructure PhononH0::populate(Points &points, bool &withVelocities,
                                     bool &withEigenvectors,
                                     bool isDistributed) {
  FullBandStructure fullBandStructure(numBands, particle, withVelocities,
                                      withEigenvectors, points, isDistributed);

  std::vector<int> ikIndices = fullBandStructure.getWavevectorIndices();
  int numIkIndices = ikIndices.size();
#pragma omp parallel for
  for (int iik = 0; iik<numIkIndices; ++iik) {
    int ik = ikIndices[iik];

    Point point = fullBandStructure.getPoint(ik);

    auto tup = diagonalize(point);
    auto ens = std::get<0>(tup);
    auto eigenVectors = std::get<1>(tup);
    fullBandStructure.setEnergies(point, ens);

    if (withVelocities) {
      auto velocities = diagonalizeVelocity(point);
      fullBandStructure.setVelocities(point, velocities);
    }
    if (withEigenvectors) {
      fullBandStructure.setEigenvectors(point, eigenVectors);
    }
  }
  return fullBandStructure;
}

void PhononH0::wsInit(const Eigen::MatrixXd &unitCell) {
  const int numMax = 2;
  int index = 0;
  const int numMaxRWS = 200;

  Eigen::MatrixXd tmpResult(3, numMaxRWS);

  for (int ir = -numMax; ir <= numMax; ir++) {
    for (int jr = -numMax; jr <= numMax; jr++) {
      for (int kr = -numMax; kr <= numMax; kr++) {
        for (int i : {0, 1, 2}) {
          tmpResult(i, index) =
              unitCell(i, 0) * ir + unitCell(i, 1) * jr + unitCell(i, 2) * kr;
        }

        if (tmpResult.col(index).squaredNorm() > 1.0e-6) {
          index += 1;
        }
        if (index > numMaxRWS) {
          Error("WSInit > numMaxRWS");
        }
      }
    }
  }
  int numRWS = index;
  Eigen::MatrixXd rws(3, numRWS);
  for (int i = 0; i < numRWS; i++) {
    rws.col(i) = tmpResult.col(i);
  }

  // now, I also prepare the wsCache, which is used to accelerate
  // the shortRange() calculation

  wsCache.resize(2 * nr3Big + 1, 2 * nr2Big + 1, 2 * nr1Big + 1, numAtoms,
                 numAtoms);
  wsCache.setZero();

  for (int na = 0; na < numAtoms; na++) {
    for (int nb = 0; nb < numAtoms; nb++) {
      double total_weight = 0.;

      // sum over r vectors in the super cell - very safe range!

      for (int n1 = -nr1Big; n1 <= nr1Big; n1++) {
        int n1ForCache = n1 + nr1Big;
        for (int n2 = -nr2Big; n2 <= nr2Big; n2++) {
          int n2ForCache = n2 + nr2Big;
          for (int n3 = -nr3Big; n3 <= nr3Big; n3++) {
            int n3ForCache = n3 + nr3Big;

            Eigen::Vector3d r_ws;
            for (int i : {0, 1, 2}) {
              // note that this cell is different from above
              r_ws(i) = double(n1) * directUnitCell(i, 0) +
                        double(n2) * directUnitCell(i, 1) +
                        double(n3) * directUnitCell(i, 2);
              r_ws(i) += atomicPositions(na, i) - atomicPositions(nb, i);
            }

            double x = wsWeight(r_ws, rws);
            wsCache(n3ForCache, n2ForCache, n1ForCache, nb, na) = x;
            total_weight += x;
          }
        }
      }

      if (abs(total_weight - qCoarseGrid(0) * qCoarseGrid(1) * qCoarseGrid(2)) >
          1.0e-8) {
        std::cout << total_weight << " " << qCoarseGrid.prod() << "\n";
        Error("wrong total_weight");
      }
    }
  }
}

double PhononH0::wsWeight(const Eigen::VectorXd &r,
                          const Eigen::MatrixXd &rws) {
  // wsWeight() assigns this weight:
  // - if a point is inside the Wigner-Seitz cell:    weight=1
  // - if a point is outside the WS cell:             weight=0
  // - if a point q is on the border of the WS cell, it finds the number N
  //   of translationally equivalent point q+G  (where G is a lattice vector)
  //   that are also on the border of the cell. Then: weight = 1/N
  //
  // I.e. if a point is on the surface of the WS cell of a cubic lattice
  // it will have weight 1/2; on the vertex of the WS it would be 1/8;
  // the K point of an hexagonal lattice has weight 1/3 and so on.

  // rws: contains the list of nearest neighbor atoms
  // r: the position of the reference point
  // rws.cols(): number of nearest neighbors

  int numREq = 1;

  for (int ir = 0; ir < rws.cols(); ir++) {
    double rrt = r.dot(rws.col(ir));
    double ck = rrt - rws.col(ir).squaredNorm() / 2.;
    if (ck > 1.0e-6) {
      return 0.;
    }
    if (abs(ck) < 1.0e-6) {
      numREq += 1;
    }
  }
  double x = 1. / (double)numREq;
  return x;
}

void PhononH0::reorderDynamicalMatrix() {
  // this part can actually be expensive to execute, so we compute it once
  // at the beginning
  // first, we compute the total number of bravais lattice vectors

  numBravaisVectors = 0;
  for (int n3 = -nr3Big; n3 <= nr3Big; n3++) {
    int n3ForCache = n3 + nr3Big;
    for (int n2 = -nr2Big; n2 <= nr2Big; n2++) {
      int n2ForCache = n2 + nr2Big;
      for (int n1 = -nr1Big; n1 <= nr1Big; n1++) {
        int n1ForCache = n1 + nr1Big;
        for (int nb = 0; nb < numAtoms; nb++) {
          for (int na = 0; na < numAtoms; na++) {
            if (wsCache(n3ForCache, n2ForCache, n1ForCache, nb, na) > 0.) {
              numBravaisVectors += 1;
            }
          }
        }
      }
    }
  }

  // next, we reorder the dynamical matrix along the bravais lattice vectors

  bravaisVectors = Eigen::MatrixXd::Zero(3, numBravaisVectors);
  weights = Eigen::VectorXd::Zero(numBravaisVectors);
  mat2R.resize(3, 3, numAtoms, numAtoms, numBravaisVectors);
  mat2R.setZero();

  int iR = 0;
  for (int n3 = -nr3Big; n3 <= nr3Big; n3++) {
    int n3ForCache = n3 + nr3Big;
    for (int n2 = -nr2Big; n2 <= nr2Big; n2++) {
      int n2ForCache = n2 + nr2Big;
      for (int n1 = -nr1Big; n1 <= nr1Big; n1++) {
        int n1ForCache = n1 + nr1Big;
        for (int nb = 0; nb < numAtoms; nb++) {
          for (int na = 0; na < numAtoms; na++) {
            double weight = wsCache(n3ForCache, n2ForCache, n1ForCache, nb, na);
            if (weight > 0.) {
              weights(iR) = weight;

              Eigen::Vector3d r;
              for (int i : {0, 1, 2}) {
                r(i) = n1 * directUnitCell(i, 0) + n2 * directUnitCell(i, 1) +
                       n3 * directUnitCell(i, 2);
              }
              bravaisVectors.col(iR) = r;

              int m1 = mod((n1 + 1), qCoarseGrid(0));
              if (m1 <= 0) {
                m1 += qCoarseGrid(0);
              }
              int m2 = mod((n2 + 1), qCoarseGrid(1));
              if (m2 <= 0) {
                m2 += qCoarseGrid(1);
              }
              int m3 = mod((n3 + 1), qCoarseGrid(2));
              if (m3 <= 0) {
                m3 += qCoarseGrid(2);
              }
              m1 += -1;
              m2 += -1;
              m3 += -1;

              for (int j : {0, 1, 2}) {
                for (int i : {0, 1, 2}) {
                  mat2R(i, j, na, nb, iR) +=
                      forceConstants(i, j, m1, m2, m3, na, nb);
                }
              }
              iR += 1;
            }
          }
        }
      }
    }
  }

  wsCache.resize(0, 0, 0, 0, 0);
  forceConstants.resize(0, 0, 0, 0, 0, 0, 0);
}

Eigen::Tensor<std::complex<double>, 3>
PhononH0::diagonalizeVelocity(Point &point) {
  Eigen::Vector3d coordinates =
      point.getCoordinates(Points::cartesianCoordinates);
  return diagonalizeVelocityFromCoordinates(coordinates);
}

Eigen::Tensor<std::complex<double>, 3>
PhononH0::diagonalizeVelocityFromCoordinates(Eigen::Vector3d &coordinates) {
  Eigen::Tensor<std::complex<double>, 3> velocity(numBands, numBands, 3);
  velocity.setZero();

  // if we are working at gamma, we set all velocities to zero.
  if (coordinates.norm() < 1.0e-6) {
    return velocity;
  }

  bool withMassScaling = false;

  // get the eigenvectors and the energies of the q-point
  auto tup = diagonalizeFromCoordinates(coordinates, withMassScaling);
  auto energies = std::get<0>(tup);
  auto eigenvectors = std::get<1>(tup);

  // now we compute the velocity operator, diagonalizing the expectation
  // value of the derivative of the dynamical matrix.
  // This works better than doing finite differences on the frequencies.
  double deltaQ = 1.0e-8;
  for (int i : {0, 1, 2}) {
    // define q+ and q- from finite differences.
    Eigen::Vector3d qPlus = coordinates;
    Eigen::Vector3d qMinus = coordinates;
    qPlus(i) += deltaQ;
    qMinus(i) -= deltaQ;

    // diagonalize the dynamical matrix at q+ and q-
    auto tup2 = diagonalizeFromCoordinates(qPlus, withMassScaling);
    auto enPlus = std::get<0>(tup2);
    auto eigPlus = std::get<1>(tup2);
    auto tup1 = diagonalizeFromCoordinates(qMinus, withMassScaling);
    auto enMinus = std::get<0>(tup1);
    auto eigMinus = std::get<1>(tup1);

    // build diagonal matrices with frequencies
    Eigen::MatrixXd enPlusMat(numBands, numBands);
    Eigen::MatrixXd enMinusMat(numBands, numBands);
    enPlusMat.setZero();
    enMinusMat.setZero();
    enPlusMat.diagonal() << enPlus;
    enMinusMat.diagonal() << enMinus;

    // build the dynamical matrix at the two wavevectors
    // since we diagonalized it before, A = M.U.M*
    Eigen::MatrixXcd sqrtDPlus(numBands, numBands);
    sqrtDPlus = eigPlus * enPlusMat * eigPlus.adjoint();
    Eigen::MatrixXcd sqrtDMinus(numBands, numBands);
    sqrtDMinus = eigMinus * enMinusMat * eigMinus.adjoint();

    // now we can build the velocity operator
    Eigen::MatrixXcd der(numBands, numBands);
    der = (sqrtDPlus - sqrtDMinus) / (2. * deltaQ);

    // and to be safe, we reimpose hermiticity
    der = 0.5 * (der + der.adjoint());

    // now we rotate in the basis of the eigenvectors at q.
    der = eigenvectors.adjoint() * der * eigenvectors;

    for (int ib2 = 0; ib2 < numBands; ib2++) {
      for (int ib1 = 0; ib1 < numBands; ib1++) {
        velocity(ib1, ib2, i) = der(ib1, ib2);
      }
    }
  }

  // turns out that the above algorithm has problems with degenerate bands
  // so, we diagonalize the velocity operator in the degenerate subspace,

  for (int ib = 0; ib < numBands; ib++) {
    // first, we check if the band is degenerate, and the size of the
    // degenerate subspace
    int sizeSubspace = 1;
    for (int ib2 = ib + 1; ib2 < numBands; ib2++) {
      // I consider bands degenerate if their frequencies are the same
      // within 0.0001 cm^-1
      if (abs(energies(ib) - energies(ib2)) > 0.0001 / ryToCmm1) {
        break;
      }
      sizeSubspace += 1;
    }

    if (sizeSubspace > 1) {
      Eigen::MatrixXcd subMat(sizeSubspace, sizeSubspace);
      // we have to repeat for every direction
      for (int iCart : {0, 1, 2}) {
        // take the velocity matrix of the degenerate subspace
        for (int j = 0; j < sizeSubspace; j++) {
          for (int i = 0; i < sizeSubspace; i++) {
            subMat(i, j) = velocity(ib + i, ib + j, iCart);
          }
        }

        // reinforce hermiticity
        subMat = 0.5 * (subMat + subMat.adjoint());

        // diagonalize the subMatrix
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigenSolver(subMat);
        const Eigen::MatrixXcd& newEigenVectors = eigenSolver.eigenvectors();

        // rotate the original matrix in the new basis
        // that diagonalizes the subspace.
        subMat = newEigenVectors.adjoint() * subMat * newEigenVectors;

        // reinforce hermiticity
        subMat = 0.5 * (subMat + subMat.adjoint());

        // substitute back
        for (int j = 0; j < sizeSubspace; j++) {
          for (int i = 0; i < sizeSubspace; i++) {
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

Eigen::Vector3i PhononH0::getCoarseGrid() { return qCoarseGrid; }

Eigen::Matrix3d PhononH0::getDielectricMatrix() { return dielectricMatrix; }

Eigen::Tensor<double, 3> PhononH0::getBornCharges() { return bornCharges; }

int PhononH0::getIndexEigenvector(const int &iAt, const int &iPol) const {
  return compress2Indices(iAt, iPol, numAtoms, 3);
}

int PhononH0::getIndexEigenvector(const int &iAt, const int &iPol,
                                  const int &nAtoms) {
  return compress2Indices(iAt, iPol, nAtoms, 3);
}

std::tuple<std::vector<Eigen::VectorXd>,
           std::vector<Eigen::MatrixXcd>> PhononH0::internalPopulate(
    const std::vector<Eigen::Vector3d>& cartesianCoordinates,
    const bool& withMassScaling) {

  int numQ = cartesianCoordinates.size();

  Eigen::MatrixXcd dyn(numBands, numBands);
  dyn.setZero();
  std::vector<Eigen::MatrixXcd> dyns(numQ, dyn);

  // first, compute the dynamical matrix from the short-range part of the
  // force constants, by doing a Fourier transform
  {
    Eigen::MatrixXcd phases(numBravaisVectors, numQ);
#pragma omp parallel for collapse(2)
    for (int iQ = 0; iQ < numQ; ++iQ) {
      for (int iR = 0; iR < numBravaisVectors; iR++) {
        double arg = cartesianCoordinates[iQ].dot(bravaisVectors.col(iR));
        phases(iR, iQ) = exp(-complexI * arg);
      }
    }

#pragma omp parallel for collapse(5) // don't collapse the iR index
    for (int iQ = 0; iQ < numQ; ++iQ) {
      for (int nb = 0; nb < numAtoms; nb++) {
        for (int j : {0, 1, 2}) {
          for (int na = 0; na < numAtoms; na++) {
            for (int i : {0, 1, 2}) {
              for (int iR = 0; iR < numBravaisVectors; iR++) {
                dyns[iQ](na * 3 + i, nb * 3 + j) +=
                    mat2R(i, j, na, nb, iR) * phases(iR,iQ) * weights(iR);
              }
            }
          }
        }
      }
    }
  }

    // then the long range term, which uses some convergence
    // tricks by X. Gonze et al.
    if (hasDielectric) {
      addLongRangeTerm(dyns, cartesianCoordinates);
    }

    // once everything is ready, here we scale by masses and diagonalize

    std::vector<Eigen::VectorXd> allEnergies(numQ);
    std::vector<Eigen::MatrixXcd> allEigenvectors(numQ);
    {
      // impose hermiticity
#pragma omp parallel for
      for (int iQ=0; iQ<numQ; ++iQ) {
        dyns[iQ] = 0.5 * (dyns[iQ] + dyns[iQ].adjoint());
      }

      //  divide by the square root of masses
      for (int jat = 0; jat < numAtoms; jat++) {
        int jType = atomicSpecies(jat);
        for (int iat = 0; iat < numAtoms; iat++) {
          int iType = atomicSpecies(iat);
#pragma omp parallel for
          for (int iQ=0; iQ<numQ; ++iQ) {
            for (int j : {0, 1, 2}) {
              for (int i : {0, 1, 2}) {
                dyns[iQ](iat * 3 + i, jat * 3 + j) /=
                    sqrt(speciesMasses(iType) * speciesMasses(jType));
              }
            }
          }
        }
      }

      // diagonalize
#pragma omp parallel for
      for (int iQ=0; iQ<numQ; ++iQ) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigenSolver(dyns[iQ]);
        Eigen::VectorXd w2 = eigenSolver.eigenvalues();
        Eigen::VectorXd energies(numBands);
        for (int i = 0; i < numBands; i++) {
          if (w2(i) < 0) {
            energies(i) = -sqrt(-w2(i));
          } else {
            energies(i) = sqrt(w2(i));
          }
        }
        allEnergies[iQ] = energies;
        allEigenvectors[iQ] = eigenSolver.eigenvectors();
      }
    }

    if (withMassScaling) {
      // we normalize with the mass.
      // In this way, the Eigenvector matrix U, doesn't satisfy (U^+) * U = I
      // but instead (U^+) * M * U = I, where M is the mass matrix
      // (M is diagonal with atomic masses on the diagonal)
      for (int iAt = 0; iAt < numAtoms; iAt++) {
        int iType = atomicSpecies(iAt);
        for (int iPol : {0, 1, 2}) {
          int ind = getIndexEigenvector(iAt, iPol);
#pragma omp parallel for
          for (int iQ=0; iQ<numQ; ++iQ) {
            for (int iBand = 0; iBand < numBands; iBand++) {
              allEigenvectors[iQ](ind, iBand) /= sqrt(speciesMasses(iType));
            }
          }
        }
      }
    }

    return {allEnergies, allEigenvectors};
}

void PhononH0::addLongRangeTerm(std::vector<Eigen::MatrixXcd> &dyns,
                                const std::vector<Eigen::Vector3d>& qs) {
  // this subroutine is the analogous of rgd_blk in QE
  // compute the rigid-ion (long-range) term for q
  // The long-range term used here, to be added to or subtracted from the
  // dynamical matrices, is exactly the same of the formula introduced in:
  // X. Gonze et al, PRB 50. 13035 (1994) . Only the G-space term is
  // implemented: the Ewald parameter alpha must be large enough to
  // have negligible r-space contribution

  //   complex(DP) :: dyn(3,3,numAtoms,numAtoms) ! dynamical matrix
  //   real(DP) &
  //        q(3),           &! q-vector
  //        sign             ! sign=+/-1.0 ==> add/subtract rigid-ion term

  // alpha, implicitly set to 1, is the Ewald parameter,
  // geg is an estimate of G^2 such that the G-space sum is convergent for alpha
  // very rough estimate: geg/4/alpha > gMax = 14
  // (exp (-14) = 10^-6)

  int numQ = qs.size();

#pragma omp parallel for
  for (int iQ = 0; iQ < numQ; ++iQ) {
    for (int na = 0; na < numAtoms; na++) {
      for (int j : {0, 1, 2}) {
        for (int i : {0, 1, 2}) {
          dyns[iQ](na*3+i, na*3+j) += longRangeCorrection1(i, j, na);
        }
      }
    }
  }

  double norm = e2 * fourPi / volumeUnitCell;

#pragma omp parallel for // not on G, because that involves a reduction
  for (int iQ = 0; iQ < numQ; ++iQ) {
    for (int ig = 0; ig < gVectors.cols(); ++ig) {
      Eigen::Vector3d gq = gVectors.col(ig) + qs[iQ];

      double geg = (gq.transpose() * dielectricMatrix * gq).value();

      if (geg > 0. && geg < 4. * gMax) {
        double normG = norm * exp(-geg * 0.25) / geg;

        Eigen::MatrixXd gqZ(3, numAtoms);
        for (int i : {0, 1, 2}) {
          for (int nb = 0; nb < numAtoms; nb++) {
            gqZ(i, nb) = gq(0) * bornCharges(nb, 0, i) + gq(1) * bornCharges(nb, 1, i) + gq(2) * bornCharges(nb, 2, i);
          }
        }

        Eigen::MatrixXcd phases(numAtoms, numAtoms);
        for (int nb = 0; nb < numAtoms; nb++) {
          for (int na = 0; na < numAtoms; na++) {
            double arg =
                (atomicPositions.row(na) - atomicPositions.row(nb)).dot(gq);
            phases(na, nb) = {cos(arg), sin(arg)};
          }
        }

        for (int nb = 0; nb < numAtoms; nb++) {
          for (int na = 0; na < numAtoms; na++) {
            for (int j : {0, 1, 2}) {
              for (int i : {0, 1, 2}) {
                dyns[iQ](na*3+i, nb*3+j) += normG * phases(na, nb) * gqZ(i, na) * gqZ(j, nb);
              }
            }
          }
        }
      }
    }
  }
}
