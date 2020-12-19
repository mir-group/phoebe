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

  if (dielectricMatrix.sum() > 0.) {
    hasDielectric = true;
  } else {
    hasDielectric = false;
  }

  directUnitCell = crystal.getDirectUnitCell();
  reciprocalUnitCell = crystal.getReciprocalUnitCell();
  volumeUnitCell = crystal.getVolumeUnitCell();
  atomicSpecies = crystal.getAtomicSpecies();
  speciesMasses = crystal.getSpeciesMasses();
  atomicPositions = crystal.getAtomicPositions();
  dielectricMatrix = dielectricMatrix_;

  Eigen::VectorXi qCoarseGrid_(3);
  qCoarseGrid_(0) = forceConstants.dimension(2);
  qCoarseGrid_(1) = forceConstants.dimension(3);
  qCoarseGrid_(2) = forceConstants.dimension(4);
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
      wscache(that.wscache), nr1Big(that.nr1Big), nr2Big(that.nr2Big),
      nr3Big(that.nr3Big), numBravaisVectors(that.numBravaisVectors),
      bravaisVectors(that.bravaisVectors), weights(that.weights),
      mat2R(that.mat2R) {}

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
    wscache = that.wscache;
    nr1Big = that.nr1Big;
    nr2Big = that.nr2Big;
    nr3Big = that.nr3Big;
    numBravaisVectors = that.numBravaisVectors;
    bravaisVectors = that.bravaisVectors;
    weights = that.weights;
    mat2R = that.mat2R;
  }
  return *this;
}

long PhononH0::getNumBands() { return numBands; }

Particle PhononH0::getParticle() { return particle; }

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
PhononH0::diagonalize(Point &point) {
  Eigen::Vector3d q = point.getCoords(Points::cartesianCoords);
  bool withMassScaling = true;
  auto tup = diagonalizeFromCoords(q, withMassScaling);
  auto energies = std::get<0>(tup);
  auto eigenvectors = std::get<1>(tup);
  return {energies, eigenvectors};
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
PhononH0::diagonalizeFromCoords(Eigen::Vector3d &q) {
  bool withMassScaling = true;
  return diagonalizeFromCoords(q, withMassScaling);
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
PhononH0::diagonalizeFromCoords(Eigen::Vector3d &q,
                                const bool &withMassScaling) {
  // to be executed at every q-point to get phonon frequencies and wavevectors

  Eigen::Tensor<std::complex<double>, 4> dyn(3, 3, numAtoms, numAtoms);
  dyn.setZero();

  Eigen::Tensor<std::complex<double>, 4> f_of_q(3, 3, numAtoms, numAtoms);
  f_of_q.setZero();

  // for now, this part is not executed
  if (na_ifc) {
    double qq = q.norm(); // q is the q-point coordinate
    if (abs(qq) < 1.0e-8) {
      qq = 1.;
    }

    Eigen::VectorXd qHat = q / qq;
    nonAnalIFC(qHat, f_of_q);
  }

  // first, the short range term, which is just a Fourier transform
  shortRangeTerm(dyn, q, f_of_q);

  // then the long range term, which uses some convergence
  // tricks by X. Gonze et al.
  if (hasDielectric && !na_ifc) {
    longRangeTerm(dyn, q, +1.);
  }

  // finally, the non-analytic term from Born charges
  if (!loto_2d && na_ifc) {
    Eigen::VectorXd qHat = q;
    if (abs(qHat(0) - round(qHat(0))) <= 1.0e-6 &&
        abs(qHat(1) - round(qHat(1))) <= 1.0e-6 &&
        abs(qHat(2) - round(qHat(2))) <= 1.0e-6) {
      // q = 0 : we need the direction q => 0 for the non-analytic part

      double qq = qHat.norm();
      if (qq != 0.) {
        qHat /= qq;
      }
      nonAnalyticTerm(qHat, dyn);
    }
  }

  // once everything is ready, here we scale by masses and diagonalize

  auto tup = dynDiag(dyn);
  auto energies = std::get<0>(tup);
  auto eigenvectors = std::get<1>(tup);

  if (withMassScaling) {
    // we normalize with the mass.
    // In this way, the Eigenvector matrix U, doesn't satisfy (U^+) * U = I
    // but instead (U^+) * M * U = I, where M is the mass matrix
    // (M is diagonal with atomic masses on the diagonal)
    for (long iBand = 0; iBand < numBands; iBand++) {
      for (long iAt = 0; iAt < numAtoms; iAt++) {
        long iType = atomicSpecies(iAt);
        for (long iPol : {0, 1, 2}) {
          long ind = getIndexEigvec(iAt, iPol);
          eigenvectors(ind, iBand) /= sqrt(speciesMasses(iType));
        }
      }
    }
  }

  return {energies, eigenvectors};
};

FullBandStructure PhononH0::populate(Points &points, bool &withVelocities,
                                     bool &withEigenvectors,
                                     bool isDistributed) {
  FullBandStructure fullBandStructure(numBands, particle, withVelocities,
                                      withEigenvectors, points, isDistributed);

  for (auto ik : fullBandStructure.getWavevectorIndices()) {
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
  const long nx = 2;
  long index = 0;
  const long nrwsx = 200;

  Eigen::MatrixXd tmpResult(3, nrwsx);

  for (int ir = -nx; ir <= nx; ir++) {
    for (int jr = -nx; jr <= nx; jr++) {
      for (int kr = -nx; kr <= nx; kr++) {
        for (int i : {0, 1, 2}) {
          tmpResult(i, index) =
              unitCell(i, 0) * ir + unitCell(i, 1) * jr + unitCell(i, 2) * kr;
        }

        if (tmpResult.col(index).squaredNorm() > 1.0e-6) {
          index += 1;
        }
        if (index > nrwsx) {
          Error e("WSInit > nrwsx", 1);
        }
      }
    }
  }
  long nrws = index;

  Eigen::MatrixXd rws(3, nrws);
  for (int i = 0; i < nrws; i++) {
    rws.col(i) = tmpResult.col(i);
  }

  // now, I also prepare the wscache, which is used to accelerate
  // the shortRange() calculation

  Eigen::Tensor<double, 5> wscache_(2 * nr3Big + 1, 2 * nr2Big + 1,
                                    2 * nr1Big + 1, numAtoms, numAtoms);
  wscache_.setZero();

  for (long na = 0; na < numAtoms; na++) {
    for (long nb = 0; nb < numAtoms; nb++) {
      double total_weight = 0.;

      // sum over r vectors in the super cell - very safe range!

      for (long n1 = -nr1Big; n1 <= nr1Big; n1++) {
        long n1ForCache = n1 + nr1Big;
        for (long n2 = -nr2Big; n2 <= nr2Big; n2++) {
          long n2ForCache = n2 + nr2Big;
          for (long n3 = -nr3Big; n3 <= nr3Big; n3++) {
            long n3ForCache = n3 + nr3Big;

            Eigen::Vector3d r_ws;
            for (long i : {0, 1, 2}) {
              // note that this cell is different from above
              r_ws(i) = double(n1) * directUnitCell(i, 0)
                      + double(n2) * directUnitCell(i, 1)
                      + double(n3) * directUnitCell(i, 2);
              if (frozenPhonon) {
                r_ws(i) =
                    r_ws(i) + atomicPositions(nb, i) - atomicPositions(na, i);
              } else {
                r_ws(i) =
                    r_ws(i) + atomicPositions(na, i) - atomicPositions(nb, i);
              }
            }

            double x = wsWeight(r_ws, rws);
            wscache_(n3ForCache, n2ForCache, n1ForCache, nb, na) = x;
            total_weight += x;
          }
        }
      }

      if (abs(total_weight - qCoarseGrid(0) * qCoarseGrid(1) * qCoarseGrid(2)) >
          1.0e-8) {
        std::cout << total_weight << " "
                  << qCoarseGrid(0) * qCoarseGrid(1) * qCoarseGrid(2) << "\n";
        Error e("wrong total_weight");
      }
    }
  }
  // save as class member
  wscache = wscache_;
}

double PhononH0::wsWeight(const Eigen::VectorXd &r,
                          const Eigen::MatrixXd &rws) {
  // wsweights assigns this weight:
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

  long nreq = 1;

  for (long ir = 0; ir < rws.cols(); ir++) {
    double rrt = r.dot(rws.col(ir));
    double ck = rrt - rws.col(ir).squaredNorm() / 2.;
    if (ck > 1.0e-6) {
      return 0.;
    }
    if (abs(ck) < 1.0e-6) {
      nreq += 1;
    }
  }
  double x = 1. / (double)nreq;
  return x;
}

void PhononH0::longRangeTerm(Eigen::Tensor<std::complex<double>, 4> &dyn,
                             const Eigen::VectorXd &q, const long &sign) {
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

  // alpha is the Ewald parameter, geg is an estimate of G^2
  // such that the G-space sum is convergent for that alpha
  // very rough estimate: geg/4/alpha > gmax = 14
  // (exp (-14) = 10^-6)

  long nr1x, nr2x, nr3x;
  Eigen::VectorXd zag(3), zbg(3), zcg(3), fnat(3);
  Eigen::MatrixXd reff(2, 2);
  double fac, facgd, arg, gp2;
  std::complex<double> facg;

  double r = 0.;
  double gmax = 14.;
  double alpha = 1.;
  double geg = gmax * alpha * 4.;

  // Estimate of nr1x,nr2x,nr3x generating all vectors up to G^2 < geg
  // Only for dimensions where periodicity is present, e.g. if nr1=1
  // and nr2=1, then the G-vectors run along nr3 only.
  // (useful if system is in vacuum, e.g. 1D or 2D)

  if (qCoarseGrid(0) == 1) {
    nr1x = 0;
  } else {
    nr1x = (long)(sqrt(geg) / reciprocalUnitCell.col(0).norm()) + 1;
  }
  if (qCoarseGrid(1) == 1) {
    nr2x = 0;
  } else {
    nr2x = (long)(sqrt(geg) / reciprocalUnitCell.col(1).norm()) + 1;
  }
  if (qCoarseGrid(2) == 1) {
    nr3x = 0;
  } else {
    nr3x = (long)(sqrt(geg) / reciprocalUnitCell.col(2).norm()) + 1;
  }

  if (abs(double(sign)) != 1.) {
    Error e("wrong value for sign", 1);
  }

  if (loto_2d) {
    fac = double(sign) * e2 / volumeUnitCell / reciprocalUnitCell(3, 3);
    reff.setZero();
    for (long i : {0, 1}) {
      for (long j : {0, 1}) {
        reff(i, j) = dielectricMatrix(i, j) * 0.5 /
                     reciprocalUnitCell(3, 3); // (eps)*c/2
      }
    }
    for (long i : {0, 1}) {
      reff(i, i) = reff(i, i) - 0.5 / reciprocalUnitCell(3, 3);
      // (-1)*c/2
    }
  } else {
    fac = double(sign) * e2 * fourPi / volumeUnitCell;
  }

  Eigen::VectorXd g(3);
  std::complex<double> phase;
  for (long m1 = -nr1x; m1 <= nr1x; m1++) {
    for (long m2 = -nr2x; m2 <= nr2x; m2++) {
      for (long m3 = -nr3x; m3 <= nr3x; m3++) {
        g(0) = double(m1) * reciprocalUnitCell(0, 0) + double(m2) * reciprocalUnitCell(0, 1) +
            double(m3) * reciprocalUnitCell(0, 2);
        g(1) = double(m1) * reciprocalUnitCell(1, 0) + double(m2) * reciprocalUnitCell(1, 1) +
            double(m3) * reciprocalUnitCell(1, 2);
        g(2) = double(m1) * reciprocalUnitCell(2, 0) + double(m2) * reciprocalUnitCell(2, 1) +
            double(m3) * reciprocalUnitCell(2, 2);

        if (loto_2d) {
          geg = g.squaredNorm();
          r = 0.;
          gp2 = g(0) * g(0) + g(1) * g(1);
          if (gp2 > 1.0e-8) {
            r = g(0) * reff(0, 0) * g(0) + g(0) * reff(0, 1) * g(1) +
                g(1) * reff(1, 0) * g(0) + g(1) * reff(1, 1) * g(1);
            r = r / gp2;
          }
        } else {
          geg = (g.transpose() * dielectricMatrix * g).value();
        }

        if (geg > 0. && geg / alpha / 4. < gmax) {
          if (loto_2d) {
            facgd =
                fac * exp(-geg / alpha / 4.) / sqrt(geg) / (1. + r * sqrt(geg));
          } else {
            facgd = fac * exp(-geg / alpha / 4.) / geg;
          }

          for (long na = 0; na < numAtoms; na++) {
            for (long i = 0; i < 3; i++) {
              zag(i) = g(0) * bornCharges(na, 0, i) +
                       g(1) * bornCharges(na, 1, i) +
                       g(2) * bornCharges(na, 2, i);
              fnat(i) = 0.;
              for (long nb = 0; nb < numAtoms; nb++) {
                arg =
                    (atomicPositions.row(na) - atomicPositions.row(nb)).dot(g);
                zcg(i) = g(0) * bornCharges(nb, 0, i) +
                         g(1) * bornCharges(nb, 1, i) +
                         g(2) * bornCharges(nb, 2, i);
                fnat(i) += zcg(i) * cos(arg);
              }
            }

            for (long i = 0; i < 3; i++) {
              for (long j = 0; j < 3; j++) {
                dyn(i, j, na, na) += -facgd * zag(i) * fnat(j);
              }
            }
          }
        }

        g += q;

        if (loto_2d) {
          geg = g.squaredNorm();
          r = 0.;
          gp2 = g(0) * g(0) + g(1) * g(1);
          if (gp2 > 1.0e-8) {
            r = g(0) * reff(0, 0) * g(0) + g(0) * reff(0, 1) * g(1) +
                g(1) * reff(1, 0) * g(0) + g(1) * reff(1, 1) * g(1);
            r = r / gp2;
          }
        } else {
          geg = (g.transpose() * dielectricMatrix * g).value();
        }

        if (geg > 0. && geg / alpha / 4. < gmax) {
          if (loto_2d) {
            facgd =
                fac * exp(-geg / alpha / 4.) / sqrt(geg) / (1. + r * sqrt(geg));
          } else {
            facgd = fac * exp(-geg / alpha / 4.) / geg;
          }

          for (long nb = 0; nb < numAtoms; nb++) {
            for (long i = 0; i < 3; i++) {
              zbg(i) = g(0) * bornCharges(nb, 0, i) +
                       g(1) * bornCharges(nb, 1, i) +
                       g(2) * bornCharges(nb, 2, i);
            }
            for (long na = 0; na < numAtoms; na++) {
              for (long i = 0; i < 3; i++) {
                zag(i) = g(0) * bornCharges(na, 0, i) +
                         g(1) * bornCharges(na, 1, i) +
                         g(2) * bornCharges(na, 2, i);
              }
              arg = (atomicPositions.row(na) - atomicPositions.row(nb)).dot(g);
              phase = {cos(arg), sin(arg)};
              facg = facgd * phase;
              for (long i = 0; i < 3; i++) {
                for (long j = 0; j < 3; j++) {
                  dyn(i, j, na, nb) += facg * zag(i) * zbg(j);
                }
              }
            }
          }
        }
      }
    }
  }
}

void PhononH0::nonAnalyticTerm(const Eigen::VectorXd &q,
                               Eigen::Tensor<std::complex<double>, 4> &dyn) {
  // add the nonAnalytical term with macroscopic electric fields
  // numAtoms: number of atoms in the cell (in the super-cell in the case
  //       of a dyn.mat. constructed in the mass approximation)
  // atomicSpecies(na): atom in the original cell corresponding to
  //                atom na in the super-cell
  // dyn(3,3,numAtoms,numAtoms) ! dynamical matrix
  // q(3), & ! polarization vector
  // dielectricMatrix(3,3),              & ! dielectric constant tensor
  // bornCharges(3,3,numAtoms),        & ! effective charges tensor
  // volumeUnitCell                      ! unit cell volume

  double qeq = (q.transpose() * dielectricMatrix * q).value();
  if (qeq < 1.e-8) {
    Warning w("A direction for q was not specified: "
              "TO-LO splitting will be absent");
    return;
  }

  Eigen::VectorXd zag(3), zbg(3);

  long iType, jType;

  for (long it = 0; it < numAtoms; it++) {
    iType = atomicSpecies(it);
    for (long jt = 0; jt < numAtoms; jt++) {
      jType = atomicSpecies(jt);

      for (long i : {0, 1, 2}) {
        zag(i) = q(0) * bornCharges(iType, 0, i) +
                 q(1) * bornCharges(iType, 1, i) +
                 q(2) * bornCharges(iType, 2, i);
        zbg(i) = q(0) * bornCharges(jType, 0, i) +
                 q(1) * bornCharges(jType, 1, i) +
                 q(2) * bornCharges(jType, 2, i);
      }

      for (long i : {0, 1, 2}) {
        for (long j : {0, 1, 2}) {
          dyn(i, j, it, jt) +=
              fourPi * e2 * zag(i) * zbg(j) / qeq / volumeUnitCell;
        }
      }
    }
  }
}

void PhononH0::nonAnalIFC(const Eigen::VectorXd &q,
                          Eigen::Tensor<std::complex<double>, 4> &f_of_q) {
  //     add the non-analytical term with macroscopic electric fields
  // q(3) : polarization vector

  Eigen::VectorXd zag(3), zbg(3); // eff. charges  times g-vector

  if (q.squaredNorm() == 0.) {
    return;
  }

  double qeq = q.transpose() * dielectricMatrix * q;

  if (qeq < 1.0e-8) {
    return;
  }

  long na_blk, nb_blk;

  double factor = fourPi * e2 / qeq / volumeUnitCell /
                  (qCoarseGrid(0) * qCoarseGrid(1) * qCoarseGrid(2));

  for (int na = 0; na < numAtoms; na++) {
    na_blk = atomicSpecies(na);
    for (int nb = 0; nb < numAtoms; nb++) {
      nb_blk = atomicSpecies(nb);
      for (int i : {0, 1, 2}) {
        zag(i) = q(0) * bornCharges(na_blk, 0, i) +
                 q(1) * bornCharges(na_blk, 1, i) +
                 q(2) * bornCharges(na_blk, 2, i);
        zbg(i) = q(0) * bornCharges(nb_blk, 0, i) +
                 q(1) * bornCharges(nb_blk, 1, i) +
                 q(2) * bornCharges(nb_blk, 2, i);
      }
      for (int i : {0, 1, 2}) {
        for (int j : {0, 1, 2}) {
          f_of_q(i, j, na, nb) = factor * zag(i) * zbg(j);
        }
      }
    }
  }
}

void PhononH0::reorderDynamicalMatrix() {
  // this part can actually be expensive to execute, so we compute it once
  // at the beginning
  // first, we compute the total number of bravais lattice vectors

  numBravaisVectors = 0;
  for (int n3 = -nr3Big; n3 < nr3Big; n3++) {
    int n3ForCache = n3 + nr3Big;
    for (int n2 = -nr2Big; n2 < nr2Big; n2++) {
      int n2ForCache = n2 + nr2Big;
      for (int n1 = -nr1Big; n1 < nr1Big; n1++) {
        int n1ForCache = n1 + nr1Big;
        for (int nb = 0; nb < numAtoms; nb++) {
          for (int na = 0; na < numAtoms; na++) {
            if (wscache(n3ForCache, n2ForCache, n1ForCache, nb, na) > 0.) {
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
  for (int n3 = -nr3Big; n3 < nr3Big; n3++) {
    int n3ForCache = n3 + nr3Big;
    for (int n2 = -nr2Big; n2 < nr2Big; n2++) {
      int n2ForCache = n2 + nr2Big;
      for (int n1 = -nr1Big; n1 < nr1Big; n1++) {
        int n1ForCache = n1 + nr1Big;
        for (int nb = 0; nb < numAtoms; nb++) {
          for (int na = 0; na < numAtoms; na++) {
            double weight = wscache(n3ForCache, n2ForCache, n1ForCache, nb, na);
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
              };
              int m2 = mod((n2 + 1), qCoarseGrid(1));
              if (m2 <= 0) {
                m2 += qCoarseGrid(1);
              };
              int m3 = mod((n3 + 1), qCoarseGrid(2));
              if (m3 <= 0) {
                m3 += qCoarseGrid(2);
              };
              m1 += -1;
              m2 += -1;
              m3 += -1;

              for (int i : {0, 1, 2}) {
                for (int j : {0, 1, 2}) {
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

  wscache.resize(0, 0, 0, 0, 0);
  forceConstants.resize(0, 0, 0, 0, 0, 0, 0);
}

void PhononH0::shortRangeTerm(Eigen::Tensor<std::complex<double>, 4> &dyn,
                              const Eigen::VectorXd &q,
                              Eigen::Tensor<std::complex<double>, 4> &f_of_q) {
  // calculates the dynamical matrix at q from the (short-range part of the)
  // force constants21, by doing the Fourier transform of the force constants

  for (int iR = 0; iR < numBravaisVectors; iR++) {
    double weight = weights(iR);
    Eigen::Vector3d r = bravaisVectors.col(iR);
    double arg = q.dot(r);
    std::complex<double> phase = {cos(arg), -sin(arg)};
    for (int nb = 0; nb < numAtoms; nb++) {
      for (int na = 0; na < numAtoms; na++) {
        for (int j : {0, 1, 2}) {
          for (int i : {0, 1, 2}) {
            dyn(i, j, na, nb) +=
                (mat2R(i, j, na, nb, iR) + f_of_q(i, j, na, nb)) *
                phase * weight;
          }
        }
      }
    }
  }
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
PhononH0::dynDiag(Eigen::Tensor<std::complex<double>, 4> &dyn) {
  // diagonalise the dynamical matrix
  // On input:  speciesMasses = masses, in amu
  // On output: w2 = energies, z = displacements

  // fill the two-indices dynamical matrix

  Eigen::MatrixXcd dyn2Tmp(numBands, numBands);
  for (int iat = 0; iat < numAtoms; iat++) {
    for (int jat = 0; jat < numAtoms; jat++) {
      for (int i : {0, 1, 2}) {
        for (int j : {0, 1, 2}) {
          dyn2Tmp(iat * 3 + i, jat * 3 + j) = dyn(i, j, iat, jat);
        }
      }
    }
  }

  // impose hermiticity

  Eigen::MatrixXcd dyn2(numBands, numBands);
  dyn2 = dyn2Tmp + dyn2Tmp.adjoint();
  dyn2 *= 0.5;

  //  divide by the square root of masses

  for (int iat = 0; iat < numAtoms; iat++) {
    int iType = atomicSpecies(iat);
    for (int jat = 0; jat < numAtoms; jat++) {
      int jType = atomicSpecies(jat);
      for (int i : {0, 1, 2}) {
        for (int j : {0, 1, 2}) {
          dyn2(iat * 3 + i, jat * 3 + j) /=
              sqrt(speciesMasses(iType) * speciesMasses(jType));
        }
      }
    }
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(dyn2);
  Eigen::VectorXd w2 = eigensolver.eigenvalues();

  Eigen::VectorXd energies(numBands);
  for (int i = 0; i < numBands; i++) {
    if (w2(i) < 0) {
      energies(i) = -sqrt(-w2(i));
    } else {
      energies(i) = sqrt(w2(i));
    }
  }
  Eigen::MatrixXcd eigenvectors = eigensolver.eigenvectors();

  return {energies, eigenvectors};
};

Eigen::Tensor<std::complex<double>, 3>
PhononH0::diagonalizeVelocity(Point &point) {
  Eigen::Vector3d coords = point.getCoords(Points::cartesianCoords);
  return diagonalizeVelocityFromCoords(coords);
}

Eigen::Tensor<std::complex<double>, 3>
PhononH0::diagonalizeVelocityFromCoords(Eigen::Vector3d &coords) {
  Eigen::Tensor<std::complex<double>, 3> velocity(numBands, numBands, 3);
  velocity.setZero();

  // if we are working at gamma, we set all velocities to zero.
  if (coords.norm() < 1.0e-6) {
    return velocity;
  }

  bool withMassScaling = false;

  // get the eigenvectors and the energies of the q-point
  auto tup = diagonalizeFromCoords(coords, withMassScaling);
  auto energies = std::get<0>(tup);
  auto eigenvectors = std::get<1>(tup);

  // now we compute the velocity operator, diagonalizing the expectation
  // value of the derivative of the dynamical matrix.
  // This works better than doing finite differences on the frequencies.
  double deltaQ = 1.0e-8;
  for (long i = 0; i < 3; i++) {
    // define q+ and q- from finite differences.
    Eigen::Vector3d qPlus = coords;
    Eigen::Vector3d qMinus = coords;
    qPlus(i) += deltaQ;
    qMinus(i) -= deltaQ;

    // diagonalize the dynamical matrix at q+ and q-
    auto tup2 = diagonalizeFromCoords(qPlus, withMassScaling);
    auto enPlus = std::get<0>(tup2);
    auto eigPlus = std::get<1>(tup2);
    auto tup1 = diagonalizeFromCoords(qMinus, withMassScaling);
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
      if (abs(energies(ib) - energies(ib2)) > 0.0001 / ryToCmm1) {
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
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigenSolver(subMat);
        Eigen::MatrixXcd newEigenVectors = eigenSolver.eigenvectors();

        // rotate the original matrix in the new basis
        // that diagonalizes the subspace.
        subMat = newEigenVectors.adjoint() * subMat * newEigenVectors;

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

Eigen::Vector3i PhononH0::getCoarseGrid() { return qCoarseGrid; }

Eigen::Matrix3d PhononH0::getDielectricMatrix() { return dielectricMatrix; }

Eigen::Tensor<double, 3> PhononH0::getBornCharges() { return bornCharges; }

long PhononH0::getIndexEigvec(const long &iAt, const long &iPol) {
  return compress2Indices(iAt, iPol, numAtoms, 3);
}

long PhononH0::getIndexEigvec(const long &iAt, const long &iPol,
                             const long &nAtoms) {
  return compress2Indices(iAt, iPol, nAtoms, 3);
}
