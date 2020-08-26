#include "interaction_elph.h"

#include <fstream>

// default constructor
InteractionElPhWan::InteractionElPhWan(
    Crystal &crystal_,
    const Eigen::Tensor<std::complex<double>, 5> &couplingWannier_,
    const Eigen::MatrixXd &elBravaisVectors_,
    const Eigen::VectorXd &elBravaisVectorsWeights_,
    const Eigen::MatrixXd &phBravaisVectors_,
    const Eigen::VectorXd &phBravaisVectorsWeights_, PhononH0 *phononH0_)
    : crystal(crystal_), phononH0(phononH0_), couplingWannier(couplingWannier_),
      elBravaisVectors(elBravaisVectors_),
      elBravaisVectorsWeights(elBravaisVectorsWeights_),
      phBravaisVectors(phBravaisVectors_),
      phBravaisVectorsWeights(phBravaisVectorsWeights_) {

  numPhBands = couplingWannier.dimension(2);
  numElBands = couplingWannier.dimension(0);
  numPhBravaisVectors = couplingWannier.dimension(3);
  numElBravaisVectors = couplingWannier.dimension(4);
  cachedK1.setZero();

  if (phononH0 != nullptr) {
    Eigen::Matrix3d epsilon = phononH0->getDielectricMatrix();
    if (epsilon.squaredNorm() > 1.0e-10) { // i.e. if epsilon wasn't computed
      if ( crystal.getNumSpecies() > 1 ) { // otherwise polar correction = 0
        usePolarCorrection = true;
      }
    }
  }
}

// copy constructor
InteractionElPhWan::InteractionElPhWan(const InteractionElPhWan &that)
    : crystal(that.crystal), phononH0(that.phononH0),
      couplingWannier(that.couplingWannier),
      elBravaisVectors(that.elBravaisVectors),
      elBravaisVectorsWeights(that.elBravaisVectorsWeights),
      phBravaisVectors(that.phBravaisVectors),
      phBravaisVectorsWeights(that.phBravaisVectorsWeights),
      numPhBands(that.numPhBands), numElBands(that.numElBands),
      numElBravaisVectors(that.numElBravaisVectors),
      numPhBravaisVectors(that.numPhBravaisVectors),
      cacheCoupling(that.cacheCoupling), elPhCached(that.elPhCached),
      cachedK1(that.cachedK1), usePolarCorrection(that.usePolarCorrection) {}

// assignment operator
InteractionElPhWan &
InteractionElPhWan::operator=(const InteractionElPhWan &that) {
  if (this != &that) {
    crystal = that.crystal;
    phononH0 = that.phononH0;
    couplingWannier = that.couplingWannier;
    elBravaisVectors = that.elBravaisVectors;
    elBravaisVectorsWeights = that.elBravaisVectorsWeights;
    phBravaisVectors = that.phBravaisVectors;
    phBravaisVectorsWeights = that.phBravaisVectorsWeights;
    numPhBands = that.numPhBands;
    numElBands = that.numElBands;
    numElBravaisVectors = that.numElBravaisVectors;
    numPhBravaisVectors = that.numPhBravaisVectors;
    cacheCoupling = that.cacheCoupling;
    elPhCached = that.elPhCached;
    cachedK1 = that.cachedK1;
    usePolarCorrection = that.usePolarCorrection;
  }
  return *this;
}

Eigen::Tensor<double, 3>
InteractionElPhWan::getCouplingSquared(const int &ik2) {
  return cacheCoupling[ik2];
}

void InteractionElPhWan::calcCouplingSquared(
    const Eigen::MatrixXcd &el1Eigenvec,
    const std::vector<Eigen::MatrixXcd> &el2Eigenvecs,
    const std::vector<Eigen::MatrixXcd> &phEigvecs, const Eigen::Vector3d &k1,
    const std::vector<Eigen::Vector3d> &k2s,
    const std::vector<Eigen::Vector3d> &q3s) {
  (void)k2s;
  int numLoops = el2Eigenvecs.size();
  cacheCoupling.resize(0);
  cacheCoupling.resize(numLoops);

  // we allow the number of bands to be different in each direction
  Eigen::MatrixXcd ev1 = el1Eigenvec;
  int nb1 = ev1.cols();

  // first, fourier transform on k2

  if (k1 != cachedK1 || elPhCached.size() == 0) {
    cachedK1 = k1;

    Eigen::Tensor<std::complex<double>, 4> tmp;
    tmp.resize(numElBands, numElBands, numPhBands, numPhBravaisVectors);
    tmp.setZero();
    // NOTE: this is a loop that must be done only once per every value of k1
    // if k2 is split in batches, the value of ElPhCached must be saved as
    // member

    std::vector<std::complex<double>> phases;
    for (int irEl = 0; irEl < numElBravaisVectors; irEl++) {
      double arg = k1.dot(elBravaisVectors.col(irEl));
      std::complex<double> phase =
          exp(complexI * arg) / elBravaisVectorsWeights(irEl);
      phases.push_back(phase);
    }

    for (int irEl = 0; irEl < numElBravaisVectors; irEl++) {
      for (int irPh = 0; irPh < numPhBravaisVectors; irPh++) {
        // As a convention, the first primitive cell in the triplet is
        // restricted to the origin, so the phase for that cell is unity.
        for (int ind3 = 0; ind3 < numPhBands; ind3++) {
          for (int ind2 = 0; ind2 < numElBands; ind2++) {
            for (int ind1 = 0; ind1 < numElBands; ind1++) {
              tmp(ind1, ind2, ind3, irPh) +=
                  couplingWannier(ind1, ind2, ind3, irPh, irEl) * phases[irEl];
            }
          }
        }
      }
    }

    elPhCached.resize(nb1, numElBands, numPhBands, numPhBravaisVectors);
    elPhCached.setZero();

    // note: nb1 can be much smaller than numElBands
    // doing this loop here in the caching makes the subsequent loops faster
    for (int irPh = 0; irPh < numPhBravaisVectors; irPh++) {
      for (int ind3 = 0; ind3 < numPhBands; ind3++) {
        for (int ind2 = 0; ind2 < numElBands; ind2++) {
          for (int ib1 = 0; ib1 < nb1; ib1++) {
            for (int ind1 = 0; ind1 < numElBands; ind1++) {
              elPhCached(ib1, ind2, ind3, irPh) +=
                  tmp(ind1, ind2, ind3, irPh) * ev1(ind1, ib1);
            }
          }
        }
      }
    }
  }

  for (int ik = 0; ik < numLoops; ik++) {
    Eigen::Vector3d q3 = q3s[ik];

    Eigen::MatrixXcd ev2 = el2Eigenvecs[ik];
    int nb2 = ev2.cols();
    Eigen::MatrixXcd ev3 = phEigvecs[ik];

    Eigen::Tensor<std::complex<double>, 3> v(nb1, nb2, numPhBands);
    v.setZero();

    // note: tmp* is a tensor over cartesian and atomic indices
    // (whose size coincides with the band numbers)
    Eigen::Tensor<std::complex<double>, 3> tmp(nb1, numElBands, numPhBands);
    tmp.setZero();
    for (int irPh = 0; irPh < numPhBravaisVectors; irPh++) {
      // As a convention, the first primitive cell in the triplet is
      // restricted to the origin, so the phase for that cell is unity.
      double arg = q3.dot(phBravaisVectors.col(irPh));
      std::complex<double> phase = exp(complexI * arg);
      for (int iac3 = 0; iac3 < numPhBands; iac3++) {
        for (int iac2 = 0; iac2 < numElBands; iac2++) {
          for (int ib1 = 0; ib1 < nb1; ib1++) {
            tmp(ib1, iac2, iac3) += elPhCached(ib1, iac2, iac3, irPh) * phase /
                                    phBravaisVectorsWeights(irPh);
          }
        }
      }
    }

    Eigen::Tensor<std::complex<double>, 3> tmp2(nb1, nb2, numPhBands);
    tmp2.setZero();
    for (int ib1 = 0; ib1 < nb1; ib1++) {
      for (int iac3 = 0; iac3 < numPhBands; iac3++) {
        for (int ib2 = 0; ib2 < nb2; ib2++) {
          for (int iac2 = 0; iac2 < numElBands; iac2++) {
            tmp2(ib1, ib2, iac3) +=
                tmp(ib1, iac2, iac3) * std::conj(ev2(iac2, ib2));
          }
        }
      }
    }
    // the last loop is split in two because nb3 is not the same for + and -
    for (int ib2 = 0; ib2 < nb2; ib2++) {
      for (int ib1 = 0; ib1 < nb1; ib1++) {
        for (int ib3 = 0; ib3 < numPhBands; ib3++) {
          for (int iac3 = 0; iac3 < numPhBands; iac3++) {
            v(ib1, ib2, ib3) += tmp2(ib1, ib2, iac3) * ev3(iac3, ib3);
          }
        }
      }
    }

    if (usePolarCorrection && q3.norm() > 1.0e-8) {
      v += getPolarCorrection(q3, ev1, ev2, ev3);
    }

    Eigen::Tensor<double, 3> coupling(nb1, nb2, numPhBands);
    for (int ib3 = 0; ib3 < numPhBands; ib3++) {
      for (int ib2 = 0; ib2 < nb2; ib2++) {
        for (int ib1 = 0; ib1 < nb1; ib1++) {
          coupling(ib1, ib2, ib3) = std::norm(v(ib1, ib2, ib3));
        }
      }
    }
    cacheCoupling[ik] = coupling;
  }
}

Eigen::Tensor<std::complex<double>, 3> InteractionElPhWan::getPolarCorrection(
    const Eigen::Vector3d &q3, const Eigen::MatrixXcd &ev1,
    const Eigen::MatrixXcd &ev2, const Eigen::MatrixXcd &ev3) {
  // doi:10.1103/physrevlett.115.176401, Eq. 4, is implemented here

  // gather variables
  double volume = crystal.getVolumeUnitCell();
  Eigen::Matrix3d reciprocalUnitCell = crystal.getReciprocalUnitCell();
  Eigen::Matrix3d epsilon = phononH0->getDielectricMatrix();
  int numAtoms = crystal.getNumAtoms();
  Eigen::Tensor<double, 3> bornCharges = phononH0->getBornCharges();
  // must be in Bohr
  Eigen::MatrixXd atomicPositions = crystal.getAtomicPositions();
  Eigen::Vector3i qCoarseMesh = phononH0->getCoarseGrid();

  // overlap = <U^+_{b2 k+q}|U_{b1 k}>
  //         = <psi_{b2 k+q}|e^{i(q+G)r}|psi_{b1 k}>
  Eigen::MatrixXcd overlap = ev2.adjoint() * ev1;
  overlap = overlap.transpose(); // matrix size (nb1,nb2)

  // auxiliary terms
  double gMax = 14.;
  double e2 = 2.; // = e^2/4/Pi/eps_0 in atomic units
  std::complex<double> factor = e2 * fourPi / volume * complexI;

  // build a list of (q+G) vectors
  std::vector<Eigen::Vector3d> gVectors; // here we insert all (q+G)
  for (int m1 = -qCoarseMesh(0); m1 <= qCoarseMesh(0); m1++) {
    for (int m2 = -qCoarseMesh(1); m2 <= qCoarseMesh(1); m2++) {
      for (int m3 = -qCoarseMesh(2); m3 <= qCoarseMesh(2); m3++) {
        Eigen::Vector3d gVector;
        gVector << m1, m2, m3;
        gVector = reciprocalUnitCell * gVector;
        gVector += q3;
        gVectors.push_back(gVector);
      }
    }
  }

  Eigen::VectorXcd x(numPhBands);
  x.setZero();
  for (Eigen::Vector3d gVector : gVectors) {
    double qEq = gVector.transpose() * epsilon * gVector;
    if (qEq > 0. && qEq / 4. < gMax) {
      std::complex<double> factor2 = factor * exp(-qEq / 4.) / qEq;
      for (int iAt = 0; iAt < numAtoms; iAt++) {
        double arg = -gVector.dot(atomicPositions.row(iAt));
        std::complex<double> phase = {cos(arg), sin(arg)};
        std::complex<double> factor3 = factor2 * phase;
        for (int iPol : {0, 1, 2}) {
          double gqDotZ = gVector(0) * bornCharges(iAt, 0, iPol) +
                          gVector(1) * bornCharges(iAt, 1, iPol) +
                          gVector(2) * bornCharges(iAt, 2, iPol);
          int k = phononH0->getIndexEigvec(iAt, iPol);
          for (int ib3 = 0; ib3 < numPhBands; ib3++) {
            x(ib3) += factor3 * gqDotZ * ev3(k, ib3);
          }
        }
      }
    }
  }

  Eigen::Tensor<std::complex<double>, 3> v(overlap.rows(),overlap.cols(),numPhBands);
  v.setZero();
  for (int ib3 = 0; ib3 < numPhBands; ib3++) {
    for (int i = 0; i < overlap.rows(); i++) {
      for (int j = 0; j < overlap.cols(); j++) {
        v(i, j, ib3) += x(ib3) * overlap(i, j);
      }
    }
  }
  return v;
}

InteractionElPhWan InteractionElPhWan::parse(const std::string &fileName,
                                             Crystal &crystal,
                                             PhononH0 *phononH0_) {
  if ( mpi->mpiHead()) {
    std::cout << "\n";
    std::cout << "Started parsing of el-ph interaction." << std::endl;
  }

  int numElBands, numElBravaisVectors, numPhBands, numPhBravaisVectors;
  Eigen::MatrixXd phBravaisVectors, elBravaisVectors;
  Eigen::VectorXd phBravaisVectorsWeights, elBravaisVectorsWeights;
  Eigen::Tensor<std::complex<double>, 5> couplingWannier;

  // Open ElPh file
  if (mpi->mpiHead()) {
    std::string line;
    std::ifstream infile(fileName);
    if (not infile.is_open()) {
      Error e("ElPh file not found");
    }

    // Read the bravais lattice vectors info for q mesh.
    int tmpI;
    infile >> tmpI >> numPhBravaisVectors;

    phBravaisVectors.resize(3, numPhBravaisVectors);
    phBravaisVectorsWeights.resize(numPhBravaisVectors);
    for (int i = 0; i < numPhBravaisVectors; i++) {
      infile >> phBravaisVectors(0, i) >> phBravaisVectors(1, i) >>
          phBravaisVectors(2, i) >> phBravaisVectorsWeights(i);
    }

    // Read the bravais lattice vectors info for k mesh.
    infile >> tmpI >> numElBravaisVectors;
    elBravaisVectors.resize(3, numElBravaisVectors);
    elBravaisVectorsWeights.resize(numElBravaisVectors);
    for (int i = 0; i < numElBravaisVectors; i++) {
      infile >> elBravaisVectors(0, i) >> elBravaisVectors(1, i) >>
          elBravaisVectors(2, i) >> elBravaisVectorsWeights(i);
    }

    // Read real space matrix elements for el-ph coupling
    infile >> numElBands >> tmpI >> tmpI >> numPhBands >> tmpI;
    couplingWannier.resize(numElBands, numElBands, numPhBands,
                           numElBravaisVectors, numPhBravaisVectors);
    couplingWannier.setZero();
    for (int i3 = 0; i3 < numElBravaisVectors; i3++) {
      double re, im;
      for (int i4 = 0; i4 < numPhBands; i4++) {
        for (int i5 = 0; i5 < numPhBravaisVectors; i5++) {
          for (int i2 = 0; i2 < numElBands; i2++) {
            for (int i1 = 0; i1 < numElBands; i1++) {
              infile >> re >> im;
              couplingWannier(i2, i1, i4, i5, i3) = {re, im};
            }
          }
        }
      }
    }
    infile.close();
  } // mpiHead done reading file

  mpi->bcast(&numElBands);
  mpi->bcast(&numPhBands);
  mpi->bcast(&numElBravaisVectors);
  mpi->bcast(&numPhBravaisVectors);

  if (!mpi->mpiHead()) { // head already allocated these
    phBravaisVectors.resize(3, numElBravaisVectors);
    phBravaisVectorsWeights.resize(numElBravaisVectors);
    elBravaisVectors.resize(3, numElBravaisVectors);
    elBravaisVectorsWeights.resize(numElBravaisVectors);
    couplingWannier.resize(numElBands, numElBands, numPhBands,
                           numElBravaisVectors, numPhBravaisVectors);
  }
  mpi->bcast(&elBravaisVectors);
  mpi->bcast(&elBravaisVectorsWeights);
  mpi->bcast(&phBravaisVectors);
  mpi->bcast(&phBravaisVectorsWeights);
  mpi->bcast(&couplingWannier);

  if (mpi->mpiHead()) {
    std::cout << "Finished parsing of el-ph interaction." << std::endl;
  }

  Eigen::Matrix3d primitiveCell = crystal.getDirectUnitCell();
  for (int i = 0; i < numElBravaisVectors; i++) {
    elBravaisVectors.col(i) = primitiveCell * elBravaisVectors.col(i);
  }
  for (int i = 0; i < numPhBravaisVectors; i++) {
    phBravaisVectors.col(i) = primitiveCell * phBravaisVectors.col(i);
  }

  InteractionElPhWan output(crystal, couplingWannier, elBravaisVectors,
                            elBravaisVectorsWeights, phBravaisVectors,
                            phBravaisVectorsWeights, phononH0_);

  return output;
}
