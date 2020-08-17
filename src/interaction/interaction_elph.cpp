#include "interaction_elph.h"

#include <fstream>

// default constructor
InteractionElPhWan::InteractionElPhWan(
    const Eigen::Tensor<std::complex<double>, 5> &couplingWannier_,
    const Eigen::MatrixXd &elBravaisVectors_,
    const Eigen::VectorXd &elBravaisVectorsWeights_,
    const Eigen::MatrixXd &phBravaisVectors_,
    const Eigen::VectorXd &phBravaisVectorsWeights_)
    : couplingWannier(couplingWannier_),
      elBravaisVectors(elBravaisVectors_),
      elBravaisVectorsWeights(elBravaisVectorsWeights_),
      phBravaisVectors(phBravaisVectors_),
      phBravaisVectorsWeights(phBravaisVectorsWeights_) {
  numPhBands = couplingWannier.dimension(2);
  numElBands = couplingWannier.dimension(0);
  numPhBravaisVectors = couplingWannier.dimension(3);
  numElBravaisVectors = couplingWannier.dimension(4);
  cachedK2 << -11.,-11.,-11.;
}

// copy constructor
InteractionElPhWan::InteractionElPhWan(const InteractionElPhWan &that)
    : couplingWannier(that.couplingWannier),
      elBravaisVectors(that.elBravaisVectors),
      elBravaisVectorsWeights(that.elBravaisVectorsWeights),
      phBravaisVectors(that.phBravaisVectors),
      phBravaisVectorsWeights(that.phBravaisVectorsWeights),
      numPhBands(that.numPhBands),
      numElBands(that.numElBands),
      numElBravaisVectors(that.numElBravaisVectors),
      numPhBravaisVectors(that.numPhBravaisVectors),
      cacheCoupling(that.cacheCoupling),
      elPhCached(that.elPhCached),
      cachedK2(that.cachedK2) {}

// assignment operator
InteractionElPhWan &InteractionElPhWan::operator=(
    const InteractionElPhWan &that) {
  if (this != &that) {
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
    cachedK2 = that.cachedK2;
  }
  return *this;
}

Eigen::Tensor<double, 3> InteractionElPhWan::getCouplingSquared(
    const int &ik1) {
  return cacheCoupling[ik1];
}

void InteractionElPhWan::calcCouplingSquared(
    Eigen::MatrixXcd &el2Eigvec, std::vector<Eigen::MatrixXcd> &el1Eigenvecs,
    std::vector<Eigen::MatrixXcd> &phEigvecs, Eigen::Vector3d &k2,
    std::vector<Eigen::Vector3d> &k1s, std::vector<Eigen::Vector3d> &q3s) {
  (void)k1s;
  int numLoops = el1Eigenvecs.size();
  cacheCoupling.resize(0);
  cacheCoupling.resize(numLoops);

  // first, fourier transform on k2

  if (k2 != cachedK2) {
    cachedK2 = k2;
    elPhCached.resize(numElBands, numElBands, numPhBands, numPhBravaisVectors);
    elPhCached.setZero();
    // NOTE: this is a loop that must be done only once per every value of k2
    // if k1 is split in batches, the value of ElPhCached must be saved as
    // member
    for (int irEl = 0; irEl < numElBravaisVectors; irEl++) {
      double arg = k2.dot(elBravaisVectors.col(irEl));
      std::complex<double> phase = exp(complexI * arg);
      for (int irPh = 0; irPh < numPhBravaisVectors; irPh++) {
        // As a convention, the first primitive cell in the triplet is
        // restricted to the origin, so the phase for that cell is unity.
        for (int ind3 = 0; ind3 < numPhBands; ind3++) {
          for (int ind2 = 0; ind2 < numElBands; ind2++) {
            for (int ind1 = 0; ind1 < numElBands; ind1++) {
              elPhCached(ind1, ind2, ind3, irPh) +=
                  couplingWannier(ind1, ind2, ind3, irPh, irEl) * phase /
                  elBravaisVectorsWeights(irEl);
            }
          }
        }
      }
    }
  }

  for (int ik = 0; ik < numLoops; ik++) {
    Eigen::Vector3d q3 = q3s[ik];

    Eigen::MatrixXcd ev1 = el1Eigenvecs[ik];
    Eigen::MatrixXcd ev2 = el2Eigvec;
    Eigen::MatrixXcd ev3 = phEigvecs[ik];

    // we allow the number of bands to be different in each direction
    int nb1 = ev1.cols();
    int nb2 = ev2.cols();

    Eigen::Tensor<std::complex<double>, 3> v(nb1, nb2, numPhBands);
    v.setZero();

    // note: tmp* is a tensor over cartesian and atomic indices
    // (whose size coincides with the band numbers)
    Eigen::Tensor<std::complex<double>, 3> tmp(numElBands, numElBands,
                                               numPhBands);
    tmp.setZero();

    for (int irPh = 0; irPh < numPhBravaisVectors; irPh++) {
      // As a convention, the first primitive cell in the triplet is
      // restricted to the origin, so the phase for that cell is unity.
      double arg = q3.dot(phBravaisVectors.col(irPh));
      std::complex<double> phase = exp(complexI * arg);
      for (int iac3 = 0; iac3 < numPhBands; iac3++) {
        for (int iac2 = 0; iac2 < numElBands; iac2++) {
          for (int iac1 = 0; iac1 < numElBands; iac1++) {
            tmp(iac1, iac2, iac3) += elPhCached(iac1, iac2, iac3, irPh) *
                                     phase / phBravaisVectorsWeights(irPh);
          }
        }
      }
    }

    Eigen::Tensor<std::complex<double>, 3> tmp1(nb1, numElBands, numPhBands);
    tmp1.setZero();
    for (int iac3 = 0; iac3 < numPhBands; iac3++) {
      for (int iac2 = 0; iac2 < numElBands; iac2++) {
        for (int ib1 = 0; ib1 < nb1; ib1++) {
          for (int iac1 = 0; iac1 < numElBands; iac1++) {
            tmp1(ib1, iac2, iac3) += tmp(iac1, iac2, iac3) * ev1(iac1, ib1);
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
                tmp1(ib1, iac2, iac3) * std::conj(ev2(iac2, ib2));
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

InteractionElPhWan InteractionElPhWan::parse(const std::string &fileName,
    Crystal &crystal) {

  // Open ElPh file

  int numElBands, numElBravaisVectors, numPhBands, numPhBravaisVectors;
  std::string line;

  std::ifstream infile(fileName);
  if (not infile.is_open()) {
    Error e("ElPh file not found");
  }

  // Read the bravais lattice vectors info for q mesh.
  int tmpI;
  infile >> tmpI >> numPhBravaisVectors;

  Eigen::MatrixXd phBravaisVectors(3, numPhBravaisVectors);
  Eigen::VectorXd phBravaisVectorsWeights(numPhBravaisVectors);
  for (int i = 0; i < numPhBravaisVectors; i++) {
    infile >> phBravaisVectors(0, i) >> phBravaisVectors(1, i) >>
        phBravaisVectors(2, i) >> phBravaisVectorsWeights(i);
  }

  // Read the bravais lattice vectors info for k mesh.
  infile >> tmpI >> numElBravaisVectors;
  Eigen::MatrixXd elBravaisVectors(3, numElBravaisVectors);
  Eigen::VectorXd elBravaisVectorsWeights(numElBravaisVectors);
  for (int i = 0; i < numElBravaisVectors; i++) {
    infile >> elBravaisVectors(0, i) >> elBravaisVectors(1, i) >>
        elBravaisVectors(2, i) >> elBravaisVectorsWeights(i);
  }

  // Read real space matrix elements for el-ph coupling
  infile >> numElBands >> tmpI >> tmpI >> numPhBands >> tmpI;
  Eigen::Tensor<std::complex<double>, 5> couplingWannier(
      numElBands, numElBands, numPhBands, numElBravaisVectors,
      numPhBravaisVectors);
  couplingWannier.setZero();
  for (int i5 = 0; i5 < numPhBravaisVectors; i5++) {
    for (int i4 = 0; i4 < numPhBands; i4++) {
      for (int i3 = 0; i3 < numElBravaisVectors; i3++) {
        for (int i2 = 0; i2 < numElBands; i2++) {
          for (int i1 = 0; i1 < numElBands; i1++) {
            infile >> couplingWannier(i1, i2, i4, i5, i3);
          }
        }
      }
    }
  }
  infile.close();


  Eigen::Matrix3d primitiveCell = crystal.getDirectUnitCell();
  for (int i=0; i<numElBravaisVectors; i++) {
    elBravaisVectors.col(i) = primitiveCell * elBravaisVectors.col(i);
  }
  for (int i=0; i<numPhBravaisVectors; i++) {
    phBravaisVectors.col(i) = primitiveCell * phBravaisVectors.col(i);
  }

  InteractionElPhWan output(couplingWannier, elBravaisVectors,
                            elBravaisVectorsWeights, phBravaisVectors,
                            phBravaisVectorsWeights);

  return output;
}
