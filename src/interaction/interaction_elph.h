#ifndef ELPHINTERACTION_H
#define ELPHINTERACTION_H

#include <complex>

#include "constants.h"
#include "crystal.h"
#include "eigen.h"
#include "points.h"
#include "utilities.h"

// * Class to calculate the coupling of a transition of the form (el-ph) -> el
class InteractionElPhWan {
 private:
  Eigen::Tensor<std::complex<double>, 5> couplingWannier;
  // numElBands,numElBands,numPhBands,numPhBravaisVectors,numElBravaisVectors);

  Eigen::MatrixXd elBravaisVectors;
  Eigen::VectorXd elBravaisVectorsWeights;
  Eigen::MatrixXd phBravaisVectors;
  Eigen::VectorXd phBravaisVectorsWeights;

  int numPhBands, numElBands, numElBravaisVectors, numPhBravaisVectors;

  std::vector<Eigen::Tensor<double, 3>> cacheCoupling;

  Eigen::Tensor<std::complex<double>, 4> elPhCached;
  Eigen::Vector3d cachedK2;

 public:
  InteractionElPhWan(
      const Eigen::Tensor<std::complex<double>, 5> &couplingWannier_,
      const Eigen::MatrixXd &elBravaisVectors_,
      const Eigen::VectorXd &elBravaisVectorsWeights_,
      const Eigen::MatrixXd &phBravaisVectors_,
      const Eigen::VectorXd &phBravaisVectorsWeights_);  // default constructor
  InteractionElPhWan(const InteractionElPhWan &that);          // copy constructor
  InteractionElPhWan &operator=(const InteractionElPhWan &that);  // assignment op

  void calcCouplingSquared(Eigen::MatrixXcd &el2Eigvec,
                           std::vector<Eigen::MatrixXcd> &el1Eigenvecs,
                           std::vector<Eigen::MatrixXcd> &phEigvecs,
                           Eigen::Vector3d &k2,
                           std::vector<Eigen::Vector3d> &k1s,
                           std::vector<Eigen::Vector3d> &q3s);

  Eigen::Tensor<double, 3> getCouplingSquared(const int &ik1);

  static InteractionElPhWan parse(const std::string &fileName,
                                  Crystal &crystal);
};

#endif
