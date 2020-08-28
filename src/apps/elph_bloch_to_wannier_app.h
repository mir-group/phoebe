#ifndef ELPHBLOCHTOWANAPP_H
#define ELPHBLOCHTOWANAPP_H

#include "app.h"

/** Main driver for the transport calculation
 */
class ElPhBlochToWannierApp : public App {
public:
  void run(Context &context);
  void checkRequirements(Context &context);
protected:
  Eigen::MatrixXd getBravaisVectors(Crystal &crystal,
                                     const Eigen::Vector3i &kMesh);
  Eigen::Tensor<std::complex<double>, 5> blochToWannier(
      const Eigen::MatrixXd &elBravaisVectors,
      const Eigen::MatrixXd &phBravaisVectors,
      const Eigen::Tensor<std::complex<double>, 5> &g_full,
      FullBandStructure &elBandStructure, FullBandStructure &phBandStructure);
};

#endif
