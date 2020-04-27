#include "harmonic.h"
#include "exceptions.h"
#include "eigen.h"
#include "points.h"

std::tuple<Eigen::VectorXd,Eigen::Tensor<std::complex<double>,3>> HarmonicHamiltonian::diagonalize(Point & point) {
  Eigen::VectorXd energies(1);
  Eigen::Tensor<std::complex<double>,3> eigvecs(1,1,1);
  energies.setZero();
  eigvecs.setZero();
  return {energies, eigvecs};
}


Eigen::Tensor<std::complex<double>,3> HarmonicHamiltonian::diagonalizeVelocity(
				Point & point) {
  Eigen::Tensor<std::complex<double>,3> c(1,1,1);
  c.setZero();
  return c;
}
