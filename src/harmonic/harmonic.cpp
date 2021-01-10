#include "harmonic.h"

#include "eigen.h"
#include "exceptions.h"
#include "points.h"

HarmonicHamiltonian::HarmonicHamiltonian() : particle(Particle::phonon) {}

int HarmonicHamiltonian::getNumBands() { return numBands; }

Particle HarmonicHamiltonian::getParticle() { return particle; }

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> HarmonicHamiltonian::diagonalize(
    Point &point) {
  (void)point;
  return {Eigen::VectorXd::Zero(1), Eigen::MatrixXcd::Zero(1, 1)};
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
HarmonicHamiltonian::diagonalizeFromCoordinates(Eigen::Vector3d &k) {
  (void)k;
  return {Eigen::VectorXd::Zero(1), Eigen::MatrixXcd::Zero(1, 1)};
}

Eigen::Tensor<std::complex<double>, 3> HarmonicHamiltonian::diagonalizeVelocity(
    Point &point) {
  Eigen::Vector3d k = point.getCoordinates(Points::cartesianCoordinates);
  return diagonalizeVelocityFromCoordinates(k);
}

Eigen::Tensor<std::complex<double>, 3>
HarmonicHamiltonian::diagonalizeVelocityFromCoordinates(Eigen::Vector3d &coordinates) {
  (void)coordinates;
  Eigen::Tensor<std::complex<double>, 3> c(1, 1, 1);
  c.setZero();
  return c;
}

FullBandStructure HarmonicHamiltonian::populate(Points &fullPoints,
                                                bool &withVelocities,
                                                bool &withEigenvectors,
                                                bool isDistributed) {
  Error("base populate not implemented");
  (void)fullPoints;
  (void)withVelocities;
  (void)withEigenvectors;
  (void)isDistributed;
  FullBandStructure t(numBands, particle, withVelocities, withEigenvectors,
                      fullPoints);
  return t;
}
