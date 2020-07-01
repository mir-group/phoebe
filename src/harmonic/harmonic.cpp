#include "harmonic.h"
#include "exceptions.h"
#include "eigen.h"
#include "points.h"

HarmonicHamiltonian::HarmonicHamiltonian() :
        particle(Particle::phonon) {
}

long HarmonicHamiltonian::getNumBands() {
    return numBands;
}

Particle HarmonicHamiltonian::getParticle() {
    return particle;
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> HarmonicHamiltonian::diagonalize(
        Point &point) {
    (void) point;
    return {Eigen::VectorXd::Zero(1), Eigen::MatrixXcd::Zero(1,1)};
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> HarmonicHamiltonian::diagonalizeFromCoords(
        Eigen::Vector3d &k) {
    (void) k;
    return {Eigen::VectorXd::Zero(1), Eigen::MatrixXcd::Zero(1,1)};
}

Eigen::Tensor<std::complex<double>, 3> HarmonicHamiltonian::diagonalizeVelocity(
        Point &point) {
    (void) point;
    Eigen::Tensor<std::complex<double>, 3> c(1, 1, 1);
    c.setZero();
    return c;
}

FullBandStructure HarmonicHamiltonian::populate(Points &fullPoints,
        bool &withVelocities, bool &withEigenvectors) {
    Error e("base populate not implemented");
    (void) fullPoints;
    (void) withVelocities;
    (void) withEigenvectors;
    FullBandStructure t(numBands, particle, withVelocities, withEigenvectors,
            fullPoints);
    return t;
}
