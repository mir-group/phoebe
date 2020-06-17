#include "state.h"
#include "exceptions.h"
#include "utilities.h"

DetachedState::DetachedState(Eigen::Vector3d &point_,
        Eigen::VectorXd &energies_,
        long numBands1_,
        long numBands2_,
        Eigen::MatrixXcd &eigenvectors_,
        Eigen::Tensor<std::complex<double>, 3> *velocities_) :
        point(point_), energies(energies_), numBands1(numBands1_),
        numBands2(numBands2_), eigenvectors(eigenvectors_) {

    if (velocities_ != nullptr) {
        velocities = *velocities_;
    }

}

DetachedState::DetachedState(const DetachedState &that) : // copy constructor
        point(that.point), energies(that.energies), numBands1(that.numBands1),
        numBands2(that.numBands2), velocities(that.velocities),
        eigenvectors(that.eigenvectors) {
}

DetachedState& DetachedState::operator=(const DetachedState &that) { // assign
    if (this != &that) {
        point = that.point;
        energies = that.energies;
        numBands1 = that.numBands1;
        numBands2 = that.numBands2;
        velocities = that.velocities;
        eigenvectors = that.eigenvectors;
    }
    return *this;
}

Eigen::Vector3d DetachedState::getCoords(const int &basis) {
    if (basis != Points::cartesianCoords) {
        Error e("Basis not supported in DetachedState");
    }
    return point;
}

double DetachedState::getEnergy(const long &bandIndex,
        double chemicalPotential) {
    return energies(bandIndex) - chemicalPotential;
}

Eigen::VectorXd DetachedState::getEnergies(double chemicalPotential) {
    Eigen::VectorXd ens;
    ens = energies;
    for (int i = 0; i < numBands2; i++) {
        ens(i) = energies(i) - chemicalPotential;
    }
    return ens;
}

Eigen::Vector3d DetachedState::getVelocity(const long &bandIndex) {
    Eigen::Vector3d groupVelocity;
    for (long j : { 0, 1, 2 }) {
        groupVelocity(j) = velocities(bandIndex, bandIndex, j).real();
    }
    return groupVelocity;
}

Eigen::Vector3cd DetachedState::getVelocity(const long &bandIndex1,
        const long &bandIndex2) {
    Eigen::Vector3cd vel;
    for (long j = 0; j < 3; j++) {
        vel(j) = velocities(bandIndex1, bandIndex2, j);
    }
    return vel;
}

Eigen::Tensor<std::complex<double>, 3> DetachedState::getVelocities() {
    return velocities;
}

void DetachedState::getEigenvectors(
        Eigen::Tensor<std::complex<double>, 3> &eigs) {
    // in this case, we have phonon eigenvectors sized (3,numAtoms,numBands)
    long numAtoms = numBands1 / 3;
    eigs = Eigen::Tensor<std::complex<double>, 3>(3, numAtoms, numBands2);
    for (long i = 0; i < numBands1; i++) {
        auto [iat,ic] = decompress2Indeces(i,numAtoms,3);
        for (long j = 0; j < numBands2; j++) {
            eigs(ic, iat, j) = eigenvectors(i, j);
        }
    }
}

void DetachedState::getEigenvectors(Eigen::MatrixXcd &eigs) {
    eigs = eigenvectors;
}

//////////////////////////////////////////////

State::State(Point &point_, double *energies_, long numBands1_, long numBands2_,
        std::complex<double> *velocities_, std::complex<double> *eigenvectors_) :
        point(point_), energies { energies_ } {
    if (velocities_ != nullptr) {
        hasVelocities = true;
        velocities = velocities_;
    }
    if (eigenvectors_ != nullptr) {
        hasEigenvectors = true;
        eigenvectors = eigenvectors_;
    }
    numBands1 = numBands1_;
    numBands2 = numBands2_;
}

State::State(const State &that) : // copy constructor
        point(that.point), energies(that.energies), numBands1(that.numBands1),
        numBands2(that.numBands2), velocities(that.velocities),
        eigenvectors(that.eigenvectors), hasVelocities(that.hasVelocities),
        hasEigenvectors(that.hasEigenvectors) {
}

State& State::operator=(const State &that) { // assignment operator
    if (this != &that) {
        point = that.point;
        energies = that.energies;
        numBands1 = that.numBands1;
        numBands2 = that.numBands2;
        velocities = that.velocities;
        eigenvectors = that.eigenvectors;
        hasVelocities = that.hasVelocities;
        hasEigenvectors = that.hasEigenvectors;
    }
    return *this;
}

Point State::getPoint() {
    return point;
}

Eigen::Vector3d State::getCoords(const int &basis) {
    return point.getCoords(basis);
}

double State::getWeight() {
    return point.getWeight();
}

double State::getEnergy(const long &bandIndex, double chemicalPotential) {
    if (bandIndex >= numBands2) {
        Error e("band index too large in getEnergy", 1);
    }
    return *(energies + bandIndex) - chemicalPotential;
}

Eigen::VectorXd State::getEnergies(double chemicalPotential) {
    Eigen::VectorXd ens(numBands2);
    for (int i = 0; i < numBands2; i++) {
        ens(i) = *(energies + i) - chemicalPotential;
    }
    return ens;
}

Eigen::Vector3d State::getVelocity(const long &bandIndex) {
    if (!hasVelocities) {
        Error e("State doesn't have velocities", 1);
    }
    if (bandIndex >= numBands2) {
        Error e("band index too large in getVelocity", 1);
    }
    std::complex<double> x;
    Eigen::Vector3d groupVelocity;
    for (long j = 0; j < 3; j++) {
        long ind = compress3Indeces(bandIndex, bandIndex, j, numBands2,
                numBands2, 3);
        x = *(velocities + ind);
        groupVelocity(j) = real(x);
    }
    return groupVelocity;
}

Eigen::Vector3cd State::getVelocity(const long &bandIndex1,
        const long &bandIndex2) {
    if (!hasVelocities) {
        Error e("State doesn't have velocities", 1);
    }
    if (bandIndex1 >= numBands2 || bandIndex2 >= numBands2) {
        Error e("band index too large in getVelocity", 1);
    }
    Eigen::Vector3cd velocity;
    for (long j = 0; j < 3; j++) {
        long ind = compress3Indeces(bandIndex1, bandIndex2, j, numBands2,
                numBands2, 3);
        velocity(j) = *(velocities + ind);
    }
    return velocity;
}

Eigen::Tensor<std::complex<double>, 3> State::getVelocities() {
    if (!hasVelocities) {
        Error e("State doesn't have velocities", 1);
    }
    Eigen::Tensor<std::complex<double>, 3> vels(numBands2, numBands2, 3);
    for (long ib1 = 0; ib1 < numBands2; ib1++) {
        for (long ib2 = 0; ib2 < numBands2; ib2++) {
            for (long j = 0; j < 3; j++) {
                long ind = compress3Indeces(ib1, ib2, j, numBands2, numBands2,
                        3);
                vels(ib1, ib2, j) = *(velocities + ind);
            }
        }
    }
    return vels;
}

Eigen::MatrixXd State::getGroupVelocities() {
    if (!hasVelocities) {
        Error e("State doesn't have velocities", 1);
    }
    std::complex<double> x;
    Eigen::MatrixXd vels(numBands2, 3);
    for (long ib1 = 0; ib1 < numBands2; ib1++) {
        for (long j = 0; j < 3; j++) {
            long ind = compress3Indeces(ib1, ib1, j, numBands2, numBands2, 3);
            x = *(velocities + ind);
            vels(ib1, j) = x.real();
        }
    }
    return vels;
}

void State::getEigenvectors(Eigen::Tensor<std::complex<double>, 3> &eigs) {
    if (!hasEigenvectors) {
        Error e("State doesn't have eigenvectors", 1);
    }
    long numAtoms = numBands1 / 3;
    Eigen::Tensor<std::complex<double>, 3> eigs_(3, numAtoms, numBands2);
    for (long ib = 0; ib < numBands2; ib++) {
        for (long ia = 0; ia < numAtoms; ia++) {
            for (long ic : { 0, 1, 2 }) {
                long ind = compress3Indeces(ia, ic, ib, numAtoms, 3, numBands2);
                eigs_(ic, ia, ib) = *(eigenvectors + ind);
            }
        }
    }
    eigs = eigs_;
}

void State::getEigenvectors(Eigen::MatrixXcd &eigs) {
    if (!hasEigenvectors) {
        Error e("State doesn't have eigenvectors", 1);
    }
    eigs = Eigen::MatrixXcd::Zero(numBands1, numBands2);
    for (long ib1 = 0; ib1 < numBands1; ib1++) {
        for (long ib2 = 0; ib2 < numBands2; ib2++) {
            long ind = compress2Indeces(ib1, ib2, numBands1, numBands2);
            eigs(ib1, ib2) = *(eigenvectors + ind);
        }
    }
}
