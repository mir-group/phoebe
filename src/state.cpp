#include "bandstructure.h"
#include "exceptions.h"

State::State(Point& point_,
		Eigen::VectorXd& energies_,
		Eigen::Tensor<std::complex<double>,3>& velocities_,
		Eigen::VectorXd& dnde_, Eigen::VectorXd& dndt_) : point{point_} {
	energies = energies_;
	velocities = velocities_;
	dnde = dnde_;
	dndt = dndt_;
	numBands = energies.size();
}

double State::getEnergy(const int bandIndex) {
	return energies(bandIndex);
}

double State::getWeight() {
	return point.getWeight();
}

Eigen::VectorXd State::getEnergies() {
	Eigen::VectorXd ens(numBands);
	for ( int i=0; i<numBands; i++ ) {
		ens(i) = energies(i);
	}
	return ens;
}

Eigen::Vector3d State::getGroupVelocity(const int bandIndex) {
	Eigen::Vector3d groupVelocity;
	for ( int j=0; j<3; j++ ) {
		groupVelocity(j) = velocities(bandIndex,bandIndex,j).real();
	}
	return groupVelocity;
}

Eigen::MatrixXd State::getGroupVelocities() {
	Eigen::VectorXd groupVelocities(numBands,3);
	for ( int i=0; i<numBands; i++ ) {
		for ( int j=0; j<3; j++ ) {
			groupVelocities(i,j) = velocities(i,i,j).real();
		}
	}
	return groupVelocities;
}

Eigen::Tensor<std::complex<double>,3> State::getVelocities() {
	return velocities;
}

Eigen::VectorXd State::getDndt() {
	return dndt;
}

Eigen::VectorXd State::getDnde() {
	return dnde;
}

ElState::ElState(Point& point_,
		Eigen::VectorXd& energies_,
		Eigen::Tensor<std::complex<double>,3>& velocities_,
		Eigen::VectorXd& dnde_,
		Eigen::VectorXd& dndt_) : State(point_, energies_,
				velocities_, dnde_, dndt_) {};

PhState::PhState(Point& point_,
		Eigen::VectorXd& energies_,
		Eigen::Tensor<std::complex<double>,3>& eigenvectors_,
		Eigen::Tensor<std::complex<double>,3>& velocities_,
		Eigen::VectorXd& dnde_,
		Eigen::VectorXd& dndt_)
               : State{point_, energies_, velocities_, dnde_, dndt_} {
	// we augment the base class initialization
	eigenvectors = eigenvectors_;
}

Eigen::Tensor<std::complex<double>,3> PhState::getEigenvectors() {
	return eigenvectors;
}
