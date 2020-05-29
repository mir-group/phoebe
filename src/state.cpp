#include "state.h"
#include "exceptions.h"
#include "utilities.h"

//DetachedState::DetachedState() {}

DetachedState::DetachedState(Eigen::Vector3d & point_,
		Eigen::VectorXd & energies_,
		long numAtoms_,
		long numBands_,
		Eigen::MatrixXcd & eigenvectors_,
		Eigen::Tensor<std::complex<double>,3> * velocities_) :
		point(point_), energies(energies_), numAtoms(numAtoms_),
		numBands(numBands_), eigenvectors(eigenvectors_) {

	if ( velocities_ != nullptr ) {
		velocities = *velocities_;
	}

}

DetachedState::DetachedState(const DetachedState & that) : // copy constructor
	point(that.point), energies(that.energies), numAtoms(that.numAtoms),
	numBands(that.numBands), velocities(that.velocities),
	eigenvectors(that.eigenvectors) {
}

DetachedState & DetachedState::operator=(const DetachedState & that) {// assign
	if ( this != &that ) {
		point = that.point;
		energies = that.energies;
		numAtoms = that.numAtoms;
		numBands = that.numBands;
		velocities = that.velocities;
		eigenvectors = that.eigenvectors;
	}
	return *this;
}

Eigen::Vector3d DetachedState::getCoords(const std::string & basis) {
	if ( basis != "cartesian" ) {
		Error e("Basis not supported in DetachedState");
	}
	return point;
}

double DetachedState::getEnergy(const long & bandIndex,
		double chemicalPotential) {
	return energies(bandIndex) - chemicalPotential;
}

Eigen::VectorXd DetachedState::getEnergies(double chemicalPotential) {
	Eigen::VectorXd ens;
	ens = energies;
	for ( int i=0; i<numBands; i++ ) {
		ens(i) = energies(i) - chemicalPotential;
	}
	return ens;
}

Eigen::Vector3d DetachedState::getVelocity(const long & bandIndex) {
	Eigen::Vector3d groupVelocity;
	for ( long j : {0,1,2} ) {
		groupVelocity(j) = velocities(bandIndex,bandIndex,j).real();
	}
	return groupVelocity;
}

Eigen::Vector3cd DetachedState::getVelocity(const long & bandIndex1,
		const long & bandIndex2) {
	Eigen::Vector3cd vel;
	for ( long j=0; j<3; j++ ) {
		vel(j) = velocities(bandIndex1,bandIndex2,j);
	}
	return vel;
}

Eigen::Tensor<std::complex<double>,3> DetachedState::getVelocities() {
	return velocities;
}

void DetachedState::getEigenvectors(
		Eigen::Tensor<std::complex<double>,3> & eigs) {
	// in this case, we have phonon eigenvectors sized (3,numAtoms,numBands)
	eigs = Eigen::Tensor<std::complex<double>,3>(3,numAtoms,numBands);
	for ( long i = 0; i<numBands; i++ ) {
		auto [iat,ic] = decompress2Indeces(i,numAtoms,3);
		for ( long j = 0; j<numBands; j++ ) {
			eigs(ic,iat,j) = eigenvectors(i,j);
		}
	}
}

void DetachedState::getEigenvectors(Eigen::MatrixXcd & eigs) {
	eigs = eigenvectors;
}
