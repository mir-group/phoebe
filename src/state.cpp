#include "bandstructure.h"
#include "exceptions.h"

State::State(Point & point_,
		double * energies_,
		long numAtoms_, long numBands_,
		std::complex<double> * velocities_,
		std::complex<double> * eigenvectors_) : point{point_},
		energies{energies_} {
	if ( velocities_ != nullptr ) {
		hasVelocities = true;
		velocities = velocities_;
	}
	if ( eigenvectors_ != nullptr ) {
		hasEigenvectors = true;
		eigenvectors = eigenvectors_;
	}
	numBands = numBands_;
	numAtoms = numAtoms_;
}

Point State::getPoint() {
	return point;
}

double State::getWeight() {
	return point.getWeight();
}

double State::getEnergy(const long & bandIndex, double chemicalPotential) {
	if ( bandIndex >= numBands ) {
		Error e("band index too large in getEnergy" ,1);
	}
	return *(energies+bandIndex) - chemicalPotential;
}

Eigen::VectorXd State::getEnergies(double chemicalPotential) {
	Eigen::VectorXd ens(numBands);
	for ( int i=0; i<numBands; i++ ) {
		ens(i) = *(energies+i) - chemicalPotential;
	}
	return ens;
}

Eigen::Vector3d State::getVelocity(const long & bandIndex) {
	if ( ! hasVelocities ) {
		Error e("State doesn't have velocities" ,1);
	}
	if ( bandIndex >= numBands ) {
		Error e("band index too large in getVelocity" ,1);
	}
	std::complex<double> x;
	Eigen::Vector3d groupVelocity;
	for ( long j=0; j<3; j++ ) {
		long ind = compressIndeces(bandIndex, bandIndex, j, numBands,
				numBands, 3);
		x = *(velocities+ind);
		groupVelocity(j) = real(x);
	}
	return groupVelocity;
}

Eigen::Vector3cd State::getVelocity(const long & bandIndex1,
		const long & bandIndex2) {
	if ( ! hasVelocities ) {
		Error e("State doesn't have velocities" ,1);
	}
	if ( bandIndex1 >= numBands || bandIndex2 >= numBands ) {
		Error e("band index too large in getVelocity" ,1);
	}
	Eigen::Vector3cd velocity;
	for ( long j=0; j<3; j++ ) {
		long ind = compressIndeces(bandIndex1,bandIndex2,j,numBands,
				numBands, 3);
		velocity(j) = *(velocities+ind);
	}
	return velocity;
}

Eigen::Tensor<std::complex<double>,3> State::getVelocities() {
	if ( ! hasEigenvectors ) {
		Error e("State doesn't have velocities" ,1);
	}
	Eigen::Tensor<std::complex<double>,3> vels(3, numAtoms, numBands);
	for ( long ib1=0; ib1<numBands; ib1++ ) {
		for ( long ib2=0; ib2<numBands; ib2++ ) {
			for ( long j=0; j<3; j++ ) {
				long ind = compressIndeces(ib1, ib2, j, numBands, numBands, 3);
				vels(ib1,ib2,j) = *(velocities+ind);
			}
		}
	}
	return vels;
}

Eigen::Tensor<std::complex<double>,3> State::getEigenvectors() {
	if ( ! hasEigenvectors ) {
		Error e("State doesn't have eigenvectors" ,1);
	}
	Eigen::Tensor<std::complex<double>,3> eigs(3, numAtoms, numBands);
	for ( long j=0; j<3; j++ ) {
		for ( long ia=0; ia<numAtoms; ia++ ) {
			for ( long ib=0; ib<numBands; ib++ ) {
				long ind = compressIndeces(j, ia, ib, 3, numAtoms, numBands);
				eigs(j,ia,ib) = *(eigenvectors+ind);
			}
		}
	}
	return eigs;
}
