#include <Eigen/Dense>
#include <Eigen/Core>
#include "exceptions.h"
#include "crystal.h"

Eigen::Matrix3d Crystal::calcReciprocalCell(
		const Eigen::Matrix3d directUnitCell)
{
	Eigen::Matrix3d reciprocalCell = directUnitCell.inverse().transpose();
	return reciprocalCell;
}

void Crystal::setDirectUnitCell(Eigen::Matrix3d directUnitCell_) {
	directUnitCell = directUnitCell_;
	reciprocalUnitCell = calcReciprocalCell(directUnitCell);
}

Eigen::Matrix3d Crystal::getDirectUnitCell() {
	return directUnitCell;
}

Eigen::Matrix3d Crystal::getReciprocalUnitCell() {
	return reciprocalUnitCell;
}

double calcVolume(const Eigen::Matrix3d& directUnitCell, const double alat)
{
	Eigen::Vector3d a1 = directUnitCell.col(0);
	Eigen::Vector3d a2 = directUnitCell.col(1);
	Eigen::Vector3d a3 = directUnitCell.col(2);
	double volume;
	volume = abs( a1.dot(( a2.cross(a3) )) );
	volume+= abs( a2.dot(( a3.cross(a1) )) );
	volume+= abs( a3.dot(( a1.cross(a2) )) );
	volume *= alat * alat * alat / 3.;
	return volume;
}

int Crystal::getNumAtoms() {
	return numAtoms;
}

double Crystal::getLatticeParameter() {
	return alat;
}

double Crystal::getVolumeUnitCell() {
	return volumeUnitCell;
}

Eigen::MatrixXd Crystal::getAtomicPositions() {
	return atomicPositions;
}

Eigen::VectorXi Crystal::getAtomicSpecies() {
	return atomicSpecies;
}

std::vector<std::string> Crystal::getAtomicNames() {
	return atomicNames;
}

Eigen::VectorXd Crystal::getAtomicMasses() {
	return atomicMasses;
}

std::vector<std::string> Crystal::getSpeciesNames() {
	return speciesNames;
}

Eigen::VectorXd Crystal::getSpeciesMasses() {
	return speciesMasses;
}

Crystal::Crystal(double& alat_, Eigen::Matrix3d& directUnitCell_,
		Eigen::MatrixXd& atomicPositions_,
		Eigen::VectorXi& atomicSpecies_,
		std::vector<std::string>& speciesNames_,
		Eigen::VectorXd& speciesMasses_) {
	alat = alat_;
	setDirectUnitCell(directUnitCell_); // sets both direct and reciprocal
	volumeUnitCell = calcVolume(directUnitCell, alat);

	if ( atomicSpecies_.size() != atomicPositions_.rows() ) {
		Error e("atomic species and positions are not aligned", 1);
	}
	if ( atomicPositions_.cols() != 3 ) {
		Error e("atomic positions need three coordinates", 1);
	}
	if ( (int)speciesMasses_.size() != (int)speciesNames_.size() ) {
		Error e("species masses and names are not aligned", 1);
	}


	atomicSpecies = atomicSpecies_;
	atomicPositions = atomicPositions_;
	speciesMasses = speciesMasses_;
	speciesNames = speciesNames_;

	numAtoms = atomicPositions.rows();

	Eigen::VectorXd atomicMasses_(numAtoms);
	std::vector<std::string> atomicNames_(numAtoms);

	for ( int i=0; i<numAtoms; i++ ) {
		atomicMasses_(i) = speciesMasses(atomicSpecies(i));
		atomicNames_[i] = speciesNames[atomicSpecies(i)];
	}
	atomicMasses = atomicMasses_;
	atomicNames = atomicNames_;
}
