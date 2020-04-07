#include <Eigen/Dense>
#include <Eigen/Core>
#include "exceptions.h"
#include "crystal.h"
#include "spglib.h"

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

const Eigen::Matrix3d& Crystal::getReciprocalUnitCell() {
	// note: reciprocalUnitCell is  in units of twoPi
	// i.e. must be multiplied by twoPi
	return reciprocalUnitCell;
}

double calcVolume(const Eigen::Matrix3d& directUnitCell)
{
	Eigen::Vector3d a1 = directUnitCell.row(0);
	Eigen::Vector3d a2 = directUnitCell.row(1);
	Eigen::Vector3d a3 = directUnitCell.row(2);
	double volume;
	volume = abs( a1.dot(( a2.cross(a3) )) );
	volume+= abs( a2.dot(( a3.cross(a1) )) );
	volume+= abs( a3.dot(( a1.cross(a2) )) );
	volume /= 3.;
	return volume;
}

int Crystal::getNumAtoms() {
	return numAtoms;
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

std::vector<Eigen::Matrix3d> Crystal::getSymmetryMatrices() {
	return symmetryRotations;
}

int Crystal::getNumSymmetries() {
	return numSymmetries;
}


Crystal::Crystal(Eigen::Matrix3d& directUnitCell_,
		Eigen::MatrixXd& atomicPositions_,
		Eigen::VectorXi& atomicSpecies_,
		std::vector<std::string>& speciesNames_,
		Eigen::VectorXd& speciesMasses_) {
	setDirectUnitCell(directUnitCell_); // sets both direct and reciprocal
	volumeUnitCell = calcVolume(directUnitCell);

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
	numSpecies = atomicSpecies.size();

	Eigen::VectorXd atomicMasses_(numAtoms);
	std::vector<std::string> atomicNames_(numAtoms);

	for ( int i=0; i<numAtoms; i++ ) {
		atomicMasses_(i) = speciesMasses(atomicSpecies(i));
		atomicNames_[i] = speciesNames[atomicSpecies(i)];
	}
	atomicMasses = atomicMasses_;
	atomicNames = atomicNames_;

	// We now look for the symmetry operations of the crystal
	// in this implementation, we rely on spglib

	double latticeSPG[3][3] = {{0.}};
	for ( int i=0; i<3; i++ ) {
		for ( int j=0; j<3; j++ ) {
			latticeSPG[i][j] = directUnitCell(i,j);
		}
	}
	double positionSPG[numAtoms][3] = {{0.}};
	Eigen::Vector3d positionCrystal;
	Eigen::Vector3d positionCartesian;
	for ( int i=0; i<numAtoms; i++ ) {
		for ( int j=0; j<3; j++ ) {
			// note: spglib wants fractional positions
			positionCartesian = atomicPositions.row(i);
			positionCrystal = reciprocalUnitCell * positionCartesian;
			positionSPG[i][j] = positionCrystal(j);
		}
	}
	int typesSPG[numAtoms] = {0};
	for ( int i=0; i<numAtoms; i++ ) {
		typesSPG[i] = atomicSpecies(i) + 1;
	}
	int maxSize = 50;
	int size;
	int rotations[maxSize][3][3];
	double translations[maxSize][3];
	double symprec = 1e-6;
	size = spg_get_symmetry(rotations,
			translations,
			maxSize,
			latticeSPG,
			positionSPG,
			typesSPG,
			numAtoms,
			symprec);

	// here we round the translational symmetries to the digits of symprec
	// this step is used to set exactly to zero some translations
	// which are slightly different from zero due to numerical noise

	int tmp=0;
	for ( int isym=0; isym<size; isym++ ) {
		for ( int i=0; i<3; i++ ) {
			tmp = (int) round(translations[isym][i]/symprec);
			translations[isym][i] = tmp * symprec;
		}
	}

	// now, we discard the symm ops with translations, more difficult to use
	// (easily usable for scalar quantities, less so for vectors)

	numSymmetries = 0;
	for ( int isym=0; isym<size; isym++ ) {
		if ( translations[isym][0] == 0. &&
				translations[isym][1] == 0. &&
				translations[isym][2] == 0. ) {
			numSymmetries += 1;
		}
	}

	// Finally, we store the remaining sym ops as a class property

	std::vector<Eigen::Matrix3d> symmetryRotations_;
	Eigen::Matrix3d thisMatrix;
	thisMatrix.setZero();
	for ( int isym=0; isym<numSymmetries; isym++ ) {
		if ( translations[isym][0] == 0 &&
				translations[isym][1] == 0 &&
				translations[isym][2] == 0 ) {
			for ( int i=0; i<3; i++ ) {
				for ( int j=0; j<3; j++ ) {
					thisMatrix(i,j) = rotations[isym][i][j];
				}
			}
			symmetryRotations_.push_back(thisMatrix);
		}
	}
	symmetryRotations = symmetryRotations_;
}





