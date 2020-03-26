#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <string>
#include <vector>
#include <Eigen/Core>

class Crystal {
private:
	Eigen::Matrix3d directUnitCell;
	Eigen::Matrix3d reciprocalUnitCell;
	double alat;
	double volumeUnitCell;
	int numAtoms;
	Eigen::Matrix3d calcReciprocalCell(const Eigen::Matrix3d directUnitCell);
	void setDirectUnitCell(Eigen::Matrix3d directUnitCell_);

	// vectors running over the number of atoms
	Eigen::MatrixXd atomicPositions;
	Eigen::VectorXi atomicSpecies;
	std::vector<std::string> atomicNames;
	Eigen::VectorXd atomicMasses;

	// vectors running over the number of species
	std::vector<std::string> speciesNames;
	Eigen::VectorXd speciesMasses;

	std::vector<Eigen::Matrix3d> symmetryRotations;
	int numSymmetries;

public:
	Crystal(double& alat_, Eigen::Matrix3d& directUnitCell_,
			Eigen::MatrixXd& atomicPositions_,
			Eigen::VectorXi& atomicSpecies_,
			std::vector<std::string>& speciesNames_,
			Eigen::VectorXd& speciesMasses_);
	//  Setter and getter for all the variables above
	Eigen::Matrix3d getDirectUnitCell();
	Eigen::Matrix3d getReciprocalUnitCell();
	int getNumAtoms();
	double getLatticeParameter();
	double getVolumeUnitCell();

	std::vector<Eigen::Matrix3d> getSymmetryMatrices();
	int getNumSymmetries();

	Eigen::MatrixXd getAtomicPositions();
	Eigen::VectorXi getAtomicSpecies();
	std::vector<std::string> getAtomicNames();
	Eigen::VectorXd getAtomicMasses();
	std::vector<std::string> getSpeciesNames();
	Eigen::VectorXd getSpeciesMasses();
};

#endif

