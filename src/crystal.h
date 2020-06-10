#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <string>
#include <vector>
#include "eigen.h"

/** Object to store the information on the crystal unit cell,
 *such as atomic positions, crystal lattice vectors, etc...
 * Note that fractional occupancies are not supported.
 */
class Crystal {
private:
	/** utility function to invert the direct unit cell
	 *
	 */
	Eigen::Matrix3d calcReciprocalCell(const Eigen::Matrix3d directUnitCell);

	/** Internal utility to set the crystal unit cell and the reciprocal cell
	 *
	 */
	void setDirectUnitCell(Eigen::Matrix3d directUnitCell_);

	/** These are the internal quantities used to store
	 * - lattice vectors
	 * - reciprocal lattice vectors
	 * - (real space) crystal unit cell volume
	 * - number of atoms in unit cell
	 * - number of atomic species in the crystal
	 * - dimensionality (to log whether we work in 1, 2, or 3D.
	 */
	Eigen::Matrix3d directUnitCell;
	Eigen::Matrix3d reciprocalUnitCell;
	double volumeUnitCell;
	int numAtoms;
	int numSpecies;
	long dimensionality;

	// these vectors/matrices  running over the number of atoms
	Eigen::MatrixXd atomicPositions; // size (numAtoms,3)
	Eigen::VectorXi atomicSpecies; // size (numAtoms)
	std::vector<std::string> atomicNames; // size (numAtoms)
	Eigen::VectorXd atomicMasses; // size (numAtoms)

	// vectors running over the number of species
	std::vector<std::string> speciesNames; // size (numSpecies)
	Eigen::VectorXd speciesMasses; // size (numSpecies)

	// Untested for now
	std::vector<Eigen::Matrix3d> symmetryRotations;
	int numSymmetries;

public:
	/** Class containing the information on the crystal structure
	 * Note: a number of important quantities used throughout the code are
	 * referenced here, so check if it does go out of scope.
	 * @param directUnitCell: a matrix containing the crystal lattice vectors
	 * in Bohr.
	 * @param atomicPositions(numAtoms,3): list of atomic positions in cartesian
	 * coordinates.
	 * @param atomicSpecies(numAtoms): list of integers identifying the species
	 * of each atom.
	 * @param speciesNames(numSpecies): list of names of species. The position should
	 * match the indices of the parameter atomicSpecies.
	 * @param speciesMasses(numSpecies): array with the masses of each species,
	 *  in rydbergs.
	 */
	Crystal(Eigen::Matrix3d& directUnitCell_,
			Eigen::MatrixXd& atomicPositions_,
			Eigen::VectorXi& atomicSpecies_,
			std::vector<std::string>& speciesNames_,
			Eigen::VectorXd& speciesMasses_,
			long & dimensionality_);

	Crystal(); // default empty constructor
	Crystal(const Crystal & obj); // copy constructor
	Crystal & operator=(const Crystal & obj); // assignment operator

	//  Setter and getter for all the variables above

	/** Returns the crystal unit cell in real space, in Bohr.
	 * Lattice vectors are rows (e.g. a0 = cell.row(0), ...)
	 */
	const Eigen::Matrix3d& getDirectUnitCell();

	/** Returns the reciprocal lattice vectors, in units of Bohr^-1/tpi.
	 * Lattice vectors are rows (e.g. a0 = cell.row(0), ...)
	 * must be multiplied by twoPi to be complete
	 */
	const Eigen::Matrix3d& getReciprocalUnitCell();

	/** get the number of atoms in the unit cell
	 *
	 */
	const int& getNumAtoms();

	/** get the volume of the crystal unit cell in Bohr^3
	 * @param dimensionality: returns the volume of the unit cell on a reduced
	 * dimensionality. If 2D, it is ASSUMED that the z-direction is the non
	 * periodic direction. If 1D, it's assumed that z is the periodic direction
	 */
	double getVolumeUnitCell(long dimensionality = 3);

	/** get the symmetry operations of the crystal, in cartesian coordinates.
	 * For the time being, we only retain symmetry operations that don't use a
	 * translation, so that the object returned by this function is simply a
	 * vector of rotations matrices.
	 */
	const std::vector<Eigen::Matrix3d>& getSymmetryMatrices();

	/** get the number of symmetries operations that are used by phoebe.
	 * For the time being, we only retain symmetry operations that don't use a
	 * translation.
	 */
	const int& getNumSymmetries();

	/** get the atomic positions, in cartesian coordinates
	 * we return an array of size (numAtoms,3)
	 */
	const Eigen::MatrixXd& getAtomicPositions();

	/** get the species of each atom in the unit cell.
	 * The atomic species are simply identified by an integer id
	 * The vector returned has size (numAtoms)
	 */
	const Eigen::VectorXi& getAtomicSpecies();

	/** get the vector of atomic names
	 * vector has size (numAtoms)
	 */
	const std::vector<std::string>& getAtomicNames();

	/** get the vector of atomic masses
	 * the vector has size (numAtoms)
	 */
	const Eigen::VectorXd& getAtomicMasses();

	/** get the vector of names of atomic species
	 * the vector has size (numSpecies)
	 */
	const std::vector<std::string>& getSpeciesNames();

	/** get the vector of masses of atomic species, in rydbergs
	 * the vector has size (numSpecies)
	 */
	const Eigen::VectorXd& getSpeciesMasses();

	/** return the dimensionality of the crystal, i.e. a number from 1 to 3
	 *
	 */
	long getDimensionality();

	/** Return the number of atomic species present in the crystal
	 *
	 */
	long getNumSpecies();
};

#endif

