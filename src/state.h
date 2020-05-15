#ifndef STATE_H
#define STATE_H

#include "points.h"
#include "bandstructure.h"
#include "exceptions.h"
#include "utilities.h"

/** Class containing harmonic information for all bands at a given k-q point.
 * State is a base class. PhState and ElState should be used in the code.
 */

template<typename T>
class State {
public:
	/** Class constructor
	 * @param point: a Point instance with the coordinates in the brillouin
	 * zone.
	 * @param energies: a vector with all the energies at a given point
	 * @param velocities: a complex tensor vector with the velocity operator at
	 * this point. Dimensions(numBands,numBands,3). The diagonal over bands is
	 * the group velocity, the off-diagonal are linked to the dipole operator.
	 */
	State(Point<T> & point_,
			double * energies_,
			long numAtoms_,
			long numBands_,
			std::complex<double> * velocities_=nullptr,
			std::complex<double> * eigenvectors_=nullptr);

	/** get the wavevector (Point object)
	 * @return point: a Point object.
	 */
	Point<T> getPoint();

	/** get the energy of a single band
	 * @param bandIndex: integer from 0 to numBands-1
	 * @return energy: Bloch state energy in rydbergs.
	 */
	double getEnergy(const long & bandIndex, double chemicalPotential = 0.);

	/** get all energies at a given point
	 * @return energies: a vector of energies in rydbergs for all bands present
	 */
	Eigen::VectorXd getEnergies(double chemicalPotential = 0.);

	/** get the group velocity of a single Bloch state.
	 * @param bandIndex: integer from 0 to numBands-1.
	 * @return velocity: the 3d-vector of the group velocity.
	 */
	Eigen::Vector3d getVelocity(const long & bandIndex);

	/** get the off-diagonal velocity operator of two Bloch states.
	 * @param bandIndex1: bloch index of the bra state,
	 * integer from 0 to numBands-1.
	 * @param bandIndex2: bloch index of the ket state,
	 * integer from 0 to numBands-1.
	 * @return velocity: the 3d-vector of the velocity.
	 */
	Eigen::Vector3cd getVelocity(const long & bandIndex1,
			const long & bandIndex2);

	/** get the group velocities of all bands for given k/q point
	 * @return groupVelocities: double matrix of size (numBands,3) with the
	 * group velocities.
	 */
	Eigen::MatrixXd getGroupVelocities();

	/** get the velocities of all bands for given k/q point
	 * @return Velocities: a complex tensor of dimensions (numBands,numBands,3)
	 * with the velocity operator.
	 */
	Eigen::Tensor<std::complex<double>,3> getVelocities();

	/** get the weight of the k/q point. Used for integrations over the
	 * brillouin zone with an irreducible mesh of points.
	 * @return weight: a vector of size (numBands).
	 */
	double getWeight();

	/** get the eigenvectors for the current Point
	 * @return eigenvectors: a complex tensor of size (3,numAtoms,numBands)
	 * with the phonon eigenvectors. Returns zeros if the eigenvectors are not
	 * set.
	 */
	Eigen::Tensor<std::complex<double>,3> getEigenvectors();
protected:
	// pointers to the bandstructure, I don't want to duplicate storage here
	Point<T> & point;
	double * energies;
	long numBands;
	long numAtoms;
	std::complex<double> * velocities = nullptr;
	std::complex<double> * eigenvectors = nullptr;
	bool hasVelocities = false;
	bool hasEigenvectors = false;
};

template<typename T>
State<T>::State(Point<T> & point_,
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

template<typename T>
Point<T> State<T>::getPoint() {
	return point;
}

template<typename T>
double State<T>::getWeight() {
	return point.getWeight();
}

template<typename T>
double State<T>::getEnergy(const long & bandIndex, double chemicalPotential) {
	if ( bandIndex >= numBands ) {
		Error e("band index too large in getEnergy" ,1);
	}
	return *(energies+bandIndex) - chemicalPotential;
}

template<typename T>
Eigen::VectorXd State<T>::getEnergies(double chemicalPotential) {
	Eigen::VectorXd ens(numBands);
	for ( int i=0; i<numBands; i++ ) {
		ens(i) = *(energies+i) - chemicalPotential;
	}
	return ens;
}

template<typename T>
Eigen::Vector3d State<T>::getVelocity(const long & bandIndex) {
	if ( ! hasVelocities ) {
		Error e("State doesn't have velocities" ,1);
	}
	if ( bandIndex >= numBands ) {
		Error e("band index too large in getVelocity" ,1);
	}
	std::complex<double> x;
	Eigen::Vector3d groupVelocity;
	for ( long j=0; j<3; j++ ) {
		long ind = compress3Indeces(bandIndex, bandIndex, j, numBands,
				numBands, 3);
		x = *(velocities+ind);
		groupVelocity(j) = real(x);
	}
	return groupVelocity;
}

template<typename T>
Eigen::Vector3cd State<T>::getVelocity(const long & bandIndex1,
		const long & bandIndex2) {
	if ( ! hasVelocities ) {
		Error e("State doesn't have velocities" ,1);
	}
	if ( bandIndex1 >= numBands || bandIndex2 >= numBands ) {
		Error e("band index too large in getVelocity" ,1);
	}
	Eigen::Vector3cd velocity;
	for ( long j=0; j<3; j++ ) {
		long ind = compress3Indeces(bandIndex1,bandIndex2,j,numBands,
				numBands, 3);
		velocity(j) = *(velocities+ind);
	}
	return velocity;
}

template<typename T>
Eigen::Tensor<std::complex<double>,3> State<T>::getVelocities() {
	if ( ! hasEigenvectors ) {
		Error e("State doesn't have velocities" ,1);
	}
	Eigen::Tensor<std::complex<double>,3> vels(3, numAtoms, numBands);
	for ( long ib1=0; ib1<numBands; ib1++ ) {
		for ( long ib2=0; ib2<numBands; ib2++ ) {
			for ( long j=0; j<3; j++ ) {
				long ind = compress3Indeces(ib1, ib2, j, numBands, numBands,3);
				vels(ib1,ib2,j) = *(velocities+ind);
			}
		}
	}
	return vels;
}

template<typename T>
Eigen::Tensor<std::complex<double>,3> State<T>::getEigenvectors() {
	if ( ! hasEigenvectors ) {
		Error e("State doesn't have eigenvectors" ,1);
	}
	Eigen::Tensor<std::complex<double>,3> eigs(3, numAtoms, numBands);
	for ( long j=0; j<3; j++ ) {
		for ( long ia=0; ia<numAtoms; ia++ ) {
			for ( long ib=0; ib<numBands; ib++ ) {
				long ind = compress3Indeces(j, ia, ib, 3, numAtoms, numBands);
				eigs(j,ia,ib) = *(eigenvectors+ind);
			}
		}
	}
	return eigs;
}

#endif
