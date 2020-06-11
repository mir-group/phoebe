#ifndef BANDSTRUCTURE_H
#define BANDSTRUCTURE_H

#include "points.h"
#include "state.h"
#include "statistics.h"
#include "exceptions.h"
#include "utilities.h"

class ActiveBandStructure; // forward declaration of friend class

/** FullBandStructure is the class that stores the energies, velocities and
 * eigenvectors of a quasiparticle computed on a set of wavevectors (as defined
 * by Points() ) in the Brillouin zone.
 */
template<typename T>
class FullBandStructure {
private:
	// stores the quasiparticle kind
	Statistics statistics;

	// link to the underlying mesh of points
	T & points; // these may be FullPoints or PathPoints

	// matrices storing the raw data
	Eigen::MatrixXd energies; // size(bands,points)
	Eigen::MatrixXcd velocities; // size(3*bands^2,points)
	Eigen::MatrixXcd eigenvectors; // size(bands^2,points)

	// pointers to the raw data, used to move
	double * rawEnergies = nullptr;
	std::complex<double> * rawVelocities = nullptr;
	std::complex<double> * rawEigenvectors = nullptr;
	// these are integers used to move the pointers in the raw data
	long energiesCols;
	long velocitiesCols;
	long eigenvectorsCols;

	// auxiliary variables
	long numBands = 0;
	long numAtoms = 0;

	bool hasEigenvectors = false;
	bool hasVelocities = false;

	// method to find the index of the kpoint, from its crystal coordinates
	long getIndex(Eigen::Vector3d & pointCoords);

	// ActiveBandStructure, which restricts the FullBandStructure to a subset
	// of wavevectors, needs low-level access to the raw data (for now)
	friend class ActiveBandStructure;

public:
	/** Constructor of the FullBandStructure
	 * @param numBands: an integer with the number of bands in the system
	 * @param statistics: a Statistics object that contains the type of
	 * quasiparticle
	 * @param withVelocities: a boolean to decide whether to store velocities
	 * @param withEigenvectors: a boolean to decide whether to store eigenvecs
	 * @param points: the underlying mesh of wavevectors.
	 */
	FullBandStructure(long numBands_, Statistics & statistics_,
			bool withVelocities, bool withEigenvectors, T & points_);

	/** Copy constructor
	 */
	FullBandStructure(const FullBandStructure & that);

	/** Assignment operator
	 */
	FullBandStructure & operator = (const FullBandStructure & that);

	/** Returns the wavevectors on which the bandstructure is computed.
	 * @return Points: the object representing the Brillouin zone wavevectors.
	 */
	T getPoints();

	/** Returns a wavevector, given a wavevector index.
	 * The wavevector index runs from 0 to numPoints-1
	 */
	Point<T> getPoint(const long & pointIndex);

	/** Returns the total number of k/q-points.
	 * @return numPoints: the total number of wavevectors of the bandStructure.
	 */
	long getNumPoints();

	/** Returns the number of bands.
	 * @return numPoints: the total number of wavevectors of the bandStructure.
	 */
	long getNumBands();

	/** Returns the total number of Bloch states, equal to numPoints*numBands.
	 * @return numStates: the total number of Bloch states in the class.
	 */
	long getNumStates();

	/** Returns a State object.
	 * The state object, defined elsewhere, is a container that holds all bands
	 * eigenvalues, eigenvectors and velocities (if available), at a fixed
	 * wavevector. This is used in particular for the construction of the
	 * scattering operator.
	 * @param point: a Point object containing the desired wavevector
	 * @return State: a State object evaluated at Point.
	 */
	State<T> getState(Point<T> & point);

	/** Returns a State object
	 * Same as getState(Point<T> & point), but the wavevector is identified
	 * with its integer index
	 * @param pointIndex: index of the wavevector, range [0,numPoints[
	 * @return State: a State object evaluated at Point.
	 */
	State<T> getState(const long & pointIndex);

	/** Returns all electronic energies for all wavevectors at fixed band index
	 * Used by the Fourier interpolation of the band structure.
	 * @param bandIndex: index in [0,numBands[ for the quasiparticle band
	 * @return energies: a vector of size numPoints with the QP energies.
	 * Energies are ordered as the underlying wavevectors in Points.
	 */
	Eigen::VectorXd getBandEnergies(long & bandIndex);

	/** Get the Statistics object associated with this class
	 * @return statistics: a Statistics object, describing e.g. whether this
	 * is a phonon or electron bandStructure
	 */
	Statistics getStatistics();

	/** Returns the energy of a quasiparticle from its Bloch index
	 * Used for accessing the bandstructure in the BTE.
	 * @param stateIndex: an integer index in range [0,numStates[
	 * @return energy: the value of the QP energy for that given Bloch index.
	 * Phonon energies are referred to zero, with negative energies being
	 * actually complex phonon frequencies. Electronic energies are not saved
	 * with any particular reference, and should be used together with the
	 * chemical potential computed by StatisticsSweep. By policy, it's in
	 * rydbergs units.
	 */
	const double & getEnergy(long & stateIndex);

	/** Returns the energy of a quasiparticle from its Bloch index
	 * Used for accessing the bandstructure in the BTE.
	 * @param stateIndex: an integer index in range [0,numStates[
	 * @return velocity: a 3d vector with velocity. By policy, we save it in
	 * the cartesian basis and in atomic rydberg units.
	 */
	Eigen::Vector3d getGroupVelocity(long & stateIndex);

	/** Returns the energy of a quasiparticle from its Bloch index
	 * Used for accessing the bandstructure in the BTE.
	 * @param stateIndex: an integer index in range [0,numStates[
	 * @return wavevector: a 3d vector with the wavevector in cartesian
	 * coordinates in units of Bohr^-1.
	 */
	Eigen::Vector3d getWavevector(long & stateIndex);

	/** Method to save quasiparticle eigenvectors inside FullBandStructure().
	 * @param point: a vector of 3 crystal coordinates. The method will look
	 * for the wavevector index.
	 * @param energies: a vector of size (numBands) with the quasiparticle
	 * energies
	 */
	void setEnergies(Eigen::Vector3d & point, Eigen::VectorXd & energies_);

	/** This method overrides setEnergies, but uses a Point object to find the
	 * k-point.
	 */
	void setEnergies(Point<T> & point, Eigen::VectorXd & energies_);

	/** Method to save quasiparticle eigenvectors inside FullBandStructure().
	 * Note that in this case, eigenvectors are passed as a rank-3 tensor, used
	 * e.g. for phonon eigenvectors Z(3,numAtoms,numBands).
	 * @param point: a Point object with the coordinates of the wavevector,
	 * which should come from the same Point class stored in FullBandStructure.
	 * The wavevector index on the mesh is a property of that Point object.
	 * @param eigenvectors: a rank-3 tensor of size (numBands,numAtoms,3)
	 */
	void setEigenvectors(Point<T> & point,
			Eigen::Tensor<std::complex<double>,3> & eigenvectors_);

	/** Method to save quasiparticle eigenvectors inside FullBandStructure().
	 * Note that in this case, eigenvectors are passed as a matrix, which is
	 * the case e.g. for the Wannier interpolation, where the eigenvectors
	 * represent the unitary transformation matrix U for Wannier localization.
	 * @param point: a Point object with the coordinates of the wavevector,
	 * which should come from the same Point class stored in FullBandStructure
	 * @param eigenvectors: a complex matrix of size (numBands,numBands)
	 */
	void setEigenvectors(Point<T> & point, Eigen::MatrixXcd & eigenvectors_);

	/** Saves in the class the velocities computed at a particular point.
	 * @param point: a Point object representing the wavevector where these
	 * velocities have been computed.
	 * @param velocities: a rank-3 tensor of size (numBands,numBands,3)
	 * containing the matrix elements of the velocity operator. Diagonal
	 * elements are the quasiparticle group velocities.
	 */
	void setVelocities(Point<T> & point,
			Eigen::Tensor<std::complex<double>,3> & velocities_);

	/** Builds a Bloch state index, which runs on both wavevector index and
	 * band index. ik runs from 0 to numPoints-1, ib from 0 to numBands-1.
	 * It's used to view the various matrices such as energy as a 1D vector,
	 * and can be used in combination with get() methods.
	 * @param wavevectorIndex: strong-typed index on wavevector
	 * @return stateIndex: integer from 0 to numStates-1=numBands*numPoints-1
	 */
	long getIndex(const WavevectorIndex & ik, const BandIndex & ib);
};

template<typename T>
FullBandStructure<T>::FullBandStructure(long numBands_,
		Statistics & statistics_,
		bool withVelocities, bool withEigenvectors, T & points_) :
			statistics{statistics_}, points(points_) {

	numBands = numBands_;
	numAtoms = numBands_ / 3;

	if ( withVelocities ) {
		hasVelocities = true;
		velocities = Eigen::MatrixXcd::Zero(numBands*numBands*3,
				getNumPoints());
	}

	if ( withEigenvectors ) {
		hasEigenvectors = true;
		eigenvectors = Eigen::MatrixXcd::Zero(3*numAtoms*numBands,
				getNumPoints());
	}

	energies = Eigen::MatrixXd::Zero(numBands,getNumPoints());

	// now, I want to manipulate the Eigen matrices at lower level
	// I create this pointer to data, so I can move it around
	rawEnergies = energies.data();
	if ( hasVelocities ) {
		rawVelocities = velocities.data();
	}
	if ( hasEigenvectors ) {
		rawEigenvectors = eigenvectors.data();
	}

	energiesCols = numBands;
	velocitiesCols = numBands * numBands * 3;
	eigenvectorsCols = numBands * numAtoms * 3;
}

// copy constructor
template<typename T>
FullBandStructure<T>::FullBandStructure(const FullBandStructure & that) :
	statistics(that.statistics), points(that.points),
	energies(that.energies), velocities(that.velocities),
	eigenvectors(that.eigenvectors), rawEnergies(that.rawEnergies),
	rawVelocities(that.rawVelocities), rawEigenvectors(that.rawEigenvectors),
	energiesCols(that.energiesCols), velocitiesCols(that.velocitiesCols),
	eigenvectorsCols(that.eigenvectorsCols), numBands(that.numBands),
	numAtoms(that.numAtoms),
	hasEigenvectors(that.hasEigenvectors), hasVelocities(that.hasVelocities) {
}

template<typename T>
FullBandStructure<T> & FullBandStructure<T>::operator = ( // copy assignment
		const FullBandStructure & that) {
	if ( this != &that ) {
		statistics = that.statistics;
		points = that.points;
		energies = that.energies;
		velocities = that.velocities;
		eigenvectors = that.eigenvectors;
		rawEnergies = that.rawEnergies;
		rawVelocities = that.rawVelocities;
		rawEigenvectors = that.rawEigenvectors;
		energiesCols = that.energiesCols;
		velocitiesCols = that.velocitiesCols;
		eigenvectorsCols = that.eigenvectorsCols;
		numBands = that.numBands;
		numAtoms = that.numAtoms;
		hasEigenvectors = that.hasEigenvectors;
		hasVelocities = that.hasVelocities;
	}
	return *this;
}

template<typename T>
Statistics FullBandStructure<T>::getStatistics() {
	return statistics;
}

template<typename T>
long FullBandStructure<T>::getNumBands() {
	return numBands;
}

template<typename T>
long FullBandStructure<T>::getNumStates() {
	return numBands*getNumPoints();
}

template<typename T>
long FullBandStructure<T>::getNumPoints() {
	return points.getNumPoints();
}

template<typename T>
long FullBandStructure<T>::getIndex(Eigen::Vector3d& pointCoords) {
	return points.getIndex(pointCoords);
}

template<typename T>
Point<T> FullBandStructure<T>::getPoint(const long& pointIndex) {
	return points.getPoint(pointIndex);
}

template<typename T>
const double & FullBandStructure<T>::getEnergy(long & stateIndex) {
	auto [ik,ib] = decompress2Indeces(stateIndex,getNumPoints(),numBands);
	return energies(ib,ik);
}

template<typename T>
Eigen::Vector3d FullBandStructure<T>::getGroupVelocity(long & stateIndex) {
	auto [ik,ib] = decompress2Indeces(stateIndex,getNumPoints(),numBands);
	Eigen::Vector3d vel;
	for ( long i=0; i<3; i++ ) {
		long ind = compress3Indeces(ib,ib,i,numBands,numBands,3);
		vel(i) = velocities(ind,ik).real();
	}
	return vel;
}

template<typename T>
Eigen::Vector3d FullBandStructure<T>::getWavevector(long & stateIndex) {
	auto [ik,ib] = decompress2Indeces(stateIndex,getNumPoints(),numBands);
	return points.getPoint(ik).getCoords(Points::cartesianCoords);
}

template<typename T>
void FullBandStructure<T>::setEnergies(Eigen::Vector3d& coords,
		Eigen::VectorXd& energies_) {
	long ik = getIndex(coords);
	energies.col(ik) = energies_;
}

template<typename T>
void FullBandStructure<T>::setEnergies(Point<T> & point,
		Eigen::VectorXd& energies_) {
	long ik = point.getIndex();
	energies.col(ik) = energies_;
}

template<typename T>
void FullBandStructure<T>::setVelocities(Point<T> & point,
		Eigen::Tensor<std::complex<double>,3>& velocities_) {
	if ( ! hasVelocities ) {
		Error e("FullBandStructure was initialized without velocities",1);
	}
	// we convert from a tensor to a vector (how it's stored in memory)
	Eigen::VectorXcd tmpVelocities_(numBands*numBands*3);
	for ( long i=0; i<numBands; i++ ) {
		for ( long j=0; j<numBands; j++ ) {
			for ( long k=0; k<3; k++ ) {
				// Note: State must know this order of index compression
				long idx = compress3Indeces(i, j, k, numBands, numBands, 3);
				tmpVelocities_(idx) = velocities_(i,j,k);
			}
		}
	}
	long ik = point.getIndex();
	velocities.col(ik) = tmpVelocities_;
}

template<typename T>
void FullBandStructure<T>::setEigenvectors(Point<T> & point,
		Eigen::Tensor<std::complex<double>,3> & eigenvectors_) {
	if ( ! hasEigenvectors ) {
		Error e("FullBandStructure was initialized without eigvecs",1);
	}
	// we convert from a tensor to a vector (how it's stored in memory)
	Eigen::VectorXcd tmp(numBands*numAtoms*3);
	for ( long i=0; i<numBands; i++ ) {
		for ( long j=0; j<numAtoms; j++ ) {
			for ( long k=0; k<3; k++ ) {
				// Note: State must know this order of index compression
				long idx = compress3Indeces(i, j, k, numBands, numAtoms, 3);
				tmp(idx) = eigenvectors_(k,j,i);
			}
		}
	}
	long ik = point.getIndex();
	eigenvectors.col(ik) = tmp;
}

template<typename T>
void FullBandStructure<T>::setEigenvectors(Point<T> & point,
		Eigen::MatrixXcd & eigenvectors_) {
	if ( ! hasEigenvectors ) {
		Error e("FullBandStructure was initialized without eigvecs",1);
	}
	// we convert from a matrix to a vector (how it's stored in memory)
	Eigen::VectorXcd tmp(numBands*numBands);
	for ( long i=0; i<numBands; i++ ) {
		for ( long j=0; j<numBands; j++ ) {
			// Note: State must know this order of index compression
			long idx = compress2Indeces(i, j, numBands, numBands);
			tmp(idx) = eigenvectors_(j,i);
		}
	}
	long ik = point.getIndex();
	eigenvectors.col(ik) = tmp;
}

template<typename T>
State<T> FullBandStructure<T>::getState(Point<T> & point) {
	long pointIndex = point.getIndex();
	State<T> state = getState(pointIndex);
	return state;
}

template<typename T>
State<T> FullBandStructure<T>::getState(const long & pointIndex) {
	Point<T> point = getPoint(pointIndex);
	// we construct the vector by defining begin() and end()
	// note that rawEnergies points at the start of the matrix
	// and pointIndex*energiesCols skips the first pointIndex-1 wavevectors
	// we are assuming column-major ordering!
	// thisEn points at the location with where the next numBands elements have
	// the desired energies for this wavevector.
	double * thisEn;
	thisEn = rawEnergies + pointIndex * energiesCols;

	// note: in some cases these are nullptr
	std::complex<double> * thisVel = nullptr;
	std::complex<double> * thisEig = nullptr;

	if ( hasVelocities ) {
		thisVel = rawVelocities + pointIndex * velocitiesCols;
	}
	if ( hasEigenvectors ) {
		thisEig = rawEigenvectors + pointIndex * eigenvectorsCols;
	}

	State<T> s(point, thisEn, numAtoms, numBands, thisVel, thisEig);
	return s;
}

template<typename T>
Eigen::VectorXd FullBandStructure<T>::getBandEnergies(long & bandIndex) {
	Eigen::VectorXd bandEnergies = energies.row(bandIndex);
	return bandEnergies;
}

template<typename T>
T FullBandStructure<T>::getPoints() {
	return points;
}

template<typename T>
long FullBandStructure<T>::getIndex(const WavevectorIndex & ik,
		const BandIndex & ib) {
	return ik.get() * numBands + ib.get();
}

#endif
