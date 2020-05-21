#ifndef BANDSTRUCTURE_H
#define BANDSTRUCTURE_H

#include "points.h"
#include "state.h"
#include "statistics.h"
#include "exceptions.h"
#include "utilities.h"

class ActiveBandStructure; // forward declaration of friend class

/** note:
 * FullBandStructure uses matrices to store datas, since they are to be
 * computed in parallel
 * ActiveBandStructure instead is meant to be kept in memory by each MPI
 * process, and we store energies in vectors to be aligned with the
 * vectors of populations.
 */

template<typename T>
class FullBandStructure {
private:
	Statistics statistics;

	T & points; // these may be FullPoints or PathPoints

	MatrixXdRowMajor energies;
	MatrixXcdRowMajor velocities;
	MatrixXcdRowMajor eigenvectors;

	double * rawEnergies = nullptr;
	std::complex<double> * rawVelocities = nullptr;
	std::complex<double> * rawEigenvectors = nullptr;

	long energiesRows;
	long velocitiesRows;
	long eigenvectorsRows;

	long numBands = 0;
	long numAtoms = 0;
	bool useIrreducible = false;

	bool hasEigenvectors = false;
	bool hasVelocities = false;

	long getIndex(Eigen::Vector3d & pointCoords);
	friend class ActiveBandStructure;


	Eigen::VectorXcd packXcd(Eigen::Tensor<std::complex<double>,3> a,
			long size1, long size2, long size3);

public:
	FullBandStructure(long numBands_, Statistics & statistics_,
			bool withVelocities, bool withEigenvectors, T & points_);
	FullBandStructure(const FullBandStructure & that);
	FullBandStructure & operator = (const FullBandStructure & that);

	T getPoints();;
	Point<T> getPoint(const long & pointIndex);
	long getNumPoints();
	State<T> getState(Point<T> & point);
	State<T> getState(const long & pointIndex);
	long getNumBands();
	long getNumStates();
	bool hasIrreduciblePoints();
	Eigen::VectorXd getBandEnergies(long & bandIndex);
	Statistics getStatistics();
	double getEnergy(long & stateIndex);
	Eigen::Vector3d getGroupVelocity(long & stateIndex);

	void setEnergies(Eigen::Vector3d & point, Eigen::VectorXd & energies_);
	void setEnergies(Point<T> & point, Eigen::VectorXd & energies_);
	void setEigenvectors(Point<T> & point,
			Eigen::Tensor<std::complex<double>,3> & eigenvectors_);
	void setVelocities(Point<T> & point,
			Eigen::Tensor<std::complex<double>,3> & velocities_);
};

//auxiliary function to work with the velocity
template<typename T>
Eigen::VectorXcd FullBandStructure<T>::packXcd(
		Eigen::Tensor<std::complex<double>,3> a,
		long size1, long size2, long size3) {
	Eigen::VectorXcd b(size1*size2*size3);
	for ( long i=0; i<size1; i++ ) {
		for ( long j=0; j<size2; j++ ) {
			for ( long k=0; k<size3; k++ ) {

				long idx = compress3Indeces(i, j, k, size1, size2, size3);
				b(idx) = a(i,j,k);
			}
		}
	}
	return b;
}

template<typename T>
FullBandStructure<T>::FullBandStructure(long numBands_,
		Statistics & statistics_,
		bool withVelocities, bool withEigenvectors, T & points_) :
			statistics{statistics_}, points(points_) {

	numBands = numBands_;
	numAtoms = numBands_ / 3;

	if ( withVelocities ) {
		hasVelocities = true;
		Eigen::MatrixXcd velocities_(getNumPoints(), numBands*numBands*3);
		velocities_.setZero();
		velocities = velocities_;
	}

	if ( withEigenvectors ) {
		hasVelocities = true;
		Eigen::MatrixXcd eigenvectors_(getNumPoints(),3*numAtoms*numBands);
		eigenvectors_.setZero();
		eigenvectors = eigenvectors_;
	}

	Eigen::MatrixXd energies_(getNumPoints(), numBands);
	energies_.setZero();
	energies = energies_;

	// now, I want to manipulate the Eigen matrices at lower level
	// I create this pointer to data, so I can move it around
	rawEnergies = energies.data();
	if ( hasVelocities ) {
		rawVelocities = velocities.data();
	}
	if ( hasEigenvectors ) {
		rawEigenvectors = eigenvectors.data();
	}

	energiesRows = numBands;
	velocitiesRows = numBands * numBands * 3;
	eigenvectorsRows = numBands * numAtoms * 3;
}

// copy constructor
template<typename T>
FullBandStructure<T>::FullBandStructure(const FullBandStructure & that) :
	statistics(that.statistics), points(that.points),
	energies(that.energies), velocities(that.velocities),
	eigenvectors(that.eigenvectors), rawEnergies(that.rawEnergies),
	rawVelocities(that.rawVelocities), rawEigenvectors(that.rawEigenvectors),
	energiesRows(that.energiesRows), velocitiesRows(that.velocitiesRows),
	eigenvectorsRows(that.eigenvectorsRows), numBands(that.numBands),
	numAtoms(that.numAtoms), useIrreducible(that.useIrreducible),
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
		energiesRows = that.energiesRows;
		velocitiesRows = that.velocitiesRows;
		eigenvectorsRows = that.eigenvectorsRows;
		numBands = that.numBands;
		numAtoms = that.numAtoms;
		useIrreducible = that.useIrreducible;
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
bool FullBandStructure<T>::hasIrreduciblePoints() {
	return useIrreducible;
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
double FullBandStructure<T>::getEnergy(long & stateIndex) {
	auto [ik,ib] = decompress2Indeces(stateIndex,getNumPoints(),numBands);
	return energies(ik,ib);
}

template<typename T>
Eigen::Vector3d FullBandStructure<T>::getGroupVelocity(long & stateIndex) {
	auto [ik,ib] = decompress2Indeces(stateIndex,getNumPoints(),numBands);
	Eigen::Vector3d vel;
	for ( long i=0; i<3; i++ ) {
		long ind = compress3Indeces(ib,ib,i,numBands,numBands,3);
		vel(i) = velocities(ik,ind).real();
	}
	return vel;
}

template<typename T>
void FullBandStructure<T>::setEnergies(Eigen::Vector3d& coords,
		Eigen::VectorXd& energies_) {
	long ik = getIndex(coords);
	energies.row(ik) = energies_;
}

template<typename T>
void FullBandStructure<T>::setEnergies(Point<T> & point,
		Eigen::VectorXd& energies_) {
	long ik = point.getIndex();
	energies.row(ik) = energies_;
}

template<typename T>
void FullBandStructure<T>::setVelocities(Point<T> & point,
		Eigen::Tensor<std::complex<double>,3>& velocities_) {
	if ( ! hasVelocities ) {
		Error e("FullBandStructure was initialized without velocities",1);
	}
	Eigen::VectorXcd tmpVelocities_(numBands*numBands*3);
	tmpVelocities_ = packXcd(velocities_, numBands, numBands, 3);
	long ik = point.getIndex();
	velocities.row(ik) = tmpVelocities_;
}

template<typename T>
void FullBandStructure<T>::setEigenvectors(Point<T> & point,
		Eigen::Tensor<std::complex<double>,3>& eigenvectors_) {
	if ( ! hasEigenvectors ) {
		Error e("FullBandStructure was initialized without eigvecs",1);
	}
	Eigen::VectorXcd tmp = packXcd(eigenvectors_, 3, numAtoms, numBands);
	long ik = point.getIndex();
	eigenvectors.row(ik) = tmp;
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
	double * thisEn;
	thisEn = rawEnergies + pointIndex * energiesRows;

	// note: in some cases these are nullptr
	std::complex<double> * thisVel = nullptr;
	std::complex<double> * thisEig = nullptr;

	if ( hasVelocities ) {
		thisVel = rawVelocities + pointIndex * velocitiesRows;
	}
	if ( hasEigenvectors ) {
		thisEig = rawEigenvectors + pointIndex * eigenvectorsRows;
	}

	State<T> s(point, thisEn, numAtoms, numBands, thisVel, thisEig);
	return s;
}

template<typename T>
Eigen::VectorXd FullBandStructure<T>::getBandEnergies(long & bandIndex) {
	Eigen::VectorXd bandEnergies = energies.col(bandIndex);
	return bandEnergies;
}

template<typename T>
T FullBandStructure<T>::getPoints() {
	return points;
}

#endif
