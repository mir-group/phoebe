#ifndef BANDSTRUCTURE_H
#define BANDSTRUCTURE_H

#include "points.h"
#include "state.h"
#include "window.h"
#include "statistics.h"
#include "harmonic.h"
#include "phonon_h0.h"

/** note:
 * FullBandStructure uses matrices to store datas, since theya are to be
 * computed in parallel
 * ActiveBandStructure instead is meant to be kept in memory by each MPI
 * process, and we store energies in vectors to be aligned with the
 * vectors of populations.
 */

class FullBandStructure {
private:
	Statistics statistics;
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

	FullPoints & fullPoints;
//	FullPoints * fullPoints = nullptr;
//	IrreduciblePoints * irreduciblePoints = nullptr;

	long getIndex(Eigen::Vector3d & pointCoords);
	friend class ActiveBandStructure;

//	double homo; // set to zero in case of phonons
//	long numValenceElectrons;
public:
	FullBandStructure(long numBands_, Statistics & statistics_,
			bool withVelocities, bool withEigenvectors,
			FullPoints & fullPoints_);
//	FullPoints * fullPoints_=nullptr,
//	IrreduciblePoints * irreduciblePoints_=nullptr);
	FullBandStructure(const FullBandStructure & that);
	FullBandStructure & operator = (const FullBandStructure & that);

	Point getPoint(const long & pointIndex);
	long getNumPoints();
	State getState(Point & point);
	void populate(PhononH0 & h0);
	void setEnergies(Eigen::Vector3d& point, Eigen::VectorXd& energies_);
	void setEnergies(Point& point, Eigen::VectorXd & energies_);
	void setEigenvectors(Point & point,
			Eigen::Tensor<std::complex<double>,3> & eigenvectors_);
	void setVelocities(Point & point,
			Eigen::Tensor<std::complex<double>,3> & velocities_);
	long getNumBands();
	bool hasIrreduciblePoints();
	Eigen::VectorXd getBandEnergies(long & bandIndex);
	Statistics getStatistics();
};

class ActiveBandStructure {
private:
	// note: we use std::vector because we want a variable number of bands
	// per each k-point
	std::vector<double> energies;
	std::vector<double> groupVelocities;
	std::vector<std::complex<double>> velocities;
	std::vector<std::complex<double>> eigenvectors;
	ActivePoints * activePoints = nullptr;
	Statistics statistics;
	bool hasEigenvectors = false;
	long numStates = 0;
	long numAtoms = 0;
	VectorXl numBands;

	// map (compressed to uncompressed)
	// inverse map (uncompressed to compressed)

	// index management
	// these are two auxiliary vectors to store indices
	MatrixXl auxBloch2Comb;
	VectorXl cumulativeKbOffset;
	VectorXl cumulativeKbbOffset;
	// this is the functionality to build the indices
	void buildIndeces(); // to be called after building the band structure
	// and these are the tools to convert indices

	// utilities to convert Bloch indices into internal indices
	long velBloch2Comb(long & ik, long & ib1, long & ib2, long & i);
	long gvelBloch2Comb(long & ik, long & ib, long & i);
	long eigBloch2Comb(long & ik, long & i, long & iat, long & ib);
	long bloch2Comb(long & k, long & b);

	long numPoints;
	bool hasPoints();
public:
	ActiveBandStructure(Statistics & statistics_);
	Statistics getStatistics();
	long getNumPoints();
	Point getPoint(const long & pointIndex);
	long getNumStates();
	State getState(Point & point);  // returns all bands at fixed k/q-point

	double getEnergy(long & stateIndex);
	Eigen::Vector3d getGroupVelocity(long & stateIndex);

	std::tuple<long,long> comb2Bloch(long & is);

	ActivePoints buildOnTheFly(Window & window, FullPoints & fullPoints,
			HarmonicHamiltonian & h0);
	ActivePoints buildAsPostprocessing(Window & window,
			FullBandStructure & fullBandStructure);
};

long compressIndeces(long i, long j, long k, long size1, long size2,
		long size3);

#endif
