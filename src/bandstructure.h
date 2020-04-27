#ifndef BANDSTRUCTURE_H
#define BANDSTRUCTURE_H

#include "points.h"
#include "phononH0.h"
#include "state.h"
#include "window.h"
#include "statistics.h"

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
//	Eigen::MatrixXd dndt;
//	Eigen::MatrixXd dnde;
//	void setOccupations();
	bool useIrreducible = false;

	bool hasEigenvectors = false;
	bool hasVelocities = false;

	FullPoints * fullPoints = nullptr;
	IrreduciblePoints * irreduciblePoints = nullptr;
	long getIndex(Eigen::Vector3d & pointCoords);

	friend class ActiveBandStructure;

//	double homo; // set to zero in case of phonons
//	long numValenceElectrons;
public:
	FullBandStructure(long numBands_, Statistics & statistics_,
			bool withVelocities, bool withEigenvectors,
			FullPoints* fullPoints_=nullptr,
			IrreduciblePoints* irreduciblePoints_=nullptr);
	Point getPoint(const long & pointIndex);
	long getNumPoints();
//	State getStateFromPointIndex(const long & index);
	State getState(Point & point);
//	void setChemicalPotential(double chemPot);
//	void setTemperature(double temp);
	void populate(HarmonicHamiltonian & h0);
	void setEnergies(Eigen::Vector3d& point, Eigen::VectorXd& energies_);
//	void setVelocities(Eigen::Vector3d& pointCoords,
//			Eigen::Tensor<std::complex<double>,3>& velocities_);
	void setEnergies(Point& point, Eigen::VectorXd & energies_);
	void setEigenvectors(Point & point,
			Eigen::Tensor<std::complex<double>,3> & eigenvectors_);
	void setVelocities(Point & point,
			Eigen::Tensor<std::complex<double>,3> & velocities_);
	long getNumBands();
	bool hasIrreduciblePoints();
//	void setNumValenceElectrons(long numElectrons);
//	void setHomo(double homo);

	Eigen::VectorXd getBandEnergies(long & bandIndex);
};

class ActiveBandStructure {
private:
	// note: we use std::vector because we want a variable number of bands
	// per each k-point
	std::vector<std::vector<double>> energies;
	std::vector<std::vector<std::complex<double>>> velocities;
	std::vector<std::vector<std::complex<double>>> eigenvectors;
//	Eigen::VectorXd * dndt = nullptr;
//	Eigen::VectorXd * dnde = nullptr;
	ActivePoints * activePoints = nullptr;
	Statistics statistics;
	bool hasEigenvectors = false;
	long numStates = 0;
	long numAtoms = 0;
	std::vector<std::vector<long>> filteredBands;
	long numPoints;
	bool hasPoints();
public:
	ActiveBandStructure(Statistics & statistics_);
	long getNumPoints();
	Point getPoint(const long & pointIndex);
	long getNumStates();
	State getState(Point & point);  // returns all bands at fixed k/q-point

	ActivePoints buildOnTheFly(Window & window, FullPoints & fullPoints,
			HarmonicHamiltonian & h0);
	ActivePoints buildAsPostprocessing(Window & window,
			FullBandStructure & fullBandStructure);
};

long compressIndeces(long i, long j, long k, long size1, long size2,
		long size3);

#endif
