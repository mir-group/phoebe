#ifndef ACT_BANDSTRUCTURE_H
#define ACT_BANDSTRUCTURE_H

#include "points.h"
#include "bandstructure.h"
#include "state.h"
#include "window.h"
#include "statistics.h"
#include "harmonic.h"

class ActiveBandStructure {
public:
	ActiveBandStructure(Statistics & statistics_);
	ActiveBandStructure(const ActiveBandStructure & that);// copy
	ActiveBandStructure & operator=(const ActiveBandStructure & that); // assign

	Statistics getStatistics();
	long getNumPoints();
	long getNumStates();

	Point<ActivePoints> getPoint(const long & pointIndex);

	State<ActivePoints> getState(Point<ActivePoints> & point);  // returns all bands at fixed k/q-point

	double getEnergy(long & stateIndex);
	Eigen::Vector3d getGroupVelocity(long & stateIndex);

	static std::tuple<ActivePoints, ActiveBandStructure, StatisticsSweep>
			builder(Context & context, HarmonicHamiltonian & h0,
					FullPoints & fullPoints);
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
	std::tuple<long,long> comb2Bloch(long & is);

	long numPoints;
	bool hasPoints();

	ActivePoints buildOnTheFly(FullPoints & fullPoints,
			HarmonicHamiltonian & h0);

	ActivePoints buildAsPostprocessing(Window & window,
			FullBandStructure<FullPoints> & fullBandStructure);
};

#endif
