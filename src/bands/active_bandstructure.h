#ifndef ACT_BANDSTRUCTURE_H
#define ACT_BANDSTRUCTURE_H

#include "points.h"
#include "bandstructure.h"
#include "state.h"
#include "window.h"
#include "particle.h"
#include "statistics_sweep.h"
#include "harmonic.h"
#include "constants.h"

class ActiveBandStructure : public BaseBandStructure {
public:
	ActiveBandStructure(Particle & particle_, ActivePoints & activePoints);
	ActiveBandStructure(const ActiveBandStructure & that);// copy
	ActiveBandStructure & operator=(const ActiveBandStructure & that); //assign

	Particle getParticle();

	Points getPoints();
	// return a point from its index. pointIndex goes over the list of points
	// in the ActivePoints object (from 0 to numActivePoints)
	Point getPoint(const long & pointIndex);
	// returns the number of active points
	long getNumPoints();
	long getNumBands(); // this only works in FullBandStructure
	long hasWindow();

	State getState(const long & pointIndex);
	State getState(Point & point);  // returns all bands at fixed k/q-point

	long getIndex(const WavevectorIndex & ik, const BandIndex & ib);
	long getNumStates();

	const double & getEnergy(const long & stateIndex);
	Eigen::Vector3d getGroupVelocity(const long & stateIndex);
	Eigen::Vector3d getWavevector(const long & stateIndex);

	void setEnergies(Point & point, Eigen::VectorXd & energies_);
	void setEnergies(Point & point, std::vector<double> & energies_);
	void setEigenvectors(Point & point,
			Eigen::MatrixXcd & eigenvectors_);
	void setVelocities(Point & point,
			Eigen::Tensor<std::complex<double>,3> & velocities_);

	static std::tuple<ActiveBandStructure, StatisticsSweep>
			builder(Context & context,
					HarmonicHamiltonian & h0,
					FullPoints & fullPoints,
					const bool & withEigenvectors=true,
					const bool & withVelocities=true);
protected:
	// stores the quasiparticle kind
	Particle particle;
	ActivePoints activePoints;

	// note: we don't store a matrix: we are storing an object (Nk,Nb),
	// with a variable number of bands Nb per point
	std::vector<double> energies;
	std::vector<std::complex<double>> velocities;
	std::vector<std::complex<double>> eigenvectors;

	bool hasEigenvectors = false;
	long numStates = 0;
	long numPoints;
	bool hasPoints();

	Eigen::VectorXi numBands;
	long numFullBands;
	long windowMethod;

	// index management
	// these are two auxiliary vectors to store indices
	Eigen::MatrixXi auxBloch2Comb;
	Eigen::VectorXi cumulativeKbOffset;
	Eigen::VectorXi cumulativeKbbOffset;
	// this is the functionality to build the indices
	void buildIndeces(); // to be called after building the band structure
	// and these are the tools to convert indices

	// utilities to convert Bloch indices into internal indices
	long velBloch2Comb(const long & ik, const long & ib1, const long & ib2,
			const long & i);
	long eigBloch2Comb(const long & ik, const long & ibFull,
			const long & ibRed);
	long bloch2Comb(const long & k, const long & b);
	std::tuple<long,long> comb2Bloch(const long & is);

	void buildOnTheFly(Window & window, FullPoints & fullPoints,
			HarmonicHamiltonian & h0, const bool & withEigenvectors=true,
			const bool & withVelocities=true);

//	ActivePoints buildAsPostprocessing(Window & window,
//			FullBandStructure & fullBandStructure);
};

#endif
