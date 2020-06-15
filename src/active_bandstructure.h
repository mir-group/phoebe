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
	ActiveBandStructure(Particle & particle_);
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

	State getState(Point & point);  // returns all bands at fixed k/q-point

	long getIndex(const WavevectorIndex & ik, const BandIndex & ib);
	long getNumStates();

	const double & getEnergy(const long & stateIndex);
	Eigen::Vector3d getGroupVelocity(const long & stateIndex);
	Eigen::Vector3d getWavevector(const long & stateIndex);

	void setEnergies(Point & point, Eigen::VectorXd & energies_);
	void setEigenvectors(Point & point,
			Eigen::MatrixXcd & eigenvectors_);
	void setVelocities(Point & point,
			Eigen::Tensor<std::complex<double>,3> & velocities_);

	template<typename T>
	static std::tuple<ActivePoints, ActiveBandStructure>
			builder(Context & context, T & h0, FullPoints & fullPoints);
protected:
	Particle particle;

	// note: we don't store a matrix: we are storing an object (Nk,Nb),
	// with a variable number of bands Nb per point
	std::vector<double> energies;
	std::vector<std::complex<double>> velocities;
	std::vector<std::complex<double>> eigenvectors;

	ActivePoints * activePoints = nullptr;
	bool hasEigenvectors = false;
	long numStates = 0;
	long numAtoms = 0;
	long numPoints;
	bool hasPoints();

	VectorXl numBands;
	long numFullBands;

	// index management
	// these are two auxiliary vectors to store indices
	MatrixXl auxBloch2Comb;
	VectorXl cumulativeKbOffset;
	VectorXl cumulativeKbbOffset;
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


	template<typename T>
	ActivePoints buildOnTheFly(Window & window, FullPoints & fullPoints,T &h0);

	ActivePoints buildAsPostprocessing(Window & window,
			FullBandStructure & fullBandStructure);
};

template<typename T>
std::tuple<ActivePoints, ActiveBandStructure> ActiveBandStructure::builder(
		Context & context, T & h0, FullPoints & fullPoints) {

	Particle particle = h0.getParticle();

	Eigen::VectorXd temperatures = context.getTemperatures();

	ActiveBandStructure activeBandStructure(particle);

	if ( particle.isPhonon() ) {
		double temperatureMin = temperatures.minCoeff();
		double temperatureMax = temperatures.maxCoeff();

		Window window(context, particle, temperatureMin, temperatureMax);

		auto aPoints = activeBandStructure.buildOnTheFly(window,fullPoints,h0);
		StatisticsSweep statisticsSweep(context);
		return {aPoints, activeBandStructure, statisticsSweep};
	} else {
		Error e("apply window for electrons not implemented");
	}
}


template<typename T>
ActivePoints ActiveBandStructure::buildOnTheFly(Window & window,
		FullPoints & fullPoints, T & h0) {
	// this function proceeds in three logical blocks:
	// 1- we find out the list of "relevant" points
	// 2- initialize internal raw buffer for energies, velocities, eigvecs
	// 3- populate the raw buffer

	// Note: eigenvectors are assumed to be phonon eigenvectors
	// might need adjustments in the future.

	numAtoms = fullPoints.getCrystal().getNumAtoms();

	// we have to build this in a way that works in parallel
	// ALGORITHM:
	// - loop over points. Diagonalize, and find if we want this k-point
	//   (while we are at it, we could save energies and the eigenvalues)
	// - find how many points each MPI rank has found
	//   MPI_GATHER
	// - use MPI_GATHERV to receive the indices of points to be considered
	// - either recollect energies, or loop over these points and compute them

	numFullBands = 0; // save the unfiltered number of bands
	std::vector<long> filteredPoints;
	std::vector<std::vector<long>> filteredBands;
	for ( long ik=0; ik<fullPoints.getNumPoints(); ik++ ) {
		Point point = fullPoints.getPoint(ik);
		// diagonalize harmonic hamiltonian
		auto [theseEnergies, theseEigenvectors] = h0.diagonalize(point);
		// ens is empty if no "relevant" energy is found.
		// bandsExtrema contains the lower and upper band index of "relevant"
		// bands at this point
		auto [ens, bandsExtrema] = window.apply(theseEnergies);
		if ( ens.empty() ) { // nothing to do
			continue;
		} else { // save point index and "relevant" band indices
			filteredPoints.push_back(ik);
			filteredBands.push_back(bandsExtrema);
		}
		numFullBands = theseEnergies.size();
	}

	// Here, we should use MPI_gather(v) to put filteredPoints in common

	numPoints = filteredPoints.size();
	// this vector mapps the indices of the new point to the old list
	VectorXl filter(numPoints);
	for ( long i=0; i<numPoints; i++ ) {
		filter(i) = filteredPoints[i];
	}

	// numBands is a book-keeping of how many bands per kpoint there are
	// this isn't a constant number.
	// on top of that, we look for the size of the arrays containing bandstruc.
	VectorXl numBands(numPoints);
	long numEnStates = 0;
	long numVelStates = 0;
	long numEigStates = 0;
	for ( long ik=0; ik<numPoints; ik++ ) {
		numBands(ik) = filteredBands[ik][1] - filteredBands[ik][0] + 1;
		//
		numEnStates += numBands(ik);
		numVelStates += 3 * numBands(ik) * numBands(ik);
		numEigStates += numBands(ik) * numFullBands;
	}

	// initialize the raw data buffers of the activeBandStructure

	ActivePoints activePoints_(fullPoints, filter);
	activePoints = &activePoints_;
	// construct the mapping from combined indices to Bloch indices
	buildIndeces();

	energies.resize(numEnStates,0.);
	velocities.resize(numVelStates,complexZero);
	if ( h0.hasEigenvectors ) {
		hasEigenvectors = true;
		eigenvectors.resize(numEigStates,complexZero);
	}

	/////////////////

	// now we can loop over the trimmed list of points
	for ( long ik=0; ik<numPoints; ik++ ) {
		Point point = activePoints->getPoint(ik);
		auto [theseEnergies, theseEigenvectors] = h0.diagonalize(point);
		// eigenvectors(3,numAtoms,numBands)
		auto [ens, bandsExtrema] = window.apply(theseEnergies);

		setEnergies(ik, ens);

		if ( h0.hasEigenvectors ) {
			// we are reducing the basis size!
			// the first index has the size of the Hamiltonian
			// the second index has the size of the filtered bands
			Eigen::MatrixXcd theseEigvecs(numBands,numBands(ik));
			long ibAct = 0;
			for ( long ibFull=filteredBands[ik][0];
					ibFull<=filteredBands[ik][1]; ibFull++) {
				theseEigvecs.col(ibAct) = theseEigenvectors.col(ibFull);
				ibAct++;
			}
			setEigenvectors(point, theseEigvecs);
		}

		// thisVelocity is a tensor of dimensions (ib, ib, 3)
		auto thisVelocity = h0.diagonalizeVelocity(point);

		// now we filter it
		Eigen::Tensor<std::complex<double>,3> thisVels(numBands(ik),
				numBands(ik),3);
		long ib1New = 0;
		for ( long ib1Old = filteredBands[ik][0];
				ib1Old<filteredBands[ik][1]+1; ib1Old++ ) {
			long ib2New = 0;
			for ( long ib2Old = filteredBands[ik][0];
					ib2Old<filteredBands[ik][1]+1; ib2Old++ ) {
				for ( long i=0; i<3; i++ ) {
					thisVels(ib1New,ib2New,i) = thisVelocity(ib1Old,ib2Old,i);
				}
				ib2New++;
			}
			ib1New++;
		}
		setVelocities(point, thisVels);
	}

	return activePoints_;
}

#endif
