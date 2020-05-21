#ifndef ACT_BANDSTRUCTURE_H
#define ACT_BANDSTRUCTURE_H

#include "points.h"
#include "bandstructure.h"
#include "state.h"
#include "window.h"
#include "statistics.h"
#include "statistics_sweep.h"
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

	template<typename T, typename S>
	static std::tuple<ActivePoints, ActiveBandStructure, S>
			builder(Context & context, T & h0, FullPoints & fullPoints);
private:
	Statistics statistics;

	// note: we don't store a matrix: we are storing an object (Nk,Nb),
	// with a variable number of bands Nb per point
	std::vector<double> energies;
	std::vector<double> groupVelocities;
	std::vector<std::complex<double>> velocities;
	std::vector<std::complex<double>> eigenvectors;

	ActivePoints * activePoints = nullptr;
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

	template<typename T>
	ActivePoints buildOnTheFly(Window & window, FullPoints & fullPoints,T &h0);

	std::vector<std::complex<double>> flattenEigenvectors(
			Eigen::MatrixXd & eigvecsIn,
			std::vector<long> bandsExtrema);

	std::vector<std::complex<double>> flattenEigenvectors(
			Eigen::Tensor<std::complex<double>,3> & eigvecsIn,
			std::vector<long> bandsExtrema);

	ActivePoints buildAsPostprocessing(Window & window,
			FullBandStructure<FullPoints> & fullBandStructure);
};

template<typename T, typename S>
std::tuple<ActivePoints, ActiveBandStructure, S>
		ActiveBandStructure::builder(Context & context,
				T & h0, FullPoints & fullPoints) {

	Statistics statistics = h0.getStatistics();

	Eigen::VectorXd temperatures = context.getTemperatures();

	ActiveBandStructure activeBandStructure(statistics);

	if ( statistics.isPhonon() ) {
		double temperatureMin = temperatures.minCoeff();
		double temperatureMax = temperatures.maxCoeff();

		Window window(context, statistics, temperatureMin, temperatureMax);

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

	// Note: eigenvectors are assumed to be phonon eigenvectors
	// might need adjustments in the future.

	numAtoms = fullPoints.getCrystal().getNumAtoms();

	std::vector<long> filteredPoints;
	std::vector<std::vector<long>> filteredBands;
	std::vector<double> filteredEnergies;
	std::vector<std::complex<double>> filteredEigenvectors;

	for ( long ik=0; ik<fullPoints.getNumPoints(); ik++ ) {
		Point point = fullPoints.getPoint(ik);

		auto [theseEnergies, theseEigenvectors] = h0.diagonalize(point);
		// eigenvectors(3,numAtoms,numBands)

		auto [ens, bandsExtrema] = window.apply(theseEnergies);

		if ( ens.empty() ) {
			continue;
		} else {

			filteredPoints.push_back(ik);
			filteredBands.push_back(bandsExtrema);

			for ( long unsigned ib=0; ib<ens.size(); ib++ ) {
				filteredEnergies.push_back(ens[ib]);
			}

			if ( h0.hasEigenvectors ) {
				auto x = flattenEigenvectors(theseEigenvectors, bandsExtrema);
				for ( auto i=0; i<x.size(); i++ ) {
					filteredEigenvectors.push_back(x);
				}
			}
		}
	}

	numPoints = filteredPoints.size();
	VectorXl filter(numPoints);
	for ( long i=0; i<numPoints; i++ ) {
		filter(i) = filteredPoints[i];
	}

	VectorXl numBands(numPoints);
	for ( long ik=0; ik<numPoints; ik++ ) {
		numBands(ik) = filteredBands[ik][1] - filteredBands[ik][0] + 1;
	}

	numStates = numBands.sum();

	ActivePoints activePoints_(fullPoints, filter);
	activePoints = &activePoints_;

	energies = filteredEnergies;

	if ( h0.hasEigenvectors ) {
		hasEigenvectors = true;
		std::vector<std::complex<double>> eigenvectors(numStates*3*numAtoms);
		eigenvectors = filteredEigenvectors;
	}

	// now we add velocities
	for ( long ik=0; ik<numPoints; ik++ ) {
		Point point = activePoints_.getPoint(ik);

		// thisVelocity is a tensor of dimensions (ib, ib, 3)
		auto thisVelocity = h0.diagonalizeVelocity(point);

		// now we filter it
		for ( long ib1Old = filteredBands[ik][0];
				ib1Old<filteredBands[ik][1]+1; ib1Old++ ) {
			for ( long i=0; i<3; i++ ) {
				groupVelocities.push_back(thisVelocity(ib1Old,ib1Old,i).real());
				for ( long ib2Old = filteredBands[ik][0];
						ib2Old<filteredBands[ik][1]+1; ib2Old++ ) {
					velocities.push_back(thisVelocity(ib1Old,ib2Old,i));
				}
			}
		}
	}
	return activePoints_;
}

#endif
