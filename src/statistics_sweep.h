#ifndef SWEEP_H
#define SWEEP_H

#include <memory>
#include "context.h"
#include "bandstructure.h"
#include "statistics.h"

struct CalcStatistics {
	double temperature;
	double chemicalPotential;
	double doping;
};

class StatisticsSweep {
public:
	StatisticsSweep(Context & context,
			FullBandStructure<FullPoints> * fullBandStructure=nullptr);
	// copy constructor
	StatisticsSweep(const StatisticsSweep & that);
	// copy assignment
	StatisticsSweep & operator = (const StatisticsSweep & that);

	struct CalcStatistics getCalcStatistics(const long & index);
	struct CalcStatistics getCalcStatistics(const long & iTemp,
			const long & iChemPot);
	long getNumCalcs();
	long getNumChemicalPotentials();
	long getNumTemperatures();

protected:
    Statistics statistics;
	long numCalcs;
	Eigen::MatrixXd infoCalcs;
	long nTemp;
	long nChemPot;
	long nDop;

	// note: these three private methods are only meant to be used
	// during the construction of this class, as they depend on energies
	double findChemicalPotentialFromDoping(const double & doping,
			const double & temperature);
	double findDopingFromChemicalPotential(const double & chemicalPotential,
			const double & temperature);
	double fPop(const double & chemPot, const double & temp);

	// this block is temporary variables for electronic calculations
	const long maxIter = 100;
	long numPoints;
	long numBands;
	long numStates;
	Eigen::VectorXd energies;
	double numElectronsDoped;
	double volume;
	double spinFactor;
	double occupiedStates;
	double fermiLevel;

};

#endif
