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
			FullBandStructure<FullPoints> & fullBandStructure);
	// copy constructor
	StatisticsSweep(const StatisticsSweep & that);
	// copy assignment
	StatisticsSweep & operator = (const StatisticsSweep & that);

	CalcStatistics getCalcStatistics(const long & index);
	struct CalcStatistics getCalcStatistics(const long & iTemp,
			const long & iChemPot);
	long getNumCalcs();
private:
	// note: these three private methods are only meant to be used
	// during the construction of this class, as they depend on energies
	double findChemicalPotentialFromDoping(const double & doping,
			const double & temperature);
	double findDopingFromChemicalPotential(const double & chemicalPotential,
			const double & temperature);
	double fPop(const double & chemPot, const double & temp);

	const long maxIter = 100;

	// temp variables, available only during the default constructor
	long numStates;
	Eigen::VectorXd energies;
	long numPoints;
	double numElectronsDoped;

    Statistics statistics;

	double volume;
	double spinFactor;
	double occupiedStates;
	double fermiLevel;

	// the relevant quantities to be outputted:
	long numCalcs;
	Eigen::MatrixXd infoCalcs;
	long nTemp;
	long nChemPot;
	long nDop;
};
