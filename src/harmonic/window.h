#ifndef WINDOW_H
#define WINDOW_H

#include "context.h"
#include "eigen.h"
#include "statistics.h"

class Window {
public:
	Window(Context & context, Statistics & statistics_);
	std::tuple<std::vector<double>,std::vector<long>> apply(
			Eigen::VectorXd & energies);
	bool canOnTheFly;

	const int nothing = 0;
	const int population = 1;
	const int energy = 2;
private:
	int method;
	double populationThreshold = 0.;
	double minEnergy = 0., maxEnergy = 0.;
	long numBands;

	double chemicalPotentialMin, chemicalPotentialMax;
	double maxTemperature;

	Statistics & statistics;

	std::tuple<std::vector<double>,std::vector<long>> internalPopWindow(
			Eigen::VectorXd& energies,
			Eigen::VectorXd& dndtMin, Eigen::VectorXd& dndtMax,
			Eigen::VectorXd& dndeMin, Eigen::VectorXd& dndeMax);
	std::tuple<std::vector<double>,std::vector<long>> internalEnWindow(
			Eigen::VectorXd energies);
};

#endif
