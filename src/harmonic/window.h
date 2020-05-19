#ifndef WINDOW_H
#define WINDOW_H

#include "context.h"
#include "eigen.h"
#include "statistics.h"

template<typename T>
class Window {
public:
	Window(Context & context, Statistics & statistics_);
	std::tuple<std::vector<double>,std::vector<long>> apply(
			Eigen::VectorXd & energies);
	bool canOnTheFly;

	static const long nothing = 0;
	static const long population = 1;
	static const long energy = 2;
private:
	T & statisticsSweep;
	int method;

	double populationThreshold = 0.;
	double minEnergy = 0., maxEnergy = 0.;
	long numBands;
	double chemicalPotentialMin, chemicalPotentialMax;
	double maxTemperature;

	std::tuple<std::vector<double>,std::vector<long>> internalPopWindow(
			Eigen::VectorXd& energies,
			Eigen::VectorXd& dndtMin, Eigen::VectorXd& dndtMax,
			Eigen::VectorXd& dndeMin, Eigen::VectorXd& dndeMax);
	std::tuple<std::vector<double>,std::vector<long>> internalEnWindow(
			Eigen::VectorXd energies);
};

template<>
Window::Window(const long & method_, T & statisticsSweep_) :
	statisticSweep(statisticsSweep_), method(method_) {

	if ( method != nothing || method != population || method != energy ) {
		Error e("Unrecognized method called in Window()", 1);
	}

	if ( method == population ) {
		canOnTheFly = false;
		populationThreshold = context.getWindowPopulationLimit();
		if ( statisticSweep.getStatistics().isFermi() ) {
			chemicalPotentialMin = context.getChemicalPotentials().minCoeff();
			chemicalPotentialMax = context.getChemicalPotentials().maxCoeff();
		} else {
			chemicalPotentialMin = 0.;
			chemicalPotentialMax = 0.;
		}
		maxTemperature = context.getTemperatures().maxCoeff();
	} else {
		canOnTheFly = true;
		minEnergy = context.getWindowEnergyLimit().minCoeff();
		maxEnergy = context.getWindowEnergyLimit().maxCoeff();
	}
}

#endif
