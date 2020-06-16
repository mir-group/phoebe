#ifndef WINDOW_H
#define WINDOW_H

#include "context.h"
#include "eigen.h"
#include "particle.h"

class Window {
public:
	Window(Context & context, Particle & particle_,
			const double & temperatureMin = 0.,
			const double & temperatureMax = 0.,
			const double & chemicalPotentialMin = 0.,
			const double & chemicalPotentialMax = 0.);

	// public interface to use the window, used by activeBandStructure
	std::tuple<std::vector<double>,std::vector<long>> apply(
			Eigen::VectorXd & energies);

	static const long nothing = 0;
	static const long population = 1;
	static const long energy = 2;

	long getMethodUsed();
private:
	Particle & particle;

	// parameters for window
	double temperatureMin, temperatureMax;
	double chemicalPotentialMin, chemicalPotentialMax;
	double populationThreshold = 0.;
	double minEnergy = 0., maxEnergy = 0.;

	// selection of window type
	int method;

	// temp variable
	long numBands;


	// internal method to apply the window
	std::tuple<std::vector<double>,std::vector<long>> internalPopWindow(
			Eigen::VectorXd& energies,
			Eigen::VectorXd& dndtMin, Eigen::VectorXd& dndtMax,
			Eigen::VectorXd& dndeMin, Eigen::VectorXd& dndeMax);
	std::tuple<std::vector<double>,std::vector<long>> internalEnWindow(
			Eigen::VectorXd energies);
};

#endif
