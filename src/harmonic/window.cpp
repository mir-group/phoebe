#include "context.h"
#include "eigen.h"
#include "window.h"
#include "statistics.h"
#include "exceptions.h"

Window::Window(Context & context, Statistics & statistics_) :
	statistics{statistics_} {

	std::string inMethod = context.getWindowType();

	if ( inMethod != "population" && inMethod != "energy" &&
			inMethod != "nothing" ){
		Error e("Unrecognized method called in Window()", 1);
	}

	if ( inMethod == "population" ) {
		method = population;
	} else if ( inMethod == "energy" ) {
		method = energy;
	} else {
		method = nothing;
	}

	if ( method == population ) {
		canOnTheFly = false;
		populationThreshold = context.getWindowPopulationLimit();
		if ( statistics.isFermi() ) {

			FIX HERE: should get chemical potentials from StatisticSweep

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

std::tuple<std::vector<double>,std::vector<long>> Window::apply(
		Eigen::VectorXd & energies) {

	if ( method == population ) {
		numBands = energies.size();
		Eigen::VectorXd dndtMin(numBands), dndtMax(numBands),
				dndeMin(numBands), dndeMax(numBands);
		for ( long ib=0; ib<numBands; ib++ ) {
			dndtMin(ib) = statistics.getDndt(energies(ib), maxTemperature,
					chemicalPotentialMin);
			dndtMax(ib) = statistics.getDndt(energies(ib), maxTemperature,
					chemicalPotentialMin);
			dndeMin(ib) = statistics.getDnde(energies(ib), maxTemperature,
					chemicalPotentialMin);
			dndeMax(ib) = statistics.getDnde(energies(ib), maxTemperature,
					chemicalPotentialMin);
		}
		auto [filteredEnergies, bandExtrema] = internalPopWindow(energies,
				dndtMin, dndtMax, dndeMin, dndeMax);
		return {filteredEnergies, bandExtrema};
	} else if ( method == energy ) {
		auto [filteredEnergies, bandExtrema] = internalEnWindow(energies);
		return {filteredEnergies, bandExtrema};
	} else { // no filter
		numBands = energies.size();
		std::vector<double> filteredEnergies(energies.size());
		for ( long ib=0; ib<numBands; ib++) {
			filteredEnergies[ib] = ib;
		}
		std::vector<long> bandExtrema(2);
		bandExtrema[0] = 0;
		bandExtrema[1] = numBands-1;
		return {filteredEnergies, bandExtrema};
	}
}

std::tuple<std::vector<double>,std::vector<long>> Window::internalPopWindow(
		Eigen::VectorXd& energies,
		Eigen::VectorXd& dndtMin, Eigen::VectorXd& dndtMax,
		Eigen::VectorXd& dndeMin, Eigen::VectorXd& dndeMax) {

	std::vector<double> filteredEnergies;
	std::vector<long> bandsExtrema(2);
	std::vector<long> bandsIndeces(energies.size());

	double thisEnergy;
	for ( long ib=0; ib<numBands; ib++ ) {
		thisEnergy = energies(ib);
		if ( ( dndtMin(ib) > populationThreshold ) ||
				( dndtMax(ib) > populationThreshold ) ||
				( dndeMin(ib) > populationThreshold ) ||
				( dndeMax(ib) > populationThreshold ) ) {
			filteredEnergies.push_back(thisEnergy);
			bandsIndeces.push_back(ib);
		}
	}
	bandsExtrema.push_back(bandsIndeces[0]);
	bandsExtrema.back();
	return {filteredEnergies,bandsExtrema};
}

std::tuple<std::vector<double>,std::vector<long>> Window::internalEnWindow(
		Eigen::VectorXd energies) {

	std::vector<double> filteredEnergies;
	std::vector<long> bandsExtrema(2);
	std::vector<long> bandsIndeces(energies.size());

	double thisEnergy;
	for ( long ib=0; ib<numBands; ib++ ) {
		thisEnergy = energies(ib);
		if ( thisEnergy < maxEnergy && thisEnergy > minEnergy ) {
			filteredEnergies.push_back(thisEnergy);
			bandsIndeces.push_back(ib);
		}
	}
	bandsExtrema.push_back(bandsIndeces[0]);
	bandsExtrema.back();
	return {filteredEnergies,bandsExtrema};
}



