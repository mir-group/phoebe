#include "context.h"
#include "eigen.h"
#include "window.h"
#include "particle.h"
#include "exceptions.h"

Window::Window(Context &context, Particle &particle_,
        const double &temperatureMin_, const double &temperatureMax_,
        const double &chemicalPotentialMin_,
        const double &chemicalPotentialMax_) :
        particle(particle_), temperatureMin(temperatureMin_),
        temperatureMax(temperatureMax_),
        chemicalPotentialMin(chemicalPotentialMin_),
        chemicalPotentialMax(chemicalPotentialMax_) {

    std::string inMethod = context.getWindowType();
    if (inMethod == "population") {
        method = population;
    } else if (inMethod == "energy") {
        method = energy;
    } else if (inMethod == "nothing") {
        method = nothing;
    } else {
        Error e("Unrecognized method called in Window()", 1);
    }

    if (method == population) {
        populationThreshold = context.getWindowPopulationLimit();

        if (std::isnan(populationThreshold)) {
            populationThreshold = 0.01;
        }

    } else if (method == energy) {
        minEnergy = context.getWindowEnergyLimit().minCoeff();
        maxEnergy = context.getWindowEnergyLimit().maxCoeff();
    }
}

std::tuple<std::vector<double>, std::vector<long>> Window::apply(
        Eigen::VectorXd &energies) {

    if (method == population) {
        numBands = energies.size();
        Eigen::VectorXd dndtMin(numBands), dndtMax(numBands), dndeMin(numBands),
                dndeMax(numBands);
        for (long ib = 0; ib < numBands; ib++) {
            dndtMin(ib) = particle.getDndt(energies(ib), temperatureMax,
                    chemicalPotentialMin);
            dndtMax(ib) = particle.getDndt(energies(ib), temperatureMax,
                    chemicalPotentialMin);
            dndeMin(ib) = particle.getDnde(energies(ib), temperatureMax,
                    chemicalPotentialMin);
            dndeMax(ib) = particle.getDnde(energies(ib), temperatureMax,
                    chemicalPotentialMin);
        }
        auto tup = internalPopWindow(energies,                dndtMin, dndtMax, dndeMin, dndeMax);
        auto filteredEnergies = std::get<0>(tup);
        auto bandExtrema = std::get<1>(tup);
        return {filteredEnergies, bandExtrema};
    } else if (method == energy) {
      auto tup = internalEnWindow(energies);
      auto filteredEnergies = std::get<0>(tup);
      auto bandExtrema = std::get<1>(tup);
        return {filteredEnergies, bandExtrema};
    } else { // no filter
        numBands = energies.size();
        std::vector<double> filteredEnergies(energies.size());
        for (long ib = 0; ib < numBands; ib++) {
            filteredEnergies[ib] = ib;
        }
        std::vector<long> bandExtrema(2);
        bandExtrema[0] = 0;
        bandExtrema[1] = numBands - 1;
        return {filteredEnergies, bandExtrema};
    }
}

std::tuple<std::vector<double>, std::vector<long>> Window::internalPopWindow(
        Eigen::VectorXd &energies, Eigen::VectorXd &dndtMin,
        Eigen::VectorXd &dndtMax, Eigen::VectorXd &dndeMin,
        Eigen::VectorXd &dndeMax) {

    std::vector<double> filteredEnergies;
    std::vector<long> bandsExtrema;
    std::vector<long> bandsIndeces(energies.size());

    double thisEnergy;
    for (long ib = 0; ib < numBands; ib++) {
        thisEnergy = energies(ib);
        if ((dndtMin(ib) > populationThreshold)
                || (dndtMax(ib) > populationThreshold)
                || (dndeMin(ib) > populationThreshold)
                || (dndeMax(ib) > populationThreshold)) {
            filteredEnergies.push_back(thisEnergy);
            bandsIndeces.push_back(ib);
        }
    }
    if (bandsIndeces.size() > 0) {
        bandsExtrema.push_back(bandsIndeces[0]);
        bandsExtrema.push_back(bandsIndeces[bandsIndeces.size() - 1]);
    }
    return {filteredEnergies,bandsExtrema};
}

std::tuple<std::vector<double>, std::vector<long>> Window::internalEnWindow(
        Eigen::VectorXd energies) {

    std::vector<double> filteredEnergies;
    std::vector<long> bandsExtrema(2);
    std::vector<long> bandsIndeces(energies.size());

    double thisEnergy;
    for (long ib = 0; ib < numBands; ib++) {
        thisEnergy = energies(ib);
        if (thisEnergy < maxEnergy && thisEnergy > minEnergy) {
            filteredEnergies.push_back(thisEnergy);
            bandsIndeces.push_back(ib);
        }
    }
    bandsExtrema.push_back(bandsIndeces[0]);
    bandsExtrema.back();
    return {filteredEnergies,bandsExtrema};
}

long Window::getMethodUsed() {
    return method;
}
