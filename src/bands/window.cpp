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
      populationThreshold = 1.0e-10;
      // note: we are filtering on the factor N(N+1) or F(1-F)
      // tested this value to be a safe one.
    }

  } else if (method == energy) {
    minEnergy = context.getWindowEnergyLimit().minCoeff();
    maxEnergy = context.getWindowEnergyLimit().maxCoeff();
    if (std::isnan(minEnergy) || std::isnan(maxEnergy) ) {
      Error e("You must set min and max energies for your energy window!"); 
    }
    else if(minEnergy == maxEnergy) {
      Error e("Your min and max window energies cannot be the same."); 
    }
  }
}

std::tuple<std::vector<double>, std::vector<int>> Window::apply(
    Eigen::VectorXd &energies) {

  // no matter the method, we must set numBands
  numBands = energies.size();

  if (method == population) {
    Eigen::VectorXd popMin(numBands), popMax(numBands);
    for (long ib = 0; ib < numBands; ib++) {
      popMin(ib) = particle.getPopPopPm1(energies(ib), temperatureMax,
                                         chemicalPotentialMin);
      popMax(ib) = particle.getPopPopPm1(energies(ib), temperatureMax,
                                         chemicalPotentialMax);
    }
    auto tup = internalPopWindow(energies, popMin, popMax);
    auto filteredEnergies = std::get<0>(tup);
    auto bandExtrema = std::get<1>(tup);
    return {filteredEnergies, bandExtrema};

  } else if (method == energy) {
    auto tup = internalEnWindow(energies);
    auto filteredEnergies = std::get<0>(tup);
    auto bandExtrema = std::get<1>(tup);
    return {filteredEnergies, bandExtrema};

  } else { // no filter
    std::vector<double> filteredEnergies(energies.size());
    for (long ib = 0; ib < numBands; ib++) {
      filteredEnergies[ib] = ib;
    }
    std::vector<int> bandExtrema(2);
    bandExtrema[0] = 0;
    bandExtrema[1] = numBands - 1;
    return {filteredEnergies, bandExtrema};
  }
}

std::tuple<std::vector<double>, std::vector<int>> Window::internalPopWindow(
    const Eigen::VectorXd &energies, const Eigen::VectorXd &popMin,
    const Eigen::VectorXd &popMax) {

  std::vector<double> filteredEnergies;
  std::vector<int> bandsExtrema;
  std::vector<int> bandsIndeces;

  double thisEnergy;
  for (int ib = 0; ib < numBands; ib++) {
    thisEnergy = energies(ib);
    if ((abs(popMin(ib)) > populationThreshold) ||
        (abs(popMax(ib)) > populationThreshold)) {
      filteredEnergies.push_back(thisEnergy);
      bandsIndeces.push_back(ib);
    }
  }
  if (bandsIndeces.size() > 0) {
    bandsExtrema.push_back(bandsIndeces[0]);
    bandsExtrema.push_back(bandsIndeces[bandsIndeces.size() - 1]);
  } // or return empty lists if nothing is found
  return {filteredEnergies, bandsExtrema};
}

std::tuple<std::vector<double>, std::vector<int>> Window::internalEnWindow(
    const Eigen::VectorXd &energies) {

  std::vector<double> filteredEnergies;
  std::vector<int> bandsExtrema(2);
  std::vector<int> bandsIndeces;
  double thisEnergy;
  for (int ib = 0; ib < numBands; ib++) {
    thisEnergy = energies(ib);
    if (thisEnergy < maxEnergy && thisEnergy > minEnergy) {
      filteredEnergies.push_back(thisEnergy);
      bandsIndeces.push_back(ib);
    }
  }
  // set the band extrema
  if (bandsIndeces.size() > 0) {
    bandsExtrema.push_back(bandsIndeces[0]);
    bandsExtrema.push_back(bandsIndeces[bandsIndeces.size() - 1]);
    bandsExtrema.back();
  } // or return empty lists if nothing is found
  return {filteredEnergies, bandsExtrema};
}

long Window::getMethodUsed() {
  return method;
}
