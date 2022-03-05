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
  } else if (inMethod == "magnetotransport") {
    method = forMagnetotransport;
  } else if (inMethod == "nothing") {
    method = nothing;
  } else {
    Error("Unrecognized method called in Window()");
  }

  if (method == population || method == forMagnetotransport) {
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
      Error("You must set min and max energies for your energy window!");
    }
    else if(minEnergy == maxEnergy) {
      Error("Your min and max window energies cannot be the same.");
    }
  }
}

std::tuple<std::vector<double>, std::vector<int>> Window::apply(
    Eigen::VectorXd &energies) {

  // no matter the method, we must set numBands
  numBands = int(energies.size());

  if (method == population) {
    Eigen::VectorXd popMin(numBands), popMax(numBands);
    for (int ib = 0; ib < numBands; ib++) {
      popMin(ib) = particle.getPopPopPm1(energies(ib), temperatureMax,
                                         chemicalPotentialMin);
      popMax(ib) = particle.getPopPopPm1(energies(ib), temperatureMax,
                                         chemicalPotentialMax);
    }
    return internalPopWindow(energies, popMin, popMax);
  }
  else if (method == forMagnetotransport) {
    Eigen::VectorXd popMin(numBands), popMax(numBands);
    for (int ib = 0; ib < numBands; ib++) {
      popMin(ib) = particle.getPopPopPm1(energies(ib), temperatureMax,
                                         chemicalPotentialMin);
      popMax(ib) = particle.getPopPopPm1(energies(ib), temperatureMax,
                                         chemicalPotentialMax);
    }
    return internalMagWindow(energies, popMin, popMax);
  } else if (method == energy) {
    return internalEnWindow(energies);
  } else { // no filter
    std::vector<double> filteredEnergies(energies.size());
    for (int ib = 0; ib < numBands; ib++) {
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
    const Eigen::VectorXd &popMax) const {

  std::vector<double> filteredEnergies;
  std::vector<int> bandsExtrema;
  std::vector<int> bandsIndices;

  double thisEnergy;
  for (int ib = 0; ib < numBands; ib++) {
    thisEnergy = energies(ib);
    if ((abs(popMin(ib)) > populationThreshold) ||
        (abs(popMax(ib)) > populationThreshold)) {
      filteredEnergies.push_back(thisEnergy);
      bandsIndices.push_back(ib);
    }
  }
  if (!bandsIndices.empty()) {
    bandsExtrema.push_back(bandsIndices[0]);
    bandsExtrema.push_back(bandsIndices[bandsIndices.size() - 1]);
  } // or return empty lists if nothing is found
  return {filteredEnergies, bandsExtrema};
}

std::tuple<std::vector<double>, std::vector<int>> Window::internalEnWindow(
    const Eigen::VectorXd &energies) const {

  std::vector<double> filteredEnergies;
  std::vector<int> bandsExtrema;
  std::vector<int> bandsIndices;
  double thisEnergy;
  for (int ib = 0; ib < numBands; ib++) {
    thisEnergy = energies(ib);
    if (thisEnergy < maxEnergy && thisEnergy > minEnergy) {
      filteredEnergies.push_back(thisEnergy);
      bandsIndices.push_back(ib);
    }
  }
  // set the band extrema
  if (!bandsIndices.empty()) {
    bandsExtrema.push_back(bandsIndices[0]);
    bandsExtrema.push_back(bandsIndices[bandsIndices.size() - 1]);
  } // or return empty lists if nothing is found
  return {filteredEnergies, bandsExtrema};
}

int Window::getMethodUsed() const {
  return method;
}

std::tuple<std::vector<double>, std::vector<int>> Window::internalMagWindow(
    const Eigen::VectorXd &energies, const Eigen::VectorXd &popMin,
    const Eigen::VectorXd &popMax) const {

  // need a window for using magnetic fields which
  // keeps all bands for the selected kpoints.
  // This is important because of the grad k
  // term we want to compute

  std::vector<double> energiesVec;
  std::vector<int> bandsExtrema;
  std::vector<int> bandsIndices;

  for (int ib = 0; ib < numBands; ib++) {
    energiesVec.push_back(energies(ib));
    if ((abs(popMin(ib)) > populationThreshold) ||
        (abs(popMax(ib)) > populationThreshold)) {
      bandsIndices.push_back(ib);
    }
  }
  if (!bandsIndices.empty()) {
    bandsExtrema.push_back(0);
    bandsExtrema.push_back(numBands - 1);
    return {energiesVec, bandsExtrema};
  } else {
    std::vector<double> filteredEnergies2;
    std::vector<int> bandsExtrema2;
    return {filteredEnergies2, bandsExtrema2};
  }
}


