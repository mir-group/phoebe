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
    method = magnetotransport;
  } else if (inMethod == "nothing") {
    method = nothing;
  } else {
    Error("Unrecognized method called in Window()");
  }

  if (method == population || method == magnetotransport) {
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
  else if (method == magnetotransport) {

    Eigen::VectorXd popMin(numBands), popMax(numBands);
    for (int ib = 0; ib < numBands; ib++) {
      popMin(ib) = particle.getPopPopPm1(energies(ib), temperatureMax,
                                         chemicalPotentialMin);
      popMax(ib) = particle.getPopPopPm1(energies(ib), temperatureMax,
                                         chemicalPotentialMax);
    }
    if(particle.isPhonon()) {
      return internalPopWindow(energies, popMin, popMax);
    } else {
      return internalMagWindow(energies);
    }
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
    return std::make_tuple(filteredEnergies, bandExtrema);
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
  return std::make_tuple(filteredEnergies, bandsExtrema);
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
  return std::make_tuple(filteredEnergies, bandsExtrema);
}

int Window::getMethodUsed() const {
  return method;
}

std::vector<int> Window::prepareMagWindow(FullBandStructure& fbs,
                                double chemPot, double temperature) {

  int numBands = fbs.getNumBands();
  std::vector<int> magBands;
  // E = ln(1/popThresh)(kbT) + mu
  double maxEn = log(1e9)*(temperature) + chemPot;
  double minEn = -log(1e9)*(temperature) + chemPot;

  for(int ib = 0; ib < numBands; ib++) {
    Eigen::VectorXd bandEns = fbs.getBandEnergies(ib);
    double maxBandEn = bandEns.maxCoeff();
    double minBandEn = bandEns.minCoeff();
    // check if ranged intersect
    if(minBandEn <= maxEn && minEn <= maxBandEn) {
      magBands.push_back(ib);
    }
  }
  std::vector<int> bandsExtrema;
  bandsExtrema.push_back(*std::min_element(magBands.begin(),magBands.end()));
  bandsExtrema.push_back(*std::max_element(magBands.begin(),magBands.end()));
  return bandsExtrema;
}

void Window::setMagBandsExtrema(std::vector<int> bandsExtrema) {
   magBandsExtrema.push_back(bandsExtrema[0]);
   magBandsExtrema.push_back(bandsExtrema[1]);
}

std::tuple<std::vector<double>, std::vector<int>> Window::internalMagWindow(
    const Eigen::VectorXd &energies) const {

  // need a window for using magnetic fields which
  // keeps all bands for the selected kpoints.
  // This is important because of the grad k
  // term we want to compute. We pre-determined which bands to
  // keep in the "prepareMagWindow" function, and set them for all
  // mpi processes in "setMagBandsExtrema". Now we simply apply them.

  std::vector<double> filteredEnergies;

  // we only save the relevant band energies
  for (int ib = magBandsExtrema[0]; ib < magBandsExtrema[1]+1; ib++) {
    filteredEnergies.push_back(energies(ib));
  }
  return {filteredEnergies,magBandsExtrema};
}


