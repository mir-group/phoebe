#include "interaction_epa.h"
#include "exceptions.h"
#include <fstream>
#include <algorithm>

// default constructor
InteractionEpa::InteractionEpa(Eigen::VectorXd &elEnergies_,
                               Eigen::VectorXd &phEnergies_,
                               Eigen::Tensor<double, 3> &elPhMatAverage_)
    : elEnergies(elEnergies_), phEnergies(phEnergies_),
      elPhMatAverage(elPhMatAverage_) {
  assert(elPhMatAverage.dimension(1) == elPhMatAverage.dimension(2));
  assert(elPhMatAverage.dimension(1) == elEnergies.size());
  assert(elPhMatAverage.dimension(0) == phEnergies.size());

  binSize = 1.;
  if (elEnergies.size() > 1) {
    binSize = elEnergies(1) - elEnergies(0);
  }
}

// copy constructor
InteractionEpa::InteractionEpa(const InteractionEpa &that)
    : elEnergies(that.elEnergies), phEnergies(that.phEnergies),
      elPhMatAverage(that.elPhMatAverage), binSize(that.binSize) {}

// overload the assignment operator
InteractionEpa &InteractionEpa::operator=(const InteractionEpa &that) {
  if (this != &that) {
    elEnergies = that.elEnergies;
    phEnergies = that.phEnergies;
    elPhMatAverage = that.elPhMatAverage;
    binSize = that.binSize;
  }
  return *this;
}

Eigen::VectorXd InteractionEpa::getElEnergies() { return elEnergies; }

Eigen::VectorXd InteractionEpa::getPhEnergies() { return phEnergies; }

double InteractionEpa::getCoupling(const int &nu, const double &enI, const double &enF) {
  // find index of the energy in the bins of the elph energies
  int j = int(std::round((enI - elEnergies(0)) / binSize));
  int k = int(std::round((enF - elEnergies(0)) / binSize));

  // check and fold indices within bounds:
  j = std::min(int(elEnergies.size())-1,j);
  k = std::min(int(elEnergies.size())-1,k);
  j = std::max(0,j);
  k = std::max(0,k);

  return elPhMatAverage(nu, j, k);
}

InteractionEpa InteractionEpa::parseEpaCoupling(Context &context) {
  // open epa.elph file for reading
  std::ifstream infile(context.getEpaFileName());
  if (!infile) {
    Error("epa.elph file not found");
  }

  int numElectrons, numSpin;
  infile >> numElectrons >> numSpin;
  if (numSpin != 1) {
    Error("Spin not supported in EPA coupling");
  }
  context.setNumOccupiedStates(numElectrons);

  // numPhFreq: number of average phonon frequencies (equal to the number of
  // phonon branches)
  int numModes;
  infile >> numModes;
  // phFreqAverage - vector containing average phonon frequencies for numPhFreq
  // branches (in cm^-1)
  Eigen::VectorXd phEnergies(numModes);
  for (int i = 0; i < numModes; i++) {
    infile >> phEnergies(i);
  }

  int numEnergies;
  infile >> numEnergies;
  // energies - vector containing the electronic energies
  // at which the el-ph coupling has been computed
  Eigen::VectorXd energies(numEnergies);
  for (int i = 0; i < numEnergies; i++) {
    infile >> energies(i);
  }

  // elPhMatAverage - tensor containing averaged squared electron-phonon matrix
  // elements for each phonon branch and each energy bin for valence and
  // conduction bands
  Eigen::Tensor<double, 3> elPhMatAverage(numModes, numEnergies, numEnergies);
  elPhMatAverage.setZero();
  for (auto i = 0; i < numModes; ++i) {
    for (auto j = 0; j < numEnergies; ++j) {
      for (auto k = 0; k < numEnergies; ++k) {
        infile >> elPhMatAverage(i, j, k);
      }
    }
  }

  InteractionEpa interactionEpa(energies, phEnergies, elPhMatAverage);
  return interactionEpa;
}
