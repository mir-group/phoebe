#ifndef SWEEP_H
#define SWEEP_H

#include <memory>

#include "bandstructure.h"
#include "context.h"
#include "particle.h"
#include "utilities.h"

/** Container for temperature, chemical potential, doping, to be used.
 * To be used for loops over temperature and chemical potentials.
 */
struct CalcStatistics {
  double temperature;
  double chemicalPotential;
  double doping;
};

/** Object for controlling the loop over temperatures and chemical potentials.
 * Given temperatures and chemical potentials (or doping), we build a list of
 * pairs of temperatures and chemical potentials (and doping). This list can be
 * accessed to repeat calculations on a range of parameters.
 */
class StatisticsSweep {
 public:
  /** Constructor of the loop over temperatures and chemical potentials.
   * For phonons, the constructor simply stores the temperatures.
   * For electrons, the class finds chemical potentials from doping (or
   * viceversa, depending on user input), and builds a table of all
   * pairs of temperatures and chemical potentials.
   * @param context: object with user input.
   */
  StatisticsSweep(Context &context,
                  FullBandStructure *fullBandStructure = nullptr);

  /** Copy constructor
   */
  StatisticsSweep(const StatisticsSweep &that);

  /** Copy assignment
   */
  StatisticsSweep &operator=(const StatisticsSweep &that);

  /** returns the CalcStatistics object containing temperature,
   * chemical potential and temperature, given an index on the number of
   * "calculations" = numTemperature*numChemicalPotentials.
   * @param index: an integer, ranging in [0,numCalcs[
   * @return calcStatistics: a struct object with temperature, chemical
   * potential and doping concentration, for this particular index.
   */
  struct CalcStatistics getCalcStatistics(const long &index);

  /** returns the CalcStatistics object containing temperature,
   * chemical potential and temperature, given indices on temperature and
   * chemical potential.
   * @param iTemp: a strong-typed index, ranging in [0,numTemperatures[
   * @param iChemPot: a strong-typed index, ranging in
   * [0,numChemicalPotentials[
   * @return calcStatistics: a struct object with temperature, chemical
   * potential and doping concentration.
   */
  struct CalcStatistics getCalcStatistics(const TempIndex &iTemp,
                                          const ChemPotIndex &iChemPot);

  /** Returns the number of "calculations" i.e. the number of temperatures
   * times the number of chemical potentials we will compute.
   */
  long getNumCalcs();

  /** Returns for how many chemical potentials we are computing properties.
   */
  long getNumChemicalPotentials();

  /** Returns for how many temperatures we are computing properties.
   */
  long getNumTemperatures();

  /** Prints to screen the information on the temperature, doping and
   * chemical potentials to be computed in the current run.
   */
  void printInfo();

 protected:
  Particle particle;
  long numCalcs;
  Eigen::MatrixXd infoCalcs;
  long nTemp;
  long nChemPot;
  long nDop;
  bool isDistributed;

  // note: these three private methods are only meant to be used
  // during the construction of this class, as they depend on energies

  /** Computes the electronic chemical potential, given the doping
   * concentration.
   */
  double findChemicalPotentialFromDoping(const double &doping,
                                         const double &temperature);

  /** The opposite of the previous method, computes the doping concentration
   * given the electronic chemical potential.
   */
  double findDopingFromChemicalPotential(const double &chemicalPotential,
                                         const double &temperature);

  // Auxiliary function for finding the chemical potential.
  double fPop(const double &chemPot, const double &temp);

  // this block is temporary variables for electronic calculations
  const long maxIter = 100;
  long numPoints;
  long numBands;
  Eigen::VectorXd energies;
  double numElectronsDoped;
  double volume;
  double spinFactor;
  double occupiedStates;
  double fermiLevel;
};

#endif