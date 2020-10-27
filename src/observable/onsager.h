#ifndef ONSAGER_H
#define ONSAGER_H

#include "bandstructure.h"
#include "context.h"
#include "crystal.h"
#include "statistics_sweep.h"
#include "vector_bte.h"

/** Class to compute the electronic transport coefficients.
 */
class OnsagerCoefficients {
public:
  /** Constructor method
   */
  OnsagerCoefficients(StatisticsSweep &statisticsSweep_, Crystal &crystal_,
                      BaseBandStructure &bandStructure_, Context &context);

  /** Copy constructor
   */
  OnsagerCoefficients(const OnsagerCoefficients &that);

  /** Copy assignment operator
   */
  OnsagerCoefficients &operator=(const OnsagerCoefficients &that);

  /** Compute the thermal conductivity from the phonon populations
   * @param n: the phonon population out-of-equilibrium. Note that this
   * method uses the absolute value of phonon populations n.
   */
  void calcFromPopulation(VectorBTE &nE, VectorBTE &nT);

  /** Compute the thermal conductivity from canonical phonon populations
   * where the "canonical" population f is defined as n = bose(bose+1)f
   * where "bose" is the Bose-Einstein population.
   * Iterative schemes are probably computing this canonical population.
   * @param f: canonical phonon population.
   */
  void calcFromCanonicalPopulation(VectorBTE &fE, VectorBTE &fT);

  /** Prints to screen the thermal conductivity at various temperatures
   * in a a nicely formatted way.
   */
  void print();

  /** Short format for printing the electrical conductivity. To be used
   * for quickly evaluate the convergence of an iterative BTE solver.
   */
  void print(const int &iter);

  /** Outputs the quantity to a json file.
  * @param outFileName: string representing the name of the json file
  */
  void outputToJSON(std::string outFileName);

  void calcTransportCoefficients();

  void calcFromEPA();

  Eigen::Tensor<double,3> getElectricalConductivity();
  Eigen::Tensor<double,3> getThermalConductivity();

protected:
  StatisticsSweep &statisticsSweep;
  Crystal &crystal;
  BaseBandStructure &bandStructure;
  Context &context;

  int dimensionality;
  double spinFactor;
  int numCalcs;

  Eigen::Tensor<double, 3> sigma, seebeck, kappa, mobility;
  Eigen::Tensor<double, 3> LEE, LET, LTE, LTT;
};

#endif
