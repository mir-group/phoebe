#ifndef ONSAGER_H
#define ONSAGER_H

#include "bandstructure.h"
#include "context.h"
#include "crystal.h"
#include "statistics_sweep.h"
#include "vector_bte.h"
#include "el_scattering.h"

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

  /** Compute the transport coefficients from the electron populations
   * @param n: the electron population out-of-equilibrium.
   */
  void calcFromPopulation(VectorBTE &nE, VectorBTE &nT);

  /** Compute the transport coefficients from the canonical electron populations
   * where the canonical population f is related to the population n as
   * n = f0 + f0 (1-f0) f, with f0 being the Fermi--Dirac population.
   *
   * @param f: the canonical electron population out-of-equilibrium.
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
  void outputToJSON(const std::string& outFileName);

  void calcTransportCoefficients();

  void calcFromEPA(BaseVectorBTE &scatteringRates, Eigen::Tensor<double,3> &energyProjVelocity, Eigen::VectorXd &energies, double &energyStep, Particle &particle);

  void calcFromRelaxons(Eigen::VectorXd &eigenvalues,
                        ParallelMatrix<double> &eigenvectors,
                        ElScatteringMatrix &scatteringMatrix);

  void calcVariational(VectorBTE &afE, VectorBTE &afT,
                       VectorBTE &fE, VectorBTE &fT,
                       VectorBTE &nE, VectorBTE &nT,
                       VectorBTE &scalingCG);

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
