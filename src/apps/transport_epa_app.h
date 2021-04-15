#ifndef EPA_TRANSPORT_APP_H
#define EPA_TRANSPORT_APP_H

#include "app.h"
#include "basevector_bte.h"
#include "delta_function.h"

/** App for computing the electron transport properties with the EPA theory
 */
class TransportEpaApp : public App {
public:
  void run(Context &context) override;
  void checkRequirements(Context &context) override;

private:
/** Auxiliary method to compute the tensor of velocities times the density of
 * states.
 *
 * @param context: object with the user-defined input parameters
 * @param bandStructure: object with the electronic band structure
 * @param energies: list of energies used to integrate transport properties.
 * @param tetrahedrons: a tetrahedronDeltaFunction object that will be used to
 * integrate the density of states.
 * @return tensor of velocities: an Eigen::Tensor of dimensions
 * (3,3,numEnergies)
 */
  Eigen::Tensor<double, 3> static calcEnergyProjVelocity(
      Context &context, FullBandStructure &bandStructure,
      const Eigen::VectorXd &energies, TetrahedronDeltaFunction &tetrahedrons);

  /** This method computes the electron lifetimes at the EPA level.
   *
   * @param context: object with the user-defined input parameters
   * @param statisticsSweep: object with values of temperature and chemical
   * potential
   * @param fullBandStructure: object with the electronic band structure on a
   * full grid of wavevectors
   * @param energies: list of energies at which the EPA lifetimes will be
   * computed
   * @param tetrahedrons: a tetrahedronDeltaFunction object that will be used to
   * integrate the Dirac-delta for energy conservation.
   * @return a BaseVectorBTE object containing the electron lifetimes.
   */
  static BaseVectorBTE getScatteringRates(Context &context,
                                   StatisticsSweep &statisticsSweep,
                                   FullBandStructure &fullBandStructure,
                                   Eigen::VectorXd &energies,
                                   TetrahedronDeltaFunction &tetrahedrons);

  /* helper function to output scattering rates as a function of energy*/
  void outputToJSON(const std::string &outFileName, BaseVectorBTE &scatteringRates,
                StatisticsSweep &statisticsSweep, int &numEnergies,
                Eigen::VectorXd &energiesEPA);

};

#endif
