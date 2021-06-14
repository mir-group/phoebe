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
   * @param crystal: the crystal used in the calcuation, to divide by volume
   * @return a BaseVectorBTE object containing the electron lifetimes.
   */
  static BaseVectorBTE getScatteringRates(Context &context,
                                          StatisticsSweep &statisticsSweep,
                                          const Eigen::VectorXd &energies,
                                          Crystal &crystal,
                                          const Eigen::VectorXd &dos);

  /** Auxiliary method to compute the tensor of velocities times the density of
   * states.
   *
   * @param context: object with the user-defined input parameters
   * @param bandStructure: object with the electronic band structure
   * @param energies: list of energies used to integrate transport properties.
   * @param tetrahedrons: a tetrahedronDeltaFunction object that will be used to
   * integrate the density of states.
   * @return tuple(0) tensor of velocities: an Eigen::Tensor of dimensions
   * (3,3,numEnergies)
   * @return tuple(1) density of states, aligned with energies
   */
  std::tuple<
      Eigen::Tensor<double, 3>,
      Eigen::VectorXd> static calcEnergyProjVelocity(Context &context,
                                                     FullBandStructure
                                                         &bandStructure,
                                                     const Eigen::VectorXd
                                                         &energies);

  /* helper function to output scattering rates as a function of energy*/
  static void outputToJSON(const std::string &outFileName,
                           BaseVectorBTE &scatteringRates,
                           StatisticsSweep &statisticsSweep, int &numEnergies,
                           Eigen::VectorXd &energiesEPA);
};

#endif
