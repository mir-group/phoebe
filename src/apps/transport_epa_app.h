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
   * @param energies: list of energies at which the EPA lifetimes will be
   * computed
   * @param crystal: the crystal used in the calculation, for integrals
   * normalization by volume
   * @param dos: an Eigen::VectorXd containing the density of states computed on
   * the vector of energies
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
  /** Helper function to output scattering rates as a function of energy
   *
   * @param outFileName: name of output JSON file.
   * @param scatteringRates: vector of scattering rates computed at the energies stored in energiesEPA
   * @param statisticsSweep: StatisticsSweep object containing info on temperature and chemical potential
   * @param energiesEPA: array with the values of energies used in the EPA calculation
   */
  static void outputToJSON(const std::string &outFileName,
                           BaseVectorBTE &scatteringRates,
                           StatisticsSweep &statisticsSweep,
                           Eigen::VectorXd &energiesEPA);
};

#endif
