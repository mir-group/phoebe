#ifndef INTERACTION_EPA_H
#define INTERACTION_EPA_H

#include "context.h"
#include "eigen.h"

/** Class used to store and access the values of the electron-phonon coupling
 * at the EPA level of theory.
 *
 * The EPA average has been already performed elsewhere (in ElPhQeToPhoebeApp)
 */
class InteractionEpa {

private:
  Eigen::VectorXd elEnergies;
  Eigen::VectorXd phEnergies;
  Eigen::Tensor<double, 3> elPhMatAverage;
  double binSize;

public:
  /** default constructor.
   *
   * @param elEnergies_ : histogram of electronic energy values on top of which
   * the electron-phonon coupling has been averaged.
   * @param phEnergies_ : histogram of phononic energy values on top of which
   * the electron-phonon coupling has been averaged.
   * @param elPhMatAverage_: tensor with the EPA averaged values.
   */
  InteractionEpa(Eigen::VectorXd &elEnergies_, Eigen::VectorXd &phEnergies_,
                 Eigen::Tensor<double, 3> &elPhMatAverage_);

  /** copy constructor
   */
  InteractionEpa(const InteractionEpa &that);

  /** Copy assignment
   */
  InteractionEpa &operator=(const InteractionEpa &that);

  /** get the values of phonon energies used to compute the EPA coupling.
   */
  Eigen::VectorXd getPhEnergies();

  /** get the values of electron energies used to compute the EPA coupling.
   */
  Eigen::VectorXd getElEnergies();

  /** Get the values of the EPA coupling
   *
   * @param nu: phonon mode index
   * @param enI: energy of initial state
   * @param enF: energy of final state (k+q)
   * @return EPA approximation for |g|^2
   */
  double getCoupling(const int &nu, const double &enI, const double &enF);

  /** Function used to parse the EPA values from file
   *
   * @param context: object with user-defined parameters.
   * @return InteractionEPA: the instance of the class with the EPA coupling
   */
  static InteractionEpa parseEpaCoupling(Context &context);
};

#endif
