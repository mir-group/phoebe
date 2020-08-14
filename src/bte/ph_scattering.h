#ifndef PHSCATTERING_H
#define PHSCATTERING_H

#include "interaction_3ph.h"
#include "phonon_h0.h"
#include "scattering.h"
#include "vector_bte.h"

/** class representing the phonon scattering matrix.
 * This class contains the logic to compute the phonon scattering matrix.
 * The parent class ScatteringMatrix instead contains the logic for managing
 * the operations with phonon distribution vectors.
 */
class PhScatteringMatrix : public ScatteringMatrix {
 public:
  /** Default constructor
   * @param context: the user-initialized variables.
   * @param statisticsSweep: the object containing the information on the
   * temperatures to be used in the calculation.
   * @param innerBandStructure: this is the bandstructure object used for
   * integrating the sum over q2 wavevectors.
   * @param outerBandStructure: this is the bandStructure object used for
   * integrating the sum over q1 wavevectors.
   * @param coupling3Ph: a pointer to the class handling the 3-phonon
   * interaction calculation.
   * @param h0: the object used for constructing phonon energies.
   *
   * Note: inner and outer bandstructures may be different, for example, if we
   * want to compute the phonon linewidths on a path, the outer bandstructure
   * contains the phonon dispersion relation along high symmetry lines, whereas
   * the inner bandstructure contains the dispersion relation on a full grid
   * of wavevectors. For transport calculations instead, inner=outer.
   * Especially in the case of computing linewidths on a path, it might be
   * necessary to compute phonon energies on a wavevector q3=q1+q2: for this
   * reason, we need to pass also the phonon h0 object.
   */
  PhScatteringMatrix(Context &context_, StatisticsSweep &statisticsSweep_,
                     BaseBandStructure &innerBandStructure_,
                     BaseBandStructure &outerBandStructure_,
                     Interaction3Ph *coupling3Ph_ = nullptr,
                     PhononH0 *h0 = nullptr);

  /** Copy constructor
   */
  PhScatteringMatrix(const PhScatteringMatrix &that);

  /** Copy assignment operator
   */
  PhScatteringMatrix &operator=(const PhScatteringMatrix &that);

 protected:
  Interaction3Ph *coupling3Ph;
  PhononH0 *h0;

  Eigen::VectorXd massVariance;
  bool doIsotopes;

  double boundaryLength;
  bool doBoundary;

  // implementation of the scattering matrix
  virtual void builder(VectorBTE *linewidth,
                       std::vector<VectorBTE> &inPopulations,
                       std::vector<VectorBTE> &outPopulations);
};

#endif
