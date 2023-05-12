#ifndef PH_SCATTERING_MATRIX_H
#define PH_SCATTERING_MATRIX_H

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
   * @param innerBandStructure: this is the band structure object used for
   * integrating the sum over q2 wavevectors.
   * @param outerBandStructure: this is the bandStructure object used for
   * integrating the sum over q1 wavevectors.
   * @param coupling3Ph: a pointer to the class handling the 3-phonon
   * interaction calculation.
   * @param h0: the object used for constructing phonon energies.
   *
   * Note: inner and outer band structures may be different, for example, if we
   * want to compute the phonon linewidths on a path, the outer band structure
   * contains the phonon dispersion relation along high symmetry lines, whereas
   * the inner band structure contains the dispersion relation on a full grid
   * of wavevectors. For transport calculations instead, inner=outer.
   * Especially in the case of computing linewidths on a path, it might be
   * necessary to compute phonon energies on a wavevector q3=q1+q2: for this
   * reason, we need to pass also the phonon h0 object.
   */
 CScatteringMatrix(Context &context_, StatisticsSweep &statisticsSweep_,
                     BaseBandStructure &innerBandStructure_,
                     BaseBandStructure &outerBandStructure_,
                     Interaction3Ph *coupling3Ph_ = nullptr,
                     InteractionElPh *couplingElPh_ = nullptr,
                     PhononH0 *phononH0 = nullptr, ElectronH0 *electronH0);

// TODO -- might be smarter to add functions to add and calculate different effects one at a time, 
 // so that we do not need to store the interaction/coupling info for both the entire time

 // TODO what are the cases where outer != inner band structure for the elph and phph cases? 
 // I think this is only for linewidths on a path, which we do not do here, so we should
 // just be able to pass directly in one phonon and one electron band structure

 protected:

  Interaction3Ph *coupling3Ph;
  InteractionElPh *couplingElPh;

  PhononH0 *phononH0;
  PhononH0 *electronH0;

  // implementation of the scattering matrix
  void builder(VectorBTE *linewidth,
               std::vector<VectorBTE> &inPopulations,
               std::vector<VectorBTE> &outPopulations) override;

 // friend functions for adding scattering rates, 
 // these live in ph_scattering.cpp
 // TODO write docstrings for these
 friend void addBoundaryScattering(PhScatteringMatrix &matrix, Context &context,
                                 std::vector<VectorBTE> &inPopulations,
                                 std::vector<VectorBTE> &outPopulations, int &switchCase);

 friend void addPhPhScattering(PhScatteringMatrix &matrix, Context &context, 
                                 std::vector<VectorBTE> &inPopulations,
                                 std::vector<VectorBTE> &outPopulations, int &switchCase, 
                                 std::vector<std::tuple<std::vector<int>, int>> qPairIterator); 
 // friend void addPhElScattering(); 
 friend void addIsotopeScattering(PhScatteringMatrix &matrix, Context &context, 
                                 std::vector<VectorBTE> &inPopulations,
                                 std::vector<VectorBTE> &outPopulations, int &switchCase, 
                                 std::vector<std::tuple<std::vector<int>, int>> qPairIterator);


};

#endif
