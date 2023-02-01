#ifndef EL_SCATTERING_H
#define EL_SCATTERING_H

#include "electron_h0_wannier.h"
#include "interaction_elph.h"
#include "phonon_h0.h"
#include "scattering.h"
#include "vector_bte.h"
#include "interaction_4el.h"

/** This class describes the construction of the electron scattering matrix.
 * The most important part is the assembly of the electron-phonon scattering.
 * We also include boundary scattering effects.
 */
class ElScatteringMatrix : public ScatteringMatrix {
public:
  /** Default constructor
   *
   * @param context_: object with user parameters for this calculation
   * @param statisticsSweep_: object with values for temperature and chemical
   * potential
   * @param innerBandStructure_: this is the band structure used for the
   * integration of lifetimes/scattering rates
   * @param outerBandStructure_: this is the band structure used to define on
   * which points to compute the lifetimes/scattering rates. For transport
   * properties outer=inner, but may differ e.g. when computing lifetimes on a
   * path
   * @param h0: phonon hamiltonian used to compute phonon energies and
   * eigenvectors.
   * @param couplingElPhWan_: object with the electron-phonon coupling.
   */
  ElScatteringMatrix(Context &context_, StatisticsSweep &statisticsSweep_,
                     BaseBandStructure &innerBandStructure_,
                     BaseBandStructure &outerBandStructure_);

  /** Copy constructor
   * @param that: object to be copied
   */
  ElScatteringMatrix(const ElScatteringMatrix &that);

  /** Copy assignment
   *
   * @param that: object to be copied
   * @return a copy of ElScatteringMatrix
   */
  ElScatteringMatrix &operator=(const ElScatteringMatrix &that);

  /** Add electron-phonon interaction object, if elph interactions are
   * to be used
   * @param context
   * @param std::shared_ptr couplingElPhWan: shared ptr to elph coupling object
   * @param std::shared_ptr phononH0: phononH0 object associated with elph int
   * */
  void addElPhInteraction(Context &context,
        std::shared_ptr<InteractionElPhWan>, PhononH0* phononH0);

  /** Add electron-electron interaction object, if elel interactions are
   * to be used
   * @param context
   * @param std::shared_ptr coupling4El: shared ptr to elel coupling object
   * */
  void add4ElInteraction(Context &context, std::shared_ptr<Interaction4El>);

protected:
  // these are shared pointers because we want to be able to optionally
  // create and add these interaction objects in apps,
  // possibly allowing them to go out of scope in the apps
  std::shared_ptr<InteractionElPhWan> couplingElPhWan;
  std::shared_ptr<Interaction4El> coupling4El;
  std::shared_ptr<PhononH0> h0;

  double boundaryLength;
  bool doBoundary;

  /** Function with the detailed calculation of the scattering matrix.
   *
   * Note: this function is computing the symmetrized scattering matrix
   * $\tilde{\Omega}$.
   * As a result, only use this with the appropriately symmetrized BTE
   *
   * @param linewidth
   * @param inPopulations
   * @param outPopulations
   */
  void builder(VectorBTE *linewidth, std::vector<VectorBTE> &inPopulations,
               std::vector<VectorBTE> &outPopulations) override;
};

#endif
