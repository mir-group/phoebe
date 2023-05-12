#ifndef PHEL_SCATTERING_MATRIX_H
#define PHEL_SCATTERING_MATRIX_H

#include "electron_h0_wannier.h"
#include "phonon_h0.h"
#include "scattering_matrix.h"
#include "vector_bte.h"
#include "phel_scattering.cpp"

/** This class describes the construction of the electron scattering matrix.
 * The most important part is the assembly of the electron-phonon scattering.
 * We also include boundary scattering effects.
 */
class PhElScatteringMatrix : public ScatteringMatrix {
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
  PhElScatteringMatrix(Context &context_, StatisticsSweep &statisticsSweep_,
                       BaseBandStructure &elBandStructure_,
                       BaseBandStructure &phBandStructure_,
                       InteractionElPhWan &couplingElPhWan_,
                       ElectronH0Wannier &electronH0_);

  /** Copy constructor
   * @param that: object to be copied
   */
  PhElScatteringMatrix(const PhElScatteringMatrix &that);

  /** Copy assignment
   *
   * @param that: object to be copied
   * @return a copy of ElScatteringMatrix
   */
  PhElScatteringMatrix &operator=(const PhElScatteringMatrix &that);

protected:
  InteractionElPhWan &couplingElPhWan;
  ElectronH0Wannier &electronH0;

  void builder(VectorBTE *linewidth, std::vector<VectorBTE> &inPopulations,
               std::vector<VectorBTE> &outPopulations) override;

  // to prevent mistakes in which these two outer and inner BS could
  // be accidentally swapped
  BaseBandStructure& getPhBandStructure() { return outerBandStructure; };
  BaseBandStructure& getElBandStructure() { return innerBandStructure; };

  std::vector<std::tuple<int, std::vector<int>>> getIteratorWavevectorPairs();
  std::vector<std::tuple<int, std::vector<int>>> getIrrWavevectorPairs();

  friend void addPhElScattering(ScatteringMatrix &matrix, Context &context); 

};

#endif
