#ifndef EL_SCATTERING_MATRIX_H
#define EL_SCATTERING_MATRIX_H

#include "electron_h0_wannier.h"
#include "interaction_elph.h"
#include "phonon_h0.h"
#include "scattering_matrix.h"
#include "vector_bte.h"
#include "general_scattering.cpp"

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
                     BaseBandStructure &outerBandStructure_, PhononH0 &h0,
                     InteractionElPhWan *couplingElPhWan_ = nullptr);

protected:
  InteractionElPhWan *couplingElPhWan;
  PhononH0 &h0;

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


  // friend functions for adding scattering rates, 
  // these live in el_scattering.cpp

 // TODO finish adding doxygen strings to these
  /** 
   * @param matrix: a el scattering matrix object
   * @param context: object with user parameters for this calculation
   * @param inPopulations: 
   * */ 
  friend void addBoundaryScattering(ScatteringMatrix &matrix, Context &context,
                                std::vector<VectorBTE> &inPopulations,
                                std::vector<VectorBTE> &outPopulations, 
                                BaseBandStructure &outerBandStructure, 
                                VectorBTE *linewidth);

  friend void addElPhScattering(ElScatteringMatrix &matrix, Context &context, 
                       std::vector<VectorBTE> &inPopulations,
                       std::vector<VectorBTE> &outPopulations, 
                       std::vector<std::tuple<std::vector<int>, int>> kPairIterator, 
                       int &switchCase,                                 
                       Eigen::MatrixXd &innerFermi, Eigen::MatrixXd &outerBose,
                       BaseBandStructure &innerBandStructure,
                       BaseBandStructure &outerBandStructure, VectorBTE *linewidth); 

};

#endif
