#ifndef COUPLED_SCATTERING_MATRIX_H
#define COUPLED_SCATTERING_MATRIX_H

// TODO I think these are likely inherited
//#include "interaction_3ph.h"
#include "ph_scattering_matrix.h"
#include "el_scattering_matrix.h"
#include "phel_scattering_matrix.h"
#include "vector_bte.h"

/** class representing the combined scattering matrix.
 * This class contains the logic to compute the combined scattering matrix.
 * For a coupled BTE solve.
 * The parent class ScatteringMatrix instead contains the logic for managing
 * the operations with distribution vectors.
 */
class CoupledScatteringMatrix : public ElScatteringMatrix,
                                public PhScatteringMatrix,
                                public PhElScatteringMatrix {
                                //public ScatteringMatrix {
 public:

// TODO update all these comments

  /** Default constructor
   * @param context: the user-initialized variables.
   * @param statisticsSweep: the object containing the information on the
   * temperatures to be used in the calculation.
   * @param innerBandStructure: this is the band structure object used for
   * integrating the sum over final state wavevectors.
   * @param outerBandStructure: this is the bandStructure object used for
   * integrating the sum over initial state wavevectors.
   * @param coupling3Ph: a pointer to the class handling the 3-phonon
   * interaction calculation.
   * @param couplingElPh: a pointer to the class handling the el-ph
   * interaction calculation.
   * @param PhononH0: the object used for constructing phonon energies.
   * @param ElectronH0: the object used for constructing electron energies.
   *
   * Note: For transport calculations inner=outer.
   * Other scattering matrices allow for the possibility that they aren't equal,
   * but this matrix will only be used for transport.
   */
 CoupledScatteringMatrix(Context &context_, StatisticsSweep &statisticsSweep_,
                    BaseBandStructure &innerBandStructure_,
                    BaseBandStructure &outerBandStructure_,
                    Interaction3Ph *coupling3Ph_ = nullptr,
                    InteractionElPhWan *couplingElPh_ = nullptr,
                    ElectronH0Wannier *electronH0_ = nullptr,
                    PhononH0 *phononH0_ = nullptr);

  // TODO we will need to override the simple scattering matrix version of this function
  // as this one will need to behave differently than the others.
  // we may want to output each kind of linewidths, etc, for testing?
  /** Outputs the quantity to a json file.
   * @param outFileName: string representing the name of the json file
   */
  void outputToJSON(const std::string &outFileName);

  /* TODO we will need to override this function to resolve the confusion of
  * multiple inheritance of this function. This override should just call
  * the ScatteringMatrix base version of this
  */
  void setup();

// TODO -- might be smarter to add functions to add and calculate different effects one at a time,
 // so that we do not need to store the interaction/coupling info for both the entire time
 // If we need the matrix-vector product or any such thing where we do not
 // want to keep the whole matrix in memory for some reason, we will need these the whole time

 // TODO what are the cases where outer != inner band structure for the elph and phph cases?
 // I think this is only for linewidths on a path, which we do not do here, so we should
 // just be able to pass directly in one phonon and one electron band structure

 // TODO check on the functions that symmetrize this matrix's components

 protected:

  // inherit these from the other El, Ph, PhEl classes
  //Interaction3Ph *coupling3Ph;
  //InteractionElPh *couplingElPh;

  //PhononH0 *phononH0;
  //ElectronH0 *electronH0;

  // implementation of the scattering matrix
  void builder(VectorBTE *linewidth,
               std::vector<VectorBTE> &inPopulations,
               std::vector<VectorBTE> &outPopulations);

  // TODO we will need to write a very smart indexing funciton for this object.
  // it should ideally take the indexing functions from the el and ph scattering matrices
  // and

  // TODO another smart idea would be to make PhEl scattering take a PhScatterign matrix,
  // and make PhElScattering matrix inherit from PhScattering, plus then have the additional
  // elph interaction object? This idea is not fully formed but it coudl be useful

  // TODO switch case should be added to all of these, basically.

  // friend functions for adding scattering rates,
  // these live in ph_scattering.cpp
  // TODO write docstrings for these

// TODO I comment these out because I speculate that they are inherited from the other smatrices
/*
  friend void addBoundaryScattering(ScatteringMatrix &matrix, Context &context,
                                std::vector<VectorBTE> &inPopulations,
                                std::vector<VectorBTE> &outPopulations,
                                BaseBandStructure &outerBandStructure);

  friend void addPhPhScattering(PhScatteringMatrix &matrix, Context &context,
                                std::vector<VectorBTE> &inPopulations,
                                std::vector<VectorBTE> &outPopulations,
                                int &switchCase,
                                std::vector<std::tuple<std::vector<int>, int>> qPairIterator,
                                Eigen::MatrixXd &innerBose, Eigen::MatrixXd &outerBose,
                                BaseBandStructure &innerBandStructure,
                                BaseBandStructure &outerBandStructure);

  // friend void addPhElScattering();
  friend void addIsotopeScattering(PhScatteringMatrix &matrix, Context &context,
                                std::vector<VectorBTE> &inPopulations,
                                std::vector<VectorBTE> &outPopulations, int &switchCase,
                                std::vector<std::tuple<std::vector<int>, int>> qPairIterator,
                                Eigen::MatrixXd &innerBose, Eigen::MatrixXd &outerBose,
                                BaseBandStructure &innerBandStructure,
                                BaseBandStructure &outerBandStructure);

  friend void addElPhScattering(ElScatteringMatrix &matrix, Context &context,
                                std::vector<VectorBTE> &inPopulations,
                                std::vector<VectorBTE> &outPopulations,
                                std::vector<std::tuple<std::vector<int>, int>> kPairIterator,
                                int &switchCase,
                                Eigen::MatrixXd &innerFermi,
                                BaseBandStructure &innerBandStructure,
                                BaseBandStructure &outerBandStructure);
*/
};

#endif
