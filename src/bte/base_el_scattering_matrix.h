#ifndef BASE_EL_SCATTERING_MATRIX_H
#define BASE_EL_SCATTERING_MATRIX_H

#include "phonon_h0.h"
#include "interaction_elph.h"
#include "electron_h0_wannier.h"
#include "scattering_matrix.h"
#include "vector_bte.h"

/** class representing the electron scattering matrix.
 * This class holds the absolute essentials (electronH0) that any
 * electron scattering matrix must have
 */
class BaseElScatteringMatrix : virtual public ScatteringMatrix {
 public:

  /** Default constructor */
  //BaseElScatteringMatrix(); 

  /** Constructor that just calls the constructor of scattering matrix */
  BaseElScatteringMatrix(Context &context_,
                                   StatisticsSweep &statisticsSweep_,
                                   BaseBandStructure &innerBandStructure_,
                                   BaseBandStructure &outerBandStructure_);

  //~BaseElScatteringMatrix() = default;

  /** Getter functions to return intial and final banstructures */
  //virtual BaseBandStructure initialElBandStructure(); 
  //virtual BaseBandStructure finalElBandStructure(); 

  //ElectronH0Wannier electronH0;

 protected:

  // TODO write docstrings for these
  // friend functions for adding scattering rates,
  // these live in el_scattering.cpp

  friend void addElPhScattering(BaseElScatteringMatrix &matrix, Context &context,
                       std::vector<VectorBTE> &inPopulations,
                       std::vector<VectorBTE> &outPopulations,
                       int &switchCase,
                       std::vector<std::tuple<std::vector<int>, int>> kPairIterator,
                       Eigen::MatrixXd &innerFermi, Eigen::MatrixXd &outerBose,
                       BaseBandStructure &innerBandStructure,
                       BaseBandStructure &outerBandStructure,
                       PhononH0 &phononH0,
                       InteractionElPhWan *couplingElPhWan,
                       VectorBTE *linewidth);

  friend void addDragTerm(BaseElScatteringMatrix &matrix, Context &context,
                  std::vector<std::tuple<std::vector<int>, int>> kqPairIterator,
                  int dragTermType,
                  ElectronH0Wannier* electronH0,
                  InteractionElPhWan *couplingElPhWan,
                  BaseBandStructure &innerBandStructure,
                  BaseBandStructure &outerBandStructure, VectorBTE *linewidth);

};

#endif
