#ifndef EL_SCATTERING_H
#define EL_SCATTERING_H

#include "el_scattering_matrix.h"

// friend functions for adding scattering rates,
// these live in el_scattering.cpp

// TODO finish adding doxygen strings to these
/**
 * @param matrix: a el scattering matrix object
 * @param context: object with user parameters for this calculation
 * @param inPopulations:
 **/

void addElPhScattering(ElScatteringMatrix &matrix, Context &context,
                  std::vector<VectorBTE> &inPopulations,
                  std::vector<VectorBTE> &outPopulations,
                  std::vector<std::tuple<std::vector<int>, int>> kPairIterator,
                  int &switchCase,
                  Eigen::MatrixXd &innerFermi, Eigen::MatrixXd &outerBose,
                  BaseBandStructure &innerBandStructure,
                  BaseBandStructure &outerBandStructure, VectorBTE *linewidth);

#endif
