#ifndef PH_SCATTERING_H
#define PH_SCATTERING_H

#include "ph_scattering_matrix.h"

 // friend functions for adding scattering rates,
 // these live in ph_scattering.cpp
 // TODO write docstrings for these

void addPhPhScattering(PhScatteringMatrix &matrix, Context &context,
                                std::vector<VectorBTE> &inPopulations,
                                std::vector<VectorBTE> &outPopulations,
                                int &switchCase,
                                std::vector<std::tuple<std::vector<int>, int>> qPairIterator,
                                Eigen::MatrixXd &innerBose, Eigen::MatrixXd &outerBose,
                                BaseBandStructure &innerBandStructure,
                                BaseBandStructure &outerBandStructure,
                                VectorBTE *linewidth);

void addIsotopeScattering(PhScatteringMatrix &matrix, Context &context,
                                std::vector<VectorBTE> &inPopulations,
                                std::vector<VectorBTE> &outPopulations, int &switchCase,
                                std::vector<std::tuple<std::vector<int>, int>> qPairIterator,
                                Eigen::MatrixXd &innerBose, Eigen::MatrixXd &outerBose,
                                BaseBandStructure &innerBandStructure,
                                BaseBandStructure &outerBandStructure,
                                VectorBTE *linewidth);

#endif
