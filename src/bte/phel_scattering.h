#ifndef PHEL_SCATTERING_H
#define PHEL_SCATTERING_H

#include "phel_scattering_matrix.h"

 // friend functions for adding scattering rates,
 // these live in ph_scattering.cpp
 // TODO write docstrings for these

 void addPhElScattering(BasePhScatteringMatrix &matrix, Context &context,
                std::vector<std::tuple<int, std::vector<int>>> kqPairIterator,
                BaseBandStructure &phBandStructure,
                BaseBandStructure &elBandStructure,
                ElectronH0Wannier* electronH0,
                InteractionElPhWan *couplingElPhWan,
                VectorBTE *linewidth) {

#endif
