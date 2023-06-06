#ifndef GENERAL_SCATTERING_H
#define GENERAL_SCATTERING_H

#include "scattering_matrix.h"
#include "context.h"
#include "vector_bte.h"

// TODO document function
void addBoundaryScattering(ScatteringMatrix &matrix, Context &context,
                                std::vector<VectorBTE> &inPopulations,
                                std::vector<VectorBTE> &outPopulations,
                                int switchCase,
                                BaseBandStructure &outerBandStructure,
                                VectorBTE *linewidth);

#endif
