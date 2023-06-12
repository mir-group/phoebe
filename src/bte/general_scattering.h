#ifndef GENERAL_SCATTERING_H
#define GENERAL_SCATTERING_H

#include "scattering_matrix.h"
#include "context.h"
#include "vector_bte.h"

// TODO document function
// NOTE: this requires an allReduceSum afterwards on the linewidhts, as is
// currently done after it is called in electron and phonon scattering matrices
void addBoundaryScattering(ScatteringMatrix &matrix, Context &context,
                                std::vector<VectorBTE> &inPopulations,
                                std::vector<VectorBTE> &outPopulations,
                                int switchCase,
                                BaseBandStructure &outerBandStructure,
                                VectorBTE *linewidth);

#endif
