#ifndef PHSCATTERING_H
#define PHSCATTERING_H

#include "interaction_3ph.h"
#include "phonon_h0.h"
#include "scattering.h"
#include "vector_bte.h"

class PhScatteringMatrix: public ScatteringMatrix {
public:
    PhScatteringMatrix(Context &context_, StatisticsSweep &statisticsSweep_,
            BaseBandStructure &innerBandStructure_,
            BaseBandStructure &outerBandStructure_,
            Interaction3Ph *coupling3Ph_ = nullptr, PhononH0 *h0 = nullptr);

    PhScatteringMatrix(const PhScatteringMatrix &that);

    PhScatteringMatrix& operator=(const PhScatteringMatrix &that);

protected:
    Interaction3Ph *coupling3Ph;
    PhononH0 *h0;

    Eigen::VectorXd massVariance;
    bool doIsotopes;

    double boundaryLength;
    bool doBoundary;

    virtual void builder(Matrix<double> &matrix, VectorBTE *linewidth,
            VectorBTE *inPopulation, VectorBTE *outPopulation);
};

#endif
