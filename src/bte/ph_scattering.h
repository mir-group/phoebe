#ifndef PHSCATTERING_H
#define PHSCATTERING_H

#include "scattering.h"
#include "vector_bte.h"

class PhScatteringMatrix : public ScatteringMatrix {
public:
	PhScatteringMatrix(Context & context_,
			StatisticsSweep & statisticsSweep_
			FullBandStructure & innerBandStructure_,
			FullBandStructure * outerBandStructure_ = nullptr,
			Interaction3Ph * coupling3Ph_=nullptr
//			InteractionIsotope * couplingIsotope_=nullptr,
//			InteractionBoundary * couplingBoundary_=nullptr
			);
	PhScatteringMatrix(const PhScatteringMatrix & that);
	PhScatteringMatrix & operator=(const PhScatteringMatrix & that);

private:
	Interaction3Ph * coupling3Ph = nullptr;
//	InteractionIsotope * couplingIsotope = nullptr;
//	InteractionBoundary * couplingBoundary = nullptr;

	virtual void builder(Eigen::MatrixXd * matrix, VectorBTE * linewidth,
			Vector3BTE * inPopulation, VectorBTE * outPopulation);
};

#endif
