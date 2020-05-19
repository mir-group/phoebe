#ifndef PHSCATTERING_H
#define PHSCATTERING_H

#include "scattering.h"
#include "vector_bte.h"
#include "interaction_3ph.h"

class PhScatteringMatrix : public ScatteringMatrix {
public:
	PhScatteringMatrix(Context & context_,
			StatisticsSweep & statisticsSweep_,
			DeltaFunction * smearing_,
			FullBandStructure<FullPoints> & innerBandStructure_,
			FullBandStructure<FullPoints> & outerBandStructure_,
			Interaction3Ph * coupling3Ph_=nullptr);
//			InteractionIsotope * couplingIsotope_=nullptr,
//			InteractionBoundary * couplingBoundary_=nullptr

	PhScatteringMatrix(const PhScatteringMatrix & that);
	PhScatteringMatrix & operator=(const PhScatteringMatrix & that);

protected:
	Interaction3Ph * coupling3Ph = nullptr;
//	InteractionIsotope * couplingIsotope = nullptr;
//	InteractionBoundary * couplingBoundary = nullptr;

	virtual void builder(Eigen::MatrixXd * matrix, VectorBTE * linewidth,
			VectorBTE * inPopulation, VectorBTE * outPopulation);
};

#endif
