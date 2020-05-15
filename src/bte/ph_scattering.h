#ifndef PHSCATTERING_H
#define PHSCATTERING_H

#include "scattering.h"
#include "vector_bte.h"

class PhScatteringMatrix : public ScatteringMatrix {
public:
	PhScatteringMatrix(Context & context_,
			ActiveBandStructure & bandStructure_,
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
	virtual VectorBTE builderManager(VectorBTE * inPopulation);
	VectorBTE builder3Ph(Eigen::MatrixXd * matrix=nullptr,
			VectorBTE * inPopulation=nullptr);
}

#endif
