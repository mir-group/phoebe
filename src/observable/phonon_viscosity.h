#ifndef PHONONVISCOSITY_H
#define PHONONVISCOSITY_H

#include "observable.h"
#include "drift.h"
#include "ph_scattering.h"

class PhononViscosity : public Observable {
public:
	PhononViscosity(StatisticsSweep & statisticsSweep_, Crystal & crystal_,
			FullBandStructure<FullPoints> & bandStructure_);
	// copy constructor
	PhononViscosity(const PhononViscosity & that);
	// copy assignment
	PhononViscosity & operator = (const PhononViscosity & that);

	virtual void calcRTA(VectorBTE & n);
	void print();

	virtual void calcFromRelaxons(Vector0 & vector0, VectorBTE & relTimes,
			PhScatteringMatrix & sMatrix, Eigen::MatrixXd & eigenvectors);
protected:
	virtual int whichType();
	FullBandStructure<FullPoints> & bandStructure;
};

#endif
