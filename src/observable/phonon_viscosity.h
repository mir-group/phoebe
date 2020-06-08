#ifndef PHONONVISCOSITY_H
#define PHONONVISCOSITY_H

#include "observable.h"

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

protected:
	virtual int whichType();
	FullBandStructure<FullPoints> & bandStructure;
};

#endif
