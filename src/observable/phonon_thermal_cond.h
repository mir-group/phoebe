#ifndef PHONONCONDUCTIVITY_H
#define PHONONCONDUCTIVITY_H

#include "observable.h"

class PhononThermalConductivity : public Observable {
public:
	PhononThermalConductivity(Context & context_, Crystal & crystal_,
			FullBandStructure<FullPoints> & bandStructure_);
	virtual void calcFromPopulation(VectorBTE & n);
//	void calcVariational(VectorBTE & af, VectorBTE & f, VectorBTE & b);
	void print();
protected:
	FullBandStructure<FullPoints> & bandStructure;
	int type = is2Tensor;
};

#endif
