#ifndef PHONONCONDUCTIVITY_H
#define PHONONCONDUCTIVITY_H

#include "observable.h"
#include "specific_heat.h"

class PhononThermalConductivity : public Observable {
public:
	PhononThermalConductivity(StatisticsSweep & statisticsSweep_,
			Crystal & crystal_,
			FullBandStructure<FullPoints> & bandStructure_);
	// copy constructor
	PhononThermalConductivity(const PhononThermalConductivity & that);
	// copy assignment
	PhononThermalConductivity & operator = (
			const PhononThermalConductivity & that);

	PhononThermalConductivity operator - (
			const PhononThermalConductivity & that);

	virtual void calcFromPopulation(VectorBTE & n);
	virtual void calcFromCanonicalPopulation(VectorBTE & f);
	void calcVariational(VectorBTE & af, VectorBTE & f, VectorBTE & scalingCG);
	void calcFromRelaxons(SpecificHeat & specificHeat, VectorBTE & relaxonV,
			VectorBTE & relaxationTimes);
	void print();
	void print(const int & iter);

protected:
	virtual int whichType();
	FullBandStructure<FullPoints> & bandStructure;
};

#endif
