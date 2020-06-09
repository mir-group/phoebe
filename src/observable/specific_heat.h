#ifndef SPECIFICHEAT_H
#define SPECIFICHEAT_H

#include "observable.h"

class SpecificHeat : public Observable {
public:
	SpecificHeat(StatisticsSweep & statisticsSweep_,
			Crystal & crystal_,
			FullBandStructure<FullPoints> & bandStructure_);
	// copy constructor
	SpecificHeat(const SpecificHeat & that);
	// copy assignment
	SpecificHeat & operator = (const SpecificHeat & that);

	virtual void calc();
	void print();
	void print(const int & iter);

	double get(const long & imu, const long & it);

protected:
	virtual int whichType();
	FullBandStructure<FullPoints> & bandStructure;
};

#endif
