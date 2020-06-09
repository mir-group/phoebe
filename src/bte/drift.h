#ifndef DRIFT_H
#define DRIFT_H

#include "vector_bte.h"
#include "specific_heat.h"

// should we have a base class?
//class BulkDriftOperator : public VectorBTE {
//};

class BulkTDrift : public VectorBTE {
public:
	BulkTDrift(StatisticsSweep & statisticsSweep_,
			FullBandStructure<FullPoints> & bandStructure_,
			const long & dimensionality_=3);
};

class BulkEDrift : public VectorBTE {
public:
	BulkEDrift(StatisticsSweep & statisticsSweep_,
			FullBandStructure<FullPoints> & bandStructure_,
			const long & dimensionality_=3);
};

class Vector0 : public VectorBTE {
public:
	Vector0(StatisticsSweep & statisticsSweep_,
			FullBandStructure<FullPoints> & bandStructure_,
			SpecificHeat & specificHeat);
};

#endif
