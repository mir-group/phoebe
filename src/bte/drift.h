#ifndef DRIFT_H
#define DRIFT_H

#include "vector_bte.h"

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

#endif
