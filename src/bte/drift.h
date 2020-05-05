#ifndef DRIFT_H
#define DRIFT_H

#include "vector_bte.h"

// should we have a base class?
//class BulkDriftOperator : public VectorBTE {
//};

class BulkTDrift : public VectorBTE {
public:
	BulkTDrift(Context & context, ActiveBandStructure & activeBandStructure);
};

class BulkEDrift : public VectorBTE {
public:
	BulkEDrift(Context & context, ActiveBandStructure & activeBandStructure);
};

#endif
