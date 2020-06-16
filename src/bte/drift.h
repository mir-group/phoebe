#ifndef DRIFT_H
#define DRIFT_H

#include "vector_bte.h"
#include "specific_heat.h"

/** Object describing the thermal diffusion operator of the BTE
 *
 */
class BulkTDrift : public VectorBTE {
public:
	BulkTDrift(StatisticsSweep & statisticsSweep_,
			BaseBandStructure & bandStructure_,
			const long & dimensionality_=3);
};

/** Object describing the electric-field drift operator of the BTE
 *
 */
class BulkEDrift : public VectorBTE {
public:
	BulkEDrift(StatisticsSweep & statisticsSweep_,
			BaseBandStructure & bandStructure_,
			const long & dimensionality_=3);
};

/** Object describing the eigenvector with zero-eigenvalue of the symmetrized
 * scattering matrix Omega (notation of Cepellotti PRX (2016)) associated with
 * the energy conservation.
 */
class Vector0 : public VectorBTE {
public:
	Vector0(StatisticsSweep & statisticsSweep_,
			BaseBandStructure & bandStructure_, SpecificHeat & specificHeat);
};

#endif
