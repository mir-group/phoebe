#ifndef DRIFT_H
#define DRIFT_H

#include "specific_heat.h"
#include "vector_bte.h"

/** Object describing the thermal diffusion operator of the BTE
 */
class BulkTDrift : public VectorBTE {
public:
  /** Default constructor
   *
   * @param statisticsSweep_: object with temperatures and chemical potentials
   * @param bandStructure_: object with the bandstructure
   * @param dimensionality_: dimensionality of space (should be 3)
   */
  BulkTDrift(StatisticsSweep &statisticsSweep_,
             BaseBandStructure &bandStructure_,
             const int &dimensionality_ = 3);
};

/** Object describing the electric-field drift operator of the BTE
 */
class BulkEDrift : public VectorBTE {
public:
  /** Default constructor
   *
   * @param statisticsSweep_: object with temperatures and chemical potentials
   * @param bandStructure_: object with the bandstructure
   * @param dimensionality_: dimensionality of space (should be 3)
   */
  BulkEDrift(StatisticsSweep &statisticsSweep_,
             BaseBandStructure &bandStructure_,
             const int &dimensionality_ = 3);
};

/** Object describing the eigenvector with zero-eigenvalue of the symmetrized
 * scattering matrix Omega (notation of Cepellotti PRX (2016)) associated with
 * the energy conservation.
 */
class Vector0 : public VectorBTE {
public:
  /** Default constructor
   *
   * @param statisticsSweep_: object with temperatures and chemical potentials
   * @param bandStructure_: object with the bandstructure
   * @param specificHeat: object with the values of specific heat
   */
  Vector0(StatisticsSweep &statisticsSweep_, BaseBandStructure &bandStructure_,
          SpecificHeat &specificHeat);
};

#endif
