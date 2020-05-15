#include "context.h"
#include "crystal.h"
#include "interaction_3ph.h"

class IFC3Parser {
 public:
  /** 
   * Parse third-order force constants (IFC3).
   *
   * Method for reading IFC3 tensor from file 
   * in various supported formats.
   * 
   * @param[in] fileName Name of the IFC3 file in disk.
   * @param[in] format IFC3 file format.
   * @param[out] numTriplets Number of triplets of atoms in displaced supercells.
   * @param[out] ifc3Tensor Rank-4 IFC3 tensor.
   * @param[out] cellPositions Cartesian coordinates of the 2nd and 3rd unitcells
   *  for each triplet. The 1st unitcell is taken to be at origin.
   * @param[out] displacedAtoms Index of the displaced atom for every triplet
   *  of unitcells.
   */
	Interaction3Ph parseFromShengBTE(Context & context, Crystal & crystal_);
//	Interaction3Ph parseFromQE(Context & context);
};
