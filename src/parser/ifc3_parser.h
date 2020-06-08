#include "interaction_3ph.h"
#include "context.h"
#include "crystal.h"

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
   * @param[out] numTriplets Number of triplets of atoms in displaced
   * supercells.
   */
  static Interaction3Ph parse(Context &context, Crystal &crystal);
  // private:
  static Interaction3Ph parseFromShengBTE(Context &context, Crystal &crystal);
  static Interaction3Ph parseFromQE(Context &context, Crystal &crystal);
};
