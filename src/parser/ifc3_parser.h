#include <string>
#include <iostream>
#include <fstream>
#include <unsupported/Eigen/CXX11/Tensor>

using namespace std;

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
  void parseIFC3(string fileName, string format, int numTriplets, \
		 Eigen::Tensor<double,4> &ifc3Tensor, \
		 Eigen::Tensor<double,3> &cellPositions,\
		 Eigen::Tensor<int,2> &displacedAtoms);
};
