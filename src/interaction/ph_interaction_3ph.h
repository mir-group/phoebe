#include <string>
#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <unsupported/Eigen/CX11/Tensor>
#include <Eigen/Core>

using namespace std;

// * Function to calculate all V_plus/V_minus processes for a given IBZ mode
// * Function to calculate single V_plus/V_minus
class PhInteraction3Ph{
 public:
  /**
   * Calculate single three-phonon matrix element (V^{+}/V^{-1}).
   * 
   * Method for calculating the matrix element of a single, plus
   * or minus, three-phonon scattering process.
   *
   * @param[in] interactingPhonons phononTriplet data structure 
   *  containing the three phonon modes (in the full Brillouin zone) 
   *  taking part in an interaction event.
   * @param[in] phononEigenvectors Eigenvectors of all phonons in the
   *  full Brillouin zone.
   * @param[in] numTriplets Number of triplets considered in the supercell
   *  for the third order force constants calculation.
   * @param[in] ifc3Tensor Third order force constants.
   * @param[in] cellPositions Cartesian coordinates of the 2nd and 3rd unitcells.
   * @param[in] displacedAtoms Index of the displaced atom for every triplet
   *  of unitcells.
   * @param[in] crystal Contains of crystal information
   * @param[in] procType Character '+' or '-' to choose the +/- type process.
   * @return The complex matrix element.
   */
  complex<double> calculateSingleV(const PhononTriplet &interactingPhonons, \
			   const Eigen::tensor<complex,3> &phononEigenvectors, \
			   const int numTriplets, const Eigen::tensor<double,4> &ifc3Tensor, \
			   const Eigen::Tensor<double,3> &cellPositions,	\
			   const Eigen::Tensor<int,2> &displacedAtoms,const Crystal &crystal,\
			   const char procType);

  /**
   * Data structure to hold three phonon modes.
   *
   * A mode \equiv (branch,wave vector).
   */
  /*
  struct PhononTriplet{
    int s1, s2, s3; //Three branches
    //Indices picking out three wave vectors
    //from the global full Brillouin zone:
    int iq1, iq2, iq3; 
  };
  */
}
