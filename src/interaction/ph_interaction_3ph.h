#include <string>
#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;

//----TESTING AREA----//

/**
 * Data structure to hold some crystal information.
 *
 * This is for testing purposes. Eventually this
 * info will be read from Crystal class.
 */
struct CrystalInfo{
  //Vector of atomic species types
  Eigen::VectorXi types;
  //Masses of the atom types
  Eigen::VectorXd masses;  
  //Number of atoms in the primitive cell
  int numAtoms;
  //Number of phonon branches
  int numBranches;
};

/**
 * Data structure to hold three phonon modes' info.
 *
 * A mode \equiv (branch,wave vector).
 */
struct PhononTriplet{
  int s1, s2, s3; //Three branches
  //Indices picking out three wave vectors
  //from the global full Brillouin zone:
  int iq1, iq2, iq3; 
  //The three eigenvectors:
  Eigen::Tensor<complex<double>,3> ev1, ev2, ev3;
};

/**
 * Data structure to hold info of a single phonon mode.
 *
 * A mode \equiv (branch,wave vector).
 */
struct PhononMode{
  int s; //Three branches
  //Index picking out a single wave vectors
  //from the global full Brillouin zone:
  int iq; 
  //Cartesian wave vector
  Eigen::Vector3d q;
};
//--------------------//

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
  /*
  double calculateSingleV(const PhononTriplet &interactingPhonons, \
				   const Eigen::MatrixXd &q, \  
				   const Eigen::tensor<complex,3> &phononEigenvectors, \
				   const int numTriplets, const Eigen::tensor<double,4> &ifc3Tensor, \
				   const Eigen::Tensor<double,3> &cellPositions, \
				   const Eigen::Tensor<int,2> &displacedAtoms,const CrystalInfo &crysInfo, \
				   const char procType);
				   }*/
  double calculateSingleV(const PhononTriplet &interactingPhonons, const Eigen::MatrixXd &q, \  
				 const int numTriplets, const Eigen::Tensor<double,4> &ifc3Tensor, \
				 const Eigen::Tensor<double,3> &cellPositions, \
				 const Eigen::Tensor<int,2> &displacedAtoms,const CrystalInfo &crysInfo, \
				 const char procType);

  void calculateIrredVminus(const int nq, const int *grid, const PhononMode &mode, \
			    const Eigen::MatrixXi &indexMesh, const Eigen::MatrixXd &qFBZ, \
			    const Eigen::Tensor<complex<double>,3> &ev, const int numTriplets, \
			    const Eigen::Tensor<double,4> &ifc3Tensor,	\
			    const Eigen::Tensor<double,3> &cellPositions, \
			    const Eigen::Tensor<int,2> &displacedAtoms,const CrystalInfo &crysInfo);
};

