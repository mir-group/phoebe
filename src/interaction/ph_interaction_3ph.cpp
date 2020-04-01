#include "ph_interaction_3ph.h"

using namespace std;

// TODO * Method to calculate single V_plus/V_minus
complex<double> calculateSingleV(PhononTriplet &interactingPhonons,		\
			 Eigen::tensor<complex,3> &phononEigenvectors,	\
			 int numTriplets, Eigen::tensor<double,4> &ifc3Tensor, \
			 Eigen::Tensor<double,3> &cellPositions,	\
			 Eigen::Tensor<int,2> &displacedAtoms, char procType){
  complex<double> phase, V, V0;

  for(int i = 0; i < numTriplets; i++){// sum over all triplets
    phase = 
  }
}

// * Function to calculate all V_plus/V_minus processes for a given IBZ mode
