#include "ph_interaction_3ph.h"

using namespace std;

// TODO * Method to calculate single V_plus/V_minus
complex<double> calculateSingleV(const PhononTriplet &interactingPhonons,		\
			 const Eigen::Tensor<complex,3> &phononEigenvectors,	\
			 const int numTriplets, const Eigen::tensor<double,4> &ifc3Tensor, \
			 const Eigen::Tensor<double,3> &cellPositions,	\
			 Eigen::Tensor<int,2> &displacedAtoms,const Crystal &crystal, \
			 char procType){

  if((procType != '+') || (procType != '-')){
    cout << "Character procType can only be '+' or '-'.";
    exit();
  }
  
  complex<double> phase, V, V0;
  const std::complex<double> I(0, 1);
  int it,ic,iDim,jDim,kDim;

  Eigen::Tensor<complex double>,3 ev1(3,numAtoms,NumBands) ev1 = state1.getEigenvectors();
  Eigen::Tensor<complex double>,3 ev2(3,numAtoms,NumBands) ev1 = state2.getEigenvectors();
  Eigen::Tensor<complex double>,3 ev3(3,numAtoms,NumBands) ev1 = state3.getEigenvectors();

  //phonon branches: s1,s2,s3
  //phonon wave vectors: q1,q2,q3
  
  //Get vector of atomic species types
  Eigen::VectorXi types = crystal.getAtomicSpecies();

  //Get masses of the atom types
  Eigen::VectorXd masses = crystal.getSpeciesMasses();

  for(it = 0; it < numTriplets; it++){// sum over all triplets
    for(ic = 0; ic < 3; ic++){
      cell2Pos(ic) = cellPositions(it,0,ic);
      cell3Pos(ic) = cellPositions(it,1,ic);
    }
    massNorm = sqrt(masses(types(displacedAtoms(it,0)))*	\
		    masses(types(displacedAtoms(it,1)))*	\
		    masses(types(displacedAtoms(it,2))));
    if(procType == '+'){
      phase = exp(I*(q2.dot(cell2Pos)-q3.dot(cell3Pos)))/massNorm;
    }else if(procType == '-'){
      phase = exp(-I*(q2.dot(cell2Pos)+q3.dot(cell3Pos)))/massNorm;
    }
    
    V0 = (0.0,0.0);
    for(iDim = 0; iDim < 3; iDim++){
      for(jDim = 0; jDim < 3; jDim++){
	for(kDim = 0; kDim < 3; kDim++){
	  if(procType == '+'){
	    V0 += ifc3Tensor(it,iDim,jDim,kDim)*	\
	      ev1(iDim,displacedAtoms(it,0),s1)*	\
	      ev2(jDim,displacedAtoms(it,1),s2)*	\
	      conj(ev3(kDim,displacedAtoms(it,2),s3));
	  }else if(procType == '-'){
	    V0 += ifc3Tensor(it,iDim,jDim,kDim)*	\
	      ev1(iDim,displacedAtoms(it,0),s1)*	\
	      conj(ev2(jDim,displacedAtoms(it,1),s2))*	\
	      conj(ev3(kDim,displacedAtoms(it,2),s3));
	  } 
	}
      }
    }
    V0 *= phase;
  }

  return V0*conj(V0);
}

// * Function to calculate all V_plus/V_minus processes for a given IBZ mode
