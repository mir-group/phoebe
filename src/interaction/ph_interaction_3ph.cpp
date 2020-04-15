#include "ph_interaction_3ph.h"
#include "ifc3_parser.h"

using namespace std;

// TODO * Method to calculate single V_plus/V_minus
/*
complex<double> PhInteraction3Ph::calculateSingleV(const PhononTriplet &interactingPhonons,		\
				 const Eigen::MatrixXd &q,		\
				 const Eigen::Tensor<complex,3> &phononEigenvectors, \
				 const int numTriplets, const Eigen::tensor<double,4> &ifc3Tensor, \
				 const Eigen::Tensor<double,3> &cellPositions, \
				 Eigen::Tensor<int,2> &displacedAtoms,const CrystalInfo &crysInfo, \
				 char procType){*/
double PhInteraction3Ph::calculateSingleV(const PhononTriplet &interactingPhonons, \
				 const Eigen::MatrixXd &q,		\
				 const int numTriplets, const Eigen::Tensor<double,4> &ifc3Tensor, \
				 const Eigen::Tensor<double,3> &cellPositions, \
				 const Eigen::Tensor<int,2> &displacedAtoms,const CrystalInfo &crysInfo, \
				 char procType){
  
  if((procType != '+') && (procType != '-')){
    cout << "Character procType can only be '+' or '-'.";
    exit(-1);
  }
  
  complex<double> phase, V, V0;
  const std::complex<double> I(0, 1);
  int s1,s2,s3,it,ic,iDim,jDim,kDim,numAtoms,numBands;
  double massNorm;
  Eigen::Vector3d q1,q2,q3;
  Eigen::Vector3d cell2Pos, cell3Pos;
  
  //Get vector of atomic species types
  //Eigen::VectorXi types = crystal.getAtomicSpecies();
  //Get masses of the atom types
  //Eigen::VectorXd masses = crystal.getSpeciesMasses();
  //numAtoms = crystal.getNumAtoms();
  //numBands = 3*numAtoms; // number of phonon branches in 3d
  //Eigen::Tensor<complex double,3> ev1(3,numAtoms,numBands) = state1.getEigenvectors();
  //Eigen::Tensor<complex double,3> ev2(3,numAtoms,numBands) = state2.getEigenvectors();
  //Eigen::Tensor<complex double,3> ev3(3,numAtoms,numBands) = state3.getEigenvectors();
  
  //For now I'm grabbing the following info from the PhononTriplet and CrystalInfo structs
  //defined below in the testing area. Later on I will connect to the Crystal and State classes. 
  Eigen::Tensor<complex <double>,3> ev1(3,crysInfo.numAtoms,crysInfo.numBands); 
  Eigen::Tensor<complex <double>,3> ev2(3,crysInfo.numAtoms,crysInfo.numBands);
  Eigen::Tensor<complex <double>,3> ev3(3,crysInfo.numAtoms,crysInfo.numBands);
  
  ev1 = interactingPhonons.ev1;
  ev2 = interactingPhonons.ev2;
  ev3 = interactingPhonons.ev3;

  //phonon branches: s1,s2,s3
  s1 = interactingPhonons.s1;
  s2 = interactingPhonons.s2;
  s3 = interactingPhonons.s3;
  
  //Cartesian phonon wave vectors: q1,q2,q3
  q1 = q.row(interactingPhonons.iq1);
  q2 = q.row(interactingPhonons.iq2);
  q3 = q.row(interactingPhonons.iq3);
  
  for(it = 0; it < numTriplets; it++){// sum over all triplets
    for(ic = 0; ic < 3; ic++){
      cell2Pos(ic) = cellPositions(it,0,ic);
      cell3Pos(ic) = cellPositions(it,1,ic);
    }

    massNorm = sqrt(crysInfo.masses(crysInfo.types(displacedAtoms(it,0)))*	\
		    crysInfo.masses(crysInfo.types(displacedAtoms(it,1)))*	\
		    crysInfo.masses(crysInfo.types(displacedAtoms(it,2))));
    
    if(procType == '+'){
      phase = exp(I*(q2.dot(cell2Pos)-q3.dot(cell3Pos)))/massNorm;
    }else if(procType == '-'){
      phase = exp(-I*(q2.dot(cell2Pos)+q3.dot(cell3Pos)))/massNorm;
    }
         
    V0 = complex<double>(0.0,0.0);
    
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

  return abs(V0*conj(V0)); //in units of 
}

// * Function to calculate all V_plus/V_minus processes for a given IBZ mode

/*
//----TESTING AREA----//
int main(){
  //const double amu=1.66053904e-27; //Kg
  
  // TODO:
  //------Set up cubic silicon-------//
  //Grid points along the 3 lattice vectors
  int grid[3] = {20,20,20};
  //Total number of wave vectors
  int nq = 20*20*20;
  //Number of atoms
  int numAtoms = 2;
  //Number of bands
  int numBands = 3*numAtoms;
  //Number of types of atoms
  int numTypes = 1;
  //Types of atoms
  Eigen::VectorXi types(numTypes);
  types << 0;
  //Masses of the types of atoms
  Eigen::VectorXd masses(numTypes);
  masses << 27.976928; //amu
  //fill CrystalInfo
  CrystalInfo crysInfo;
  crysInfo.numAtoms = numAtoms;
  crysInfo.numBands = numBands;
  crysInfo.types = types;
  crysInfo.masses = masses;
  //------------end setup-----------//
  
  // Read ifc3 from file
  string ifc3FileName = "./silicon_test/FORCE_CONSTANTS_3RD";
  string format = "shengbte";
  int numTriplets;
  Eigen::Tensor<double,4> ifc3Tensor;
  Eigen::Tensor<double,3> cellPositions;
  Eigen::Tensor<int,2> displacedAtoms;

  IFC3Parser parser;
  parser.parseIFC3(ifc3FileName, format, numTriplets, ifc3Tensor, cellPositions, displacedAtoms);

  // Read full BZ wave vectors (Cartesian) from file
  string qFileName = "./silicon_test/qpoints_full_cart";
  Eigen::MatrixXd q(nq,3);
  ifstream qfile(qFileName);
  for(int iq = 0; iq < nq; iq++){
    for(int ip = 0; ip < numBands; ip++){
      qfile >> q(iq,ip);
    }
  }
  qfile.close();

  // Read full BZ eigenvectors from file
  // Note that this is read as outputted by ShengBTE
  // and has to be reshaped to conform to PHOEBE's format
  // of the eigenvector.
  string evFileName = "./silicon_test/evecs_full";
  Eigen::Tensor<complex<double>,3> ev(nq,numBands,numBands);
  ifstream evfile(evFileName);
  for(int iq = 0; iq < nq; iq++){
    for(int ib = 0; ib < numBands; ib++){
      for(int jb = 0; jb < numBands; jb++){
	evfile >> ev(iq,ib,jb);
      }
    }
  }
  evfile.close();

  // Form a triplet to test vertex calculator
  int s1 = 0; //TA1
  int s2 = 1; //TA2
  int s3 = 2; //LA
  int iq1 = 0; //Gamma
  int iq2 = 5;
  int iq3 = 9;
  // Reshape the eigenvectors read from file
  Eigen::Tensor<complex<double>,3> ev1(3,numAtoms,numBands);
  Eigen::Tensor<complex<double>,3> ev2(3,numAtoms,numBands);
  Eigen::Tensor<complex<double>,3> ev3(3,numAtoms,numBands);
  for(int idim = 0; idim < 3; idim++){
    for(int iat = 0; iat < numAtoms; iat++){
      for(int ib = 0; ib < numBands; ib++){
	ev1(idim,iat,ib) = ev(iq1,ib,iat*numAtoms+idim);
	ev2(idim,iat,ib) = ev(iq2,ib,iat*numAtoms+idim);
	ev3(idim,iat,ib) = ev(iq3,ib,iat*numAtoms+idim);
      }
    }
  }
  PhononTriplet interactingPhonons;
  interactingPhonons.s1 = s1;
  interactingPhonons.s2 = s2;
  interactingPhonons.s3 = s3;
  interactingPhonons.iq1 = iq1;
  interactingPhonons.iq2 = iq2;
  interactingPhonons.iq3 = iq3;
  interactingPhonons.ev1 = ev1;
  interactingPhonons.ev2 = ev2;
  interactingPhonons.ev3 = ev3;

  // Calculate single process vertex
  complex<double> Vp2 = calculateSingleV(interactingPhonons, q, numTriplets, ifc3Tensor, cellPositions, \
				 displacedAtoms, crysInfo, '+');
  complex<double> Vm2 = calculateSingleV(interactingPhonons, q, numTriplets, ifc3Tensor, cellPositions, \
				 displacedAtoms, crysInfo, '-');

  // Compare to result from ShengBTE
  // Test permutation symmetry
}*/

/*
int main(){
  string ifc3FileName = "FORCE_CONSTANTS_3RD";
  string format = "shengbte";
  int numTriplets;
  //vector<vector<vector<vector<double>>>> ifc3Tensor;
  //vector<vector<vector<double>>> cellPositions;
  //vector<vector<int>> displacedAtoms;

  Eigen::Tensor<double,4> ifc3Tensor;
  Eigen::Tensor<double,3> cellPositions;
  Eigen::Tensor<int,2> displacedAtoms;

  IFC3Parser parser;
  parser.parseIFC3(ifc3FileName, format, numTriplets, ifc3Tensor, cellPositions, displacedAtoms);
  
  int iTrip = 31;
  cout << cellPositions(iTrip,0,0) << " " << cellPositions(iTrip,0,1) << " " << cellPositions(iTrip,0,2) << " "  << endl;
  cout << cellPositions(iTrip,1,0) << " " << cellPositions(iTrip,1,1) << " " << cellPositions(iTrip,1,2) << " "  << endl;
  cout << displacedAtoms(iTrip,0) << " " << displacedAtoms(iTrip,1) << " " << displacedAtoms(iTrip,2) << " " << endl;
  cout << ifc3Tensor(iTrip,1,2,1) << endl;
  
  return 0;
  }*/
