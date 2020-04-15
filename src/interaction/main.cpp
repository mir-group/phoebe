#include "ifc3_parser.h"
#include "ph_interaction_3ph.h"

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
  Eigen::VectorXi types(numAtoms);
  types << 0,0;
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
    for(int ip = 0; ip < 3; ip++){
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
  double re,im;
  complex<double> z;
  for(int iq = 0; iq < nq; iq++){
    for(int ib = 0; ib < numBands; ib++){
      for(int jb = 0; jb < numBands; jb++){
	evfile >> re >> im;
	z = complex<double>(re,im);
	ev(iq,ib,jb) = z;
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
  PhInteraction3Ph phInt;
  double Vp2 = phInt.calculateSingleV(interactingPhonons, q, numTriplets, ifc3Tensor, cellPositions, \
				 displacedAtoms, crysInfo, '+');
  double Vm2 = phInt.calculateSingleV(interactingPhonons, q, numTriplets, ifc3Tensor, cellPositions, \
				 displacedAtoms, crysInfo, '-');
    
  // Compare to result from ShengBTE
  cout << Vp2 << " " << Vm2 << "\n";
  // Test permutation symmetry

  return 0;
}

/*
int main(){
  string ifc3FileName = "./silicon_test/FORCE_CONSTANTS_3RD";
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
}
*/
