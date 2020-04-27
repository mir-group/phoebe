#include "ifc3_parser.h"
#include "ph_interaction_3ph.h"

int main(){
  //const double amu=1.66053904e-27; //Kg
  
  //------Set up cubic silicon-------//
  //Grid points along the 3 lattice vectors
  int grid[3] = {20,20,20};
  //Total number of wave vectors
  int nq = 20*20*20;
  //Number of atoms
  int numAtoms = 2;
  //Number of bands
  int numBranches = 3*numAtoms;
  //Number of types of atoms
  int numTypes = 1;
  //Types of atoms
  Eigen::VectorXi types(numAtoms);
  types << 0,0;
  //Masses of the types of atoms
  Eigen::VectorXd masses(numTypes);
  masses << 28.085509989600006; //amu
  //fill CrystalInfo
  CrystalInfo crysInfo;
  crysInfo.numAtoms = numAtoms;
  crysInfo.numBranches = numBranches;
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
  Eigen::Tensor<complex<double>,3> ev(nq,numBranches,numBranches);
  ifstream evfile(evFileName);
  double re,im;
  complex<double> z;
  for(int iq = 0; iq < nq; iq++){
    for(int ib = 0; ib < numBranches; ib++){
      for(int jb = 0; jb < numBranches; jb++){
	evfile >> re >> im;
	z = complex<double>(re,im);
	ev(iq,ib,jb) = z;
      }
    }
  }
  evfile.close();
  
  // For a single first phonon mode, calculate all V- processes:
  int s1 = 0; //TA1
  int iq1 = 0; //Gamma point in both FBZ and IBZ
  //int iq1 = 1; 

  // Create 1st phonon mode
  PhononMode ph1;
  ph1.s = s1;
  ph1.iq = iq1;

  // Create indexMesh
  // For testing purposes this is in ShengBTE ordering
  Eigen::MatrixXi indexMesh(nq,3);
  int count = 0;
  for(int k = 0; k < grid[2]; k++){
    for(int j = 0; j < grid[1]; j++){
      for(int i = 0; i < grid[0]; i++){
	indexMesh(count,0) = i;
	indexMesh(count,1) = j;
	indexMesh(count,2) = k;
	count++;
      }
    }
  }

  cout << "Calculating irreducible set of V- for mode: (" << s1 << "," << iq1 << ")\n";
  PhInteraction3Ph phInt;
  phInt.calculateIrredVminus(nq, grid, ph1, indexMesh, q, ev, numTriplets, \
			     ifc3Tensor, cellPositions, displacedAtoms, crysInfo);
  			     
  /*
  // Form a triplet to test vertex calculator
  int s1 = 0;
  int s2 = 0;
  int s3 = 0;
  int iq1 = 0;
  int iq2 = 0;
  int iq3 = 0;
  cout << "triplet: (" << s1 << "," << iq1 << ") (" << s2 << "," << iq2 << ") (" << s3 << "," << iq3 << "):\n";
  
  // Reshape the eigenvectors read from file
  Eigen::Tensor<complex<double>,3> ev1(3,numAtoms,numBranches);
  Eigen::Tensor<complex<double>,3> ev2(3,numAtoms,numBranches);
  Eigen::Tensor<complex<double>,3> ev3(3,numAtoms,numBranches);

  for(int idim = 0; idim < 3; idim++){
    for(int iat = 0; iat < numAtoms; iat++){
      for(int ib = 0; ib < numBranches; ib++){
	ev1(idim,iat,ib) = ev(iq1,idim+3*iat,ib);
	ev2(idim,iat,ib) = ev(iq2,idim+3*iat,ib);
	ev3(idim,iat,ib) = ev(iq3,idim+3*iat,ib);
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
  //double Vp2 = phInt.calculateSingleV(interactingPhonons, q, numTriplets, ifc3Tensor, cellPositions, \
				 displacedAtoms, crysInfo, '+');
  double Vm2 = phInt.calculateSingleV(interactingPhonons, q, numTriplets, ifc3Tensor, cellPositions, \
				 displacedAtoms, crysInfo, '-');
 
  //cout << "|V+|^2 = " << Vp2 << " |V-|^2 = " << Vm2 << "\n";
  cout << "..............\n";
  */
  /*
  // TODO: Compare to result from ShengBTE
  
  // Test permutation symmetry: V^{-}(lam1,lam2,lam3) = V^{-}(lam1,lam3,lam2)
  cout << "Test exchange symmetry of - process.\n";
  cout << "Is V^{-}(lam1,lam2,lam3) = V^{-}(lam1,lam3,lam2)?\n";
  // Form a triplet to test vertex calculator 
  s1 = 0; //TA1
  s3 = 1; //TA2
  s2 = 2; //LA
  iq1 = 6; //Gamma
  iq3 = 19;
  iq2 = 13;
  cout << "triplet: (" << s1 << "," << iq1 << ") (" << s2 << "," << iq2 << ") (" << s3 << "," << iq3 << "):\n";
  
  // Reshape the eigenvectors read from file
  for(int idim = 0; idim < 3; idim++){
    for(int iat = 0; iat < numAtoms; iat++){
      for(int ib = 0; ib < numBranches; ib++){
	ev1(idim,iat,ib) = ev(iq1,ib,iat*numAtoms+idim);
	ev2(idim,iat,ib) = ev(iq2,ib,iat*numAtoms+idim);
	ev3(idim,iat,ib) = ev(iq3,ib,iat*numAtoms+idim);
      }
    }
  }
  
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
  Vp2 = phInt.calculateSingleV(interactingPhonons, q, numTriplets, ifc3Tensor, cellPositions, \
				 displacedAtoms, crysInfo, '+');
  Vm2 = phInt.calculateSingleV(interactingPhonons, q, numTriplets, ifc3Tensor, cellPositions, \
				 displacedAtoms, crysInfo, '-');

  cout << "|V+|^2 = " << Vp2 << " |V-|^2 = " << Vm2 << "\n";
  cout << "..............\n";
  // TODO: Compare to result from ShengBTE
  /////////////////////////////////////

  // Test time-reversal symmetry:
  // 1. V^{-}(lam1,-lam2,lam3) = V^{+}(lam1,lam2,lam3)
  // 1. V^{+}(lam1,-lam2,lam3) = V^{-}(lam1,lam2,lam3)

  cout << "Test time-reversal mapping of + and - processes.\n";
  cout << "1) Is V^{-}(lam1,-lam2,lam3) = V^{+}(lam1,lam2,lam3)?\n";
  cout << "2) Is V^{+}(lam1,-lam2,lam3) = V^{-}(lam1,lam2,lam3)?\n";
  // Form a triplet to test vertex calculator
  s1 = 0; //TA1
  s2 = 1; //TA2
  s3 = 2; //LA
  iq1 = 6; //Gamma
  iq2 = 19;
  //Demux assuming shengBTE ordering
  int N = grid[0];
  int qx = iq2%N;
  int qy = (iq2/N)%N;
  int qz = iq2/N/N;
  //Negate and Umklapp
  qx = (-qx+N)%N;
  qy = (-qy+N)%N;
  qz = (-qz+N)%N;
  //Mux assuming shengBTE ordering
  int TR_iq2 = (qz*N + qy)*N + qx;
  iq2 = TR_iq2;
  iq3 = 13;
  
  cout << "triplet: (" << s1 << "," << iq1 << ") (" << s2 << "," << iq2 << ") (" << s3 << "," << iq3 << "):\n";
  
  // Reshape the eigenvectors read from file
  for(int idim = 0; idim < 3; idim++){
    for(int iat = 0; iat < numAtoms; iat++){
      for(int ib = 0; ib < numBranches; ib++){
	ev1(idim,iat,ib) = ev(iq1,ib,iat*numAtoms+idim);
	ev2(idim,iat,ib) = ev(iq2,ib,iat*numAtoms+idim);
	ev3(idim,iat,ib) = ev(iq3,ib,iat*numAtoms+idim);
      }
    }
  }
  
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
  Vp2 = phInt.calculateSingleV(interactingPhonons, q, numTriplets, ifc3Tensor, cellPositions, \
				 displacedAtoms, crysInfo, '+');
  Vm2 = phInt.calculateSingleV(interactingPhonons, q, numTriplets, ifc3Tensor, cellPositions, \
				 displacedAtoms, crysInfo, '-');
    
  cout << "|V+|^2 = " << Vp2 << " |V-|^2 = " << Vm2 << "\n";
  cout << "..............\n";
  // TODO: Compare to result from ShengBTE
  /////////////////////////////////////
  
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
  */

  return 0;
}
