#include "ph_interaction_3ph.h"

int main(){
  //const double amu=1.66053904e-27; //Kg
  
  //------Set up cubic silicon-------//
  //Grid points along the 3 lattice vectors
  int grid[3] = {8,8,8};
  //Total number of wave vectors
  int nq = grid[2]*grid[1]*grid[0];
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
  string ifc3FileName = "./silicon_test_8c/FORCE_CONSTANTS_3RD";
  string format = "shengbte";
  int numTriplets;
  Eigen::Tensor<double,4> ifc3Tensor;
  Eigen::Tensor<double,3> cellPositions;
  Eigen::Tensor<int,2> displacedAtoms;
  IFC3Parser parser;
  parser.parseIFC3(ifc3FileName, format, numTriplets, ifc3Tensor, cellPositions, displacedAtoms);

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
  
  // Read full BZ wave vectors (Cartesian) from file
  Eigen::MatrixXd q(nq,3);
  string qFileName = "./silicon_test_8c/qpoints_cart_full";
  ifstream qfile(qFileName);
  for(int iq = 0; iq < nq; iq++){
    for(int ip = 0; ip < 3; ip++){
      qfile >> q(iq,ip); //1/nm
    }
  }
  qfile.close();

  // Read full BZ angular frequencies in 2pi.THz from file
  Eigen::MatrixXd omega(nq,numBranches);
  ifstream infile("./silicon_test_8c/omega_full");
  for(int iq = 0; iq < nq; iq++){
    for(int ib = 0; ib < numBranches; ib++){
      infile >> omega(iq,ib);
    }
  }

  // Read full BZ eigenvectors from file
  // Note that this is read as outputted by ShengBTE
  // and has to be reshaped to conform to PHOEBE's format
  // of the eigenvector.
  string evFileName = "./silicon_test_8c/evecs_full";
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

  //TEST:
  double T = 300.0; //temperature in K
  int s1 = 2; //LA
  int iq1 = 2; //3rd q vector in both FBZ and IBZ
  
  // Create 1st phonon mode
  PhononMode ph1;
  ph1.s = s1;
  ph1.iq = iq1;

  // Set up tetrahedra
  TetraData tetra;
  cout << "Setting up the tetrahedra...\n";
  formTets(grid,tetra);
  fillTetsEigs(numBranches,omega,tetra);
  cout << "Done!\n";
  
  cout << "Calculating all V- for mode: (" << s1 << "," << iq1 << ")...\n";
  PhInteraction3Ph phInt;
  phInt.calculateAllVminus(grid, ph1, q, ev, numTriplets,ifc3Tensor,	\
			   cellPositions, displacedAtoms, crysInfo);
  cout << "Done!\n";

  cout << "Calculating W+ and W- for mode: (" << s1 << "," << iq1 << ")...\n";
  phInt.calculateAllW(T,grid,ph1,indexMesh,crysInfo,omega,tetra);
  cout << "Done!\n";

  //Read W+/W- from file and form scattering rates
  int plusProcessCount, minusProcessCount;
  
  string counterFileName = "WCounter.iq"+to_string(iq1)+".s"+to_string(s1);
  ifstream counterFile;
  counterFile.open(counterFileName);
  counterFile >> plusProcessCount;
  counterFile >> minusProcessCount;
  counterFile.close();

  cout << plusProcessCount << " " << minusProcessCount << "\n";
  
  int dummy;
  double Wplus,WplusSum,Wminus,WminusSum;

  WplusSum = 0.0;
  string WplusFileName = "Wp2.iq"+to_string(iq1)+".s"+to_string(s1);
  ifstream WplusFile;
  WplusFile.open(WplusFileName);
  for(int l = 0; l < plusProcessCount; l++){
    WplusFile >> dummy >> dummy >> dummy >> dummy >> Wplus;
    //cout << Wplus << "\n";
    WplusSum += Wplus;
  }
  WplusFile.close();
  
  WminusSum = 0.0;
  string WminusFileName = "Wm2.iq"+to_string(iq1)+".s"+to_string(s1);
  ifstream WminusFile;
  WminusFile.open(WminusFileName);
  for(int l = 0; l < minusProcessCount; l++){
    WminusFile >> dummy >> dummy >> dummy >> dummy >> Wminus;
    //cout << Wminus << "\n";
    WminusSum += Wminus;
  }
  WminusFile.close();

  cout << "Scattering rate for this mode is " << (WplusSum+0.5*WminusSum) << " THz \n";
  
  return 0;
}

