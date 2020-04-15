#include "ifc3_parser.h"

using namespace std;

void IFC3Parser::parseIFC3(string fileName, string format, int &numTriplets, \
	       Eigen::Tensor<double,4> &ifc3Tensor,		\
	       Eigen::Tensor<double,3> &cellPositions, \
	       Eigen::Tensor<int,2> &displacedAtoms){
  int tripletCount; // triplet counter 
  int tmp[3]; // temporary int holder
  
  // Open IFC3 file
  ifstream infile(fileName);
  
  // Check format of input IFC3 file
  if(format.compare("shengbte") == 0){
    // Number of triplets
    infile >> numTriplets;

    // Allocate readables
    ifc3Tensor = Eigen::Tensor<double,4>(numTriplets,3,3,3);
    cellPositions = Eigen::Tensor<double,3>(numTriplets,2,3);
    displacedAtoms = Eigen::Tensor<int,2>(numTriplets,3);
    
    for(int i = 0; i < numTriplets; i++){// loop over all triplets
      // Triplet counter
      infile >> tripletCount;
      
      // Read position of 2nd cell
      infile >> cellPositions(i,0,0) >> cellPositions(i,0,1) >> cellPositions(i,0,2);

      // Read position of 3rd cell
      infile >> cellPositions(i,1,0) >> cellPositions(i,1,1) >> cellPositions(i,1,2);
      
      // Read triplet atom indices
      infile >> displacedAtoms(i,0) >> displacedAtoms(i,1) >> displacedAtoms(i,2);
      for(int a  = 0; a < 3; a++)
	displacedAtoms(i,a) = displacedAtoms(i,a) - 1; // go to zero based indexing

      // Read the 3x3x3 force constants tensor
      for(int a = 0; a < 3; a++){
	for(int b = 0; b < 3; b++){
	  for(int c = 0; c < 3; c++){
	    infile >> tmp[0] >> tmp[1] >> tmp[2] >> ifc3Tensor(i,a,b,c);
	  }
	}
      }
    }
  }else{
    cout << "parseIFC3: Only the shengbte format is currently supported.\n";
    exit(-1);
  }
  
  // Close IFC3 file
  infile.close();

  //TODO Round cellPositions to the nearest lattice vectors
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

  parseIFC3(ifc3FileName, format, numTriplets, ifc3Tensor, cellPositions, displacedAtoms);

  int iTrip = 31;
  cout << cellPositions(iTrip,0,0) << " " << cellPositions(iTrip,0,1) << " " << cellPositions(iTrip,0,2) << " "  << endl;
  cout << cellPositions(iTrip,1,0) << " " << cellPositions(iTrip,1,1) << " " << cellPositions(iTrip,1,2) << " "  << endl;
  cout << displacedAtoms(iTrip,0) << " " << displacedAtoms(iTrip,1) << " " << displacedAtoms(iTrip,2) << " " << endl;
  cout << ifc3Tensor(iTrip,1,2,1) << endl;
  
  return 0;
  }*/
