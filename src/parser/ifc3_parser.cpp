#include "ifc3_parser.h"

using namespace std;

void parseIFC3(string fileName, string format, int numTriplets, \
	       vector<vector<vector<vector<double>>>> &ifc3Tensor,	\
	       vector<vector<vector<double>>> &cellPositions,\
	       vector<vector<int>> &displacedAtoms){
  string line; // line read from file
  int tripletCount; // triplet counter 
  int tmp[3]; // temporary int holder
  
  // Open IFC3 file
  ifstream infile(fileName);
  
  // Check format of input IFC3 file
  if(format.compare("shengbte") == 0){
    
    // Number of triplets
    infile >> numTriplets;

    // Allocate readables
    ifc3Tensor.resize(numTriplets,vector<vector<vector<double>>>(3));
    for(int i = 0; i < numTriplets; i++)
      for(int j = 0; j < 3; j++)
	ifc3Tensor[i][j].resize(3);
    for(int i = 0; i < numTriplets; i++)
      for(int j = 0; j < 3; j++)
	for(int k = 0; k < 3; k++)
	  ifc3Tensor[i][j][k].resize(3);

    cellPositions.resize(numTriplets,vector<vector<double>>(2));
    for(int i = 0; i < numTriplets; i++)
      for(int j = 0; j < 2; j++)
	cellPositions[i][j].resize(3);

    displacedAtoms.resize(numTriplets,vector<int>(3));
    
    for(int i = 0; i < numTriplets; i++){// loop over all triplets
      // Triplet counter
      infile >> tripletCount;
      
      // Read position of 2nd cell
      infile >> cellPositions[i][0][0] >> cellPositions[i][0][1] >> cellPositions[i][0][2];

      // Read position of 3rd cell
      infile >> cellPositions[i][1][0] >> cellPositions[i][1][1] >> cellPositions[i][1][2];
      
      // Read triplet atom indices
      infile >> displacedAtoms[i][0] >> displacedAtoms[i][1] >> displacedAtoms[i][2];
      for(int a  = 0; a < 3; a++)
	displacedAtoms[i][a] = displacedAtoms[i][a] - 1; // go to zero based indexing

      // Read the 3x3x3 force constants tensor
      for(int a = 0; a < 3; a++){
	for(int b = 0; b < 3; b++){
	  for(int c = 0; c < 3; c++){
	    infile >> tmp[0] >> tmp[1] >> tmp[2] >> ifc3Tensor[i][a][b][c];
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
  string ifc3FileName = "FORCE_CONSTANTS_3RD";
  string format = "shengbte";
  int numTriplets;
  vector<vector<vector<vector<double>>>> ifc3Tensor;
  vector<vector<vector<double>>> cellPositions;
  vector<vector<int>> displacedAtoms;

  parseIFC3(ifc3FileName, format, numTriplets, ifc3Tensor, cellPositions, displacedAtoms);

  cout << cellPositions[615][0][1];
  cout << displacedAtoms[615][0];
  cout << ifc3Tensor[615][2][2][1];
  
  return 0;
}
*/

