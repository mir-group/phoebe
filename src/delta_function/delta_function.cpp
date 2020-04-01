#include "delta_function.h"

using namespace std;

void formTets(int grid[3], TetraData &tetra){
  //Some internal variables
  int counter,aux,ip1,jp1,kp1,ii,jj,kk;
  
  // Number of grid points
  int n1 = grid[0]; int n2 = grid[1]; int n3 = grid[2];

  // Number of tetrahedra
  int numTetra = 6*n1*n2*n3;

  // Allocate tetrahedron data holders 
  tetra.tetrahedra = Eigen::MatrixXi(numTetra,6);
  tetra.qToTetCount = Eigen::VectorXi::Zero(n1*n2*n3);
  tetra.qToTet = Eigen::Tensor<int,3>(n1*n2*n3,24,2);
    
  // Label the vertices of each tetrahedron in a subcell
  Eigen::MatrixXi verticesLabels(6,4);
  verticesLabels <<
    0,1,2,5,
    0,2,4,5,
    2,4,5,6,
    2,5,6,7,
    2,3,5,7,
    1,2,3,5;

  // 8 corners of a subcell
  Eigen::MatrixXi subcellCorners(8,3);
  counter = 0;
  //tetra.qToTetCount = ;
  for(int i = 0; i < n1; i++){
    ip1 = (i+1)%n1;
    for(int j = 0; j < n2; j++){
      jp1 = (j+1)%n2;
      for(int k = 0; k < n3; k++){
	kp1 = (k+1)%n3;
	subcellCorners <<
	  i,j,k,
	  ip1,j,k,
	  i,jp1,k,
	  ip1,jp1,k,
	  i,j,kp1,
	  ip1,j,kp1,
	  i,jp1,kp1,
	  ip1,jp1,kp1;

	//cout << i << j << k << i*n2*n3 + j*n3 + k << "\n";
	
	for(int it = 0; it < 6; it++){ //over 6 tetrahedra
	  for(int iv = 0; iv < 4; iv++){ //over 4 vertices
	    // Grab a label
	    aux = verticesLabels(it,iv);
	    // Grab a corner of subcell
	    ii = subcellCorners(aux,0);
	    jj = subcellCorners(aux,1);
	    kk = subcellCorners(aux,2);
	    // Get muxed index of corner
	    aux = ii*n2*n3 + jj*n3 + kk;
	    // Save corner as a tetrahedron vertex
	    tetra.tetrahedra(counter,iv) = aux;
	    // Save mapping of a wave vector index
	    // to the ordered pair (tetrahedron,vertex)
	    tetra.qToTetCount(aux) = tetra.qToTetCount(aux) + 1;
	    tetra.qToTet(aux,tetra.qToTetCount(aux)-1,0) = counter;
	    tetra.qToTet(aux,tetra.qToTetCount(aux)-1,1) = iv;
	  }
	  counter++;
	}
      }
    }
  }  
}

int main(){
  
  int grid[3] = {10,10,10};

  TetraData tetra;
  formTets(grid,tetra);

  cout << tetra.qToTetCount(0) << "\n";
  cout << tetra.qToTet(1,3,0) << "\n";
  cout << tetra.qToTet(1,3,1) << "\n";
  
  return 0;
}
