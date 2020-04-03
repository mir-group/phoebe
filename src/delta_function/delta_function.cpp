#include "delta_function.h"

using namespace std;

void formTets(const int grid[3], TetraData &tetra){
  // Some internal variables
  int counter,aux,ip1,jp1,kp1,ii,jj,kk;
  
  // Number of grid points
  int n1 = grid[0]; int n2 = grid[1]; int n3 = grid[2];

  // Number of tetrahedra
  tetra.numTetra = 6*n1*n2*n3;

  // Allocate tetrahedron data holders
  tetra.tetrahedra = Eigen::MatrixXi(tetra.numTetra,6);
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

void fillTetsEigs(const int numBands, const Eigen::MatrixXd &energy, TetraData &tetra){
  // Internal variables
  int ik;
  double temp[4];
  size_t size = sizeof(temp)/sizeof(temp[0]);

  // Allocate tetraEigVals
  tetra.tetraEigVals = Eigen::Tensor<double,3>(tetra.numTetra,numBands,4);
  
  for(int it = 0; it < tetra.numTetra; it++){ //over tetrahedra
    for(int ib = 0; ib < numBands; ib++){ //over bands
      for(int iv = 0; iv < 4; iv++){ //over vertices
	// Index of wave vector
	ik = tetra.tetrahedra(it,iv);
	// Fill tetrahedron vertex
	tetra.tetraEigVals(it,ib,iv) = energy(ik,ib);
	temp[iv] = energy(ik,ib); //save for later
      }
      //sort energies in the vertex
      sort(temp,temp+size);
      //refill tetrahedron vertex
      for(int iv = 0; iv < 4; iv++) tetra.tetraEigVals(it,ib,iv) = temp[iv];
    }
  }
}

double fillTetsWeights(const double energy, const int ib, const int iq, const TetraData &tetra){
  // Tetrahedron weight
  double weight = 0.0;
  // Internal variables
  int it,iv,num;
  double e1,e2,e3,e4,e1e,e2e,e3e,e4e,e21,e31,e41,e32,e42,e43;
  bool c1,c2,c3,c4;
  double tmp = 0.0;
  
  //Number of tetrahedra in which the wave vector belongs 
  num = tetra.qToTetCount(iq);

  for(int i = 0; i < num; i++){//over all tetrahedra
    it = tetra.qToTet(iq,i,0); //get index of tetrahedron
    iv = tetra.qToTet(iq,i,1); //get index of vertex

    //Sorted energies at the 4 vertices
    e1 = tetra.tetraEigVals(it,ib,0);
    e2 = tetra.tetraEigVals(it,ib,1);
    e3 = tetra.tetraEigVals(it,ib,2);
    e4 = tetra.tetraEigVals(it,ib,3);
    
    //Define the shorthands
    // Refer to Lambin and Vigneron prb 29.6 (1984): 3430 to understand
    // what these mean.
    e1e = e1-energy;
    e2e = e2-energy;
    e3e = e3-energy;
    e4e = e4-energy;
    e21 = e2-e1;
    e31 = e3-e1;
    e41 = e4-e1;
    e32 = e3-e2;
    e42 = e4-e2;
    e43 = e4-e3;

    //Check the inequalities
    c1 = (e1 <= energy) && (energy <= e2);
    c2 = (e2 <= energy) && (energy <= e3);
    c3 = (e3 <= energy) && (energy <= e4);
    
    if( !((energy < e1) || (energy > e4)) ){
      switch(iv)
	{//switch over 4 vertices
	case 0:
	  if(c1){
	    tmp = (e2e/e21 + e3e/e31 + e4e/e41)*pow(e1e,2)/e41/e31/e21;

	    if(e1 == e2) tmp = 0.0;
	  }else if(c2){
	    tmp = -0.5*(e3e/pow(e31,2)*(e3e*e2e/e42/e32 + e4e*e1e/e41/e42 + e3e*e1e/e32/e41)
			+ e4e/pow(e41,2)*(e4e*e1e/e42/e31 + e4e*e2e/e42/e32 + e3e*e1e/e31/e32));

	    if(e2 == e3){
	      tmp = -0.5*(e4e*e1e/e41/e42 + e1e/e41
			  + e4e/pow(e41,2)*(e4e*e1e/e42/e31 + e4e/e42 + e1e/e31));
	    }
	  }else if(c3){
	    tmp = pow(e4e,3)/pow(e41,2)/e42/e43;

	    if(e3 == e4){
	      tmp = pow(e4e,2)/pow(e41,2)/e42;
	    }
	  }
	case 1:
	  if(c1){
	    tmp = -pow(e1e,3)/pow(e21,2)/e31/e41;

	    if(e1 == e2) tmp = 0.0;
	  }else if(c2){
	    tmp = -0.5*(e3e/pow(e32,2)*(e3e*e2e/e42/e31 + e4e*e2e/e42/e41 + e3e*e1e/e31/e41)
                        + e4e/pow(e42,2)*(e3e*e2e/e32/e31 + e4e*e1e/e41/e31 + e4e*e2e/e32/e41));

	    if(e2 == e3){
	      tmp = -0.5*(e4e/e42/e41
			  + e4e/pow(e42,2)*(e4e*e1e/e41/e31 + 1.0));
	    }
	  }else if(c3){
	    tmp = pow(e4e,3)/e41/pow(e42,2)/e43;

	    if(e3 == e4) tmp = 0.0;
	  }
	case 2:
	  if(c1){
	    tmp = -pow(e1e,3)/e21/pow(e31,2)/e41;

	    if(e1 == e2) tmp = 0.0;
	  }else if(c2){
	    tmp = 0.5*(e2e/pow(e32,2)*(e3e*e2e/e42/e31 + e4e*e2e/e42/e41 + e3e*e1e/e31/e41)
		       + e1e/pow(e31,2)*(e3e*e2e/e42/e32 + e4e*e1e/e41/e42 + e3e*e1e/e32/e41));
                   
	    if(e2 == e3){
	      tmp = 0.5*(e4e/e42/e41 + e1e/e31/e41
			 + e1e/pow(e31,2)*(e4e*e1e/e41/e42 + e1e/e41));
	    }
	  }else if(c3){
	    tmp = pow(e4e,3)/e41/e42/pow(e43,2);

	    if(e3 == e4) tmp = 0.0;
	  }
	case 3:
	  if(c1){
	    tmp = -pow(e1e,3)/e21/e31/pow(e41,2);
	    
	    if(e1 == e2) tmp = 0.0;
	  }else if(c2){
	    tmp = 0.5*(e2e/pow(e42,2)*(e3e*e2e/e32/e31 + e4e*e1e/e41/e31 + e4e*e2e/e32/e41)
		       + e1e/pow(e41,2)*(e4e*e1e/e42/e31 + e4e*e2e/e42/e32 + e3e*e1e/e31/e32));

	    if(e2 == e3){
	      tmp = 0.5*e1e/pow(e41,2)*(e4e*e1e/e42/e31 + e4e/e42 + e1e/e31);
	    }
	  }else if(c3){
	    tmp = -(e3e/e43 + e2e/e42 + e1e/e41)*pow(e4e,2)/e41/e42/e43;
                   
	    if(e3 == e4) tmp = 0.0;
	  }
	} //switch over 4 vertices

      if( (e1 == e2) && (e1 == e3) && (e1 == e4) & (energy == e1) ) tmp = 0.25;

      weight = weight + tmp;
    } //!((energy < e1) || (energy > e4))
  } //over all tetrahedra

  // Zero out extremely small weights
  if(weight < 1.0e-12) weight = 0.0;

  // Normalize by number of tetrahedra
  weight = weight/tetra.numTetra;

  return weight;
}

/*
//A test for the tetrahedron method
int main(){
  //Number of bands
  int numBands = 12;
  //Grid points along the 3 lattice vectors
  int grid[3] = {18,18,18};
  //Total number of wave vectors
  int nq = 18*18*18;
  //Band structure
  Eigen::MatrixXd energy(nq,numBands);
  //Tetrahedron data
  TetraData tetra;

  //Form tetrahedra
  formTets(grid,tetra);
  
  //Read test phonon angular frequency from file
  // I'm abusing notation here as these are not energy but angular frequency...
  ifstream infile("./phonon_dos_test/angfreq_wBAs18c");
  for(int iq = 0; iq < nq; iq++){
    for(int ib = 0; ib < numBands; ib++){
      infile >> energy(iq,ib);
    }
  }

  //Convert to frequency by dividing by 2pi
  energy = energy/2.0/M_PI;

  //cout << energy.row(0) << "\n";
  //cout << energy.row(nq-1) << "\n";

  //Fill tetrahedra with energies
  fillTetsEigs(numBands,energy,tetra);

  //Array of uniform frequency to sample
  int numFreq = 201;
  double freq[numFreq];
  double del = 0.125; //frequency spacing, THz
  for(int i = 0; i < numFreq; i++){
    freq[i] = i*del;
  }
  
  //Calculate phonon density of states (DOS) [1/THz]
  // I'll implement a parallel version of this as a separate method later...
  double weight; // tetrahedron weight
  double dos[numFreq] = {0.0}; // phonon DOS initialized to zero
  
  for(int i = 0; i < numFreq; i++){
    for(int iq = 0; iq < nq; iq++){
      for(int ib = 0; ib < numBands; ib++){
	weight = fillTetsWeights(freq[i],ib,iq,tetra);
	dos[i] = dos[i] + weight;
      }
    }
  }
  
  //Save phonon DOS to file
  ofstream outfile("./phonon_dos_test/phonon_dos");
  for(int i = 0; i < numFreq; i++){
    outfile << freq[i] << "\t" << dos[i] << "\n";
  }
  
  return 0;
}
*/
