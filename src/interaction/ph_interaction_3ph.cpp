#include "ph_interaction_3ph.h"

using namespace std;

//Some constants.
//These should be globally available.
const double kb = 1.380648813e-23; // J/K
const double hbar = 1.05457172647e-22; // J/THz

//Bose distribution
//(there should be a globally accessible version of this)
double bose(const double omega, const double T){
  return 1.0/(exp(hbar*omega/(kb*T)) - 1.0);
}

//Time reversal operation on a wave vector
//identified by a muxed index.
//This too should be an utility method that is
//globally available.
int timeReverse(const int iq, const int grid[3]){
  int ix,iy,iz;
  //!!WARNING: For testing purposes using ShengBTE ordering!!!
  // demux index
  ix = iq%grid[0];
  iy = (iq/grid[0])%grid[1];
  iz = iq/grid[0]/grid[1];
  // take negative and Umklapp
  ix = (-ix+grid[0])%grid[0];
  iy = (-iy+grid[1])%grid[1];
  iz = (-iz+grid[2])%grid[2];
  // return muxed index
  return (iz*grid[1] + iy)*grid[0] + ix;
}

double PhInteraction3Ph::calculateSingleV(const PhononTriplet &interactingPhonons,\
				 const Eigen::MatrixXd &q,\
				 const int numTriplets, const Eigen::Tensor<double,4> &ifc3Tensor,\
				 const Eigen::Tensor<double,3> &cellPositions,\
				 const Eigen::Tensor<int,2> &displacedAtoms,const CrystalInfo &crysInfo,\
				 const char procType){
  
  if((procType != '+') && (procType != '-')){
    cout << "Character procType can only be '+' or '-'.";
    exit(-1);
  }
  
  complex<double> phase, V0, V;
  const std::complex<double> I(0, 1);
  int s1,s2,s3,it,ic,iDim,jDim,kDim,numAtoms,numBranches;
  double massNorm;
  Eigen::Vector3d q1,q2,q3;
  Eigen::Vector3d cell2Pos, cell3Pos;
  
  //Get vector of atomic species types
  //Eigen::VectorXi types = crystal.getAtomicSpecies();
  //Get masses of the atom types
  //Eigen::VectorXd masses = crystal.getSpeciesMasses();
  //numAtoms = crystal.getNumAtoms();
  //numBranches = 3*numAtoms; // number of phonon branches in 3d
  //Eigen::Tensor<complex double,3> ev1(3,numAtoms,numBranches) = state1.getEigenvectors();
  //Eigen::Tensor<complex double,3> ev2(3,numAtoms,numBranches) = state2.getEigenvectors();
  //Eigen::Tensor<complex double,3> ev3(3,numAtoms,numBranches) = state3.getEigenvectors();
  
  //For now I'm grabbing the following info from the PhononTriplet and CrystalInfo structs
  //defined below in the testing area. Later on I will connect to the Crystal and State classes. 
  Eigen::Tensor<complex <double>,3> ev1(3,crysInfo.numAtoms,crysInfo.numBranches); 
  Eigen::Tensor<complex <double>,3> ev2(3,crysInfo.numAtoms,crysInfo.numBranches);
  Eigen::Tensor<complex <double>,3> ev3(3,crysInfo.numAtoms,crysInfo.numBranches);
  
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

  V = complex<double>(0.0,0.0);
  for(it = 0; it < numTriplets; it++){// sum over all triplets
    for(ic = 0; ic < 3; ic++){
      cell2Pos(ic) = cellPositions(it,0,ic);
      cell3Pos(ic) = cellPositions(it,1,ic);
    }

    massNorm = sqrt(crysInfo.masses(crysInfo.types(displacedAtoms(it,0)))*\
		    crysInfo.masses(crysInfo.types(displacedAtoms(it,1)))*\
		    crysInfo.masses(crysInfo.types(displacedAtoms(it,2))));

    //Recall that the first primitive cell in the triplet is restricted to
    //the origin, so the phase for that cell is unity.
    if(procType == '+'){
      phase = exp(I*(q2.dot(cell2Pos)-q3.dot(cell3Pos)))/massNorm;
    }else if(procType == '-'){
      phase = exp(-I*(q2.dot(cell2Pos)+q3.dot(cell3Pos)))/massNorm;
    }
         
    V0 = complex<double>(0.0,0.0);
    for(kDim = 0; kDim < 3; kDim++){
      for(jDim = 0; jDim < 3; jDim++){
	for(iDim = 0; iDim < 3; iDim++){
	  if(procType == '+'){
	    V0 += ifc3Tensor(it,iDim,jDim,kDim)*\
	      ev1(iDim,displacedAtoms(it,0),s1)*\
	      ev2(jDim,displacedAtoms(it,1),s2)*\
	      conj(ev3(kDim,displacedAtoms(it,2),s3));	    
	  }else if(procType == '-'){
	    V0 += ifc3Tensor(it,iDim,jDim,kDim)*\
	      ev1(iDim,displacedAtoms(it,0),s1)*\
	      conj(ev2(jDim,displacedAtoms(it,1),s2))*\
	      conj(ev3(kDim,displacedAtoms(it,2),s3));
	  }
	}
      }
    }
    V += V0*phase;
  }
  return pow(abs(V),2); //in units of ?
}

// Function to calculate the full set of V_minus processes for a given IBZ mode
void PhInteraction3Ph::calculateAllVminus(const int grid[3], const PhononMode &mode,\
					  const Eigen::MatrixXd &qFBZ,	\
					  const Eigen::Tensor<complex<double>,3> &ev,const int numTriplets, \
					  const Eigen::Tensor<double,4> &ifc3Tensor, \
					  const Eigen::Tensor<double,3> &cellPositions,	\
					  const Eigen::Tensor<int,2> &displacedAtoms, const CrystalInfo &crysInfo){

  int iq1,iq2,iq3,ib,jb,s1,s2,s3,idim,iat;
  int i1x,i1y,i1z,i2x,i2y,i2z,i3x,i3y,i3z;
 
  Eigen::Tensor<complex <double>,3> ev1(3,crysInfo.numAtoms,crysInfo.numBranches); 
  Eigen::Tensor<complex <double>,3> ev2(3,crysInfo.numAtoms,crysInfo.numBranches);
  Eigen::Tensor<complex <double>,3> ev3(3,crysInfo.numAtoms,crysInfo.numBranches);
  
  // Edge lengths of BZ
  int nx = grid[0];
  int ny = grid[1];
  int nz = grid[2];
  int gridSize = nx*ny*nz;

  PhononTriplet interactingPhonons;
  PhInteraction3Ph phInt;

  int numAtoms = crysInfo.numAtoms;
  int numBranches = crysInfo.numBranches;
  Eigen::Tensor<double,4> Vm2(gridSize,numBranches,gridSize,numBranches);

  // Grab irred phonon mode info:
  iq1 = mode.iq; //index of wave vector in the full BZ
  s1 = mode.s; //branch
   
  //Demux 1st phonon wave vector
  //!!WARNING: For testing purposes using ShengBTE ordering!!!
  i1x = iq1%nx;
  i1y = (iq1/nx)%ny;
  i1z = iq1/nx/ny;
  
  interactingPhonons.s1 = s1;
  interactingPhonons.iq1 = iq1; 
  for(idim = 0; idim < 3; idim++){
    for(iat = 0; iat < numAtoms; iat++){
      for(ib = 0; ib < numBranches; ib++){
	ev1(idim,iat,ib) = ev(iq1,idim+3*iat,ib);
      }
    }
  }
  interactingPhonons.ev1 = ev1;

  int count = 0;
  
  // Loop over all 2nd phonon wave vectors
  for(i2z = 0; i2z < nz; i2z++){	
    for(i2y = 0; i2y < ny; i2y++){
      for(i2x = 0; i2x < nx; i2x++){
	//Muxed index of 2nd phonon wave vector
	//!!WARNING: For testing purposes using ShengBTE ordering!!!
	iq2 = (i2z*ny + i2y)*nx + i2x;

	interactingPhonons.iq2 = iq2;

	for(idim = 0; idim < 3; idim++){
	  for(iat = 0; iat < numAtoms; iat++){
	    for(ib = 0; ib < numBranches; ib++){
	      ev2(idim,iat,ib) = ev(iq2,idim+3*iat,ib);
	    }
	  }
	}
	interactingPhonons.ev2 = ev2;
	
	// Third phonon wave vector (Umklapped, if needed)
	i3x = (i1x - i2x + nx)%nx;
	i3y = (i1y - i2y + ny)%ny;
	i3z = (i1z - i2z + nz)%nz;
	//!!WARNING: For testing purposes using ShengBTE ordering!!!
	iq3 = (i3z*ny + i3y)*nx + i3x;
	
	interactingPhonons.iq3 = iq3;
	for(idim = 0; idim < 3; idim++){
	  for(iat = 0; iat < numAtoms; iat++){
	    for(ib = 0; ib < numBranches; ib++){
	      ev3(idim,iat,ib) = ev(iq3,idim+3*iat,ib);
	    }
	  }
	}
	interactingPhonons.ev3 = ev3;
	
	// Sum over 2nd phonon branches
	for(ib = 0; ib < numBranches; ib++){
	  interactingPhonons.s2 = ib;
	  // Sum over 3rd phonon branches
	  for(jb = 0; jb < numBranches; jb++){
	    interactingPhonons.s3 = jb;
	    
	    // Call calculateSingleV
	    Vm2(iq2,ib,iq3,jb) = phInt.calculateSingleV(interactingPhonons, qFBZ, numTriplets, ifc3Tensor, \
					 cellPositions, displacedAtoms, crysInfo, '-');

	    //cout << iq2+1 << " " << ib+1 << " " << iq3+1 << " " << jb+1 << " " << Vm2(iq2,ib,iq3,jb) << "\n";
	  }
	}
      }    
    }
    cout << ++count*2304.0/18432.0*100 << " % done.\n";
  }

  
  // Write to disk
  string fileName = "Vm2.iq"+to_string(iq1)+".s"+to_string(s1);
  ofstream outFile;
  //outFile.open(fileName, ios::out | ios::trunc | ios::binary);
  outFile.open(fileName,ios::trunc);
  for(iq2 = 0; iq2 < gridSize; iq2++){
    for(iq3 = 0; iq3 < gridSize; iq3++){
      for(ib = 0; ib < numBranches; ib++){
	for(jb = 0; jb < numBranches; jb++){
	  outFile << Vm2(iq2,ib,iq3,jb) << "\n";
	}
      }
    }
  }
  outFile.close();
}

//Transition probabilities for a given irreducible phonon mode
void PhInteraction3Ph::calculateAllW(const double T,const int grid[3], const PhononMode &mode,\
				     const Eigen::MatrixXi &indexMesh,const CrystalInfo &crysInfo,\
				     const Eigen::MatrixXd omega,const TetraData tetra){

  int iq1,iq2,iq3,iq3Plus,iq3Minus,s1,s2,s3,iDim,plusProcessCount,minusProcessCount;
  int numBranches = crysInfo.numBranches;
  int nq = grid[0]*grid[1]*grid[2];

  double omega1,omega2,omega3Plus,omega3Minus,n01,n02,n03Plus,n03Minus,Vplus2,Vminus2,\
    Wplus,Wminus,tetWeightPlus,tetWeightMinus;
  const double a = M_PI*hbar/4.0*5.60626442*1.0e30;

  Eigen::Vector3i q1,q2,q3Plus,q3Minus;
  
  // Grab irred phonon mode info:
  iq1 = mode.iq; //index of wave vector in the full BZ
  s1 = mode.s; //branch
  omega1 = omega(iq1,s1); //irred mode energy
  q1 = indexMesh.row(iq1);

  //Read full set of V-(q2,s2;q3,s3) for given mode from file
  Eigen::Tensor<double,4> Vm2(nq,numBranches,nq,numBranches);
  string fileName = "Vm2.iq"+to_string(iq1)+".s"+to_string(s1);
  ifstream inFile;
  inFile.open(fileName);
  for(iq2 = 0; iq2 < nq; iq2++){
    for(iq3 = 0; iq3 < nq; iq3++){
      for(s2 = 0; s2 < numBranches; s2++){
	for(s3 = 0; s3 < numBranches; s3++){
	  inFile >> Vm2(iq2,s2,iq3,s3);
	}
      }
    }
  }
  inFile.close();
  
  //Open W+, W-, and process counter files
  string WplusFileName = "Wp2.iq"+to_string(iq1)+".s"+to_string(s1);
  string WminusFileName = "Wm2.iq"+to_string(iq1)+".s"+to_string(s1);
  ofstream WplusFile, WminusFile;
  WplusFile.open(WplusFileName, ios::trunc);
  WminusFile.open(WminusFileName, ios::trunc);

  plusProcessCount = 0;
  minusProcessCount = 0;
  if(omega1 > 0){ //skip zero energy phonon
    //Bose distribution for first phonon mode
    n01 = bose(omega1,T);

    //Sum over second phonon wave vector in full BZ
    for(iq2 = 0; iq2 < nq; iq2++){
      q2 = indexMesh.row(iq2);

      //Get final phonon wave vector location
      //modulo reciprocal lattice vector
      for(iDim = 0; iDim < 3; iDim++){
	//plus process
	q3Plus(iDim) = (q1(iDim)+q2(iDim))%grid[iDim];
	//minus process
	q3Minus(iDim) = (q1(iDim)-q2(iDim)+grid[iDim])%grid[iDim];
      }
      //!!WARNING: For testing purposes using ShengBTE ordering!!!
      iq3Plus = (q3Plus(2)*grid[1] + q3Plus(1))*grid[0] + q3Plus(0);
      iq3Minus = (q3Minus(2)*grid[1] + q3Minus(1))*grid[0] + q3Minus(0);
      
      //Sum over second phonon branches
      for(s2 = 0; s2 < numBranches; s2++){
	omega2 = omega(iq2,s2); //second phonon ang. freq.

	if(omega2 > 0){ //skip zero energy phonon
	
	  //Bose distribution for second phonon mode
	  n02 = bose(omega2,T);
	  
	  //Sum over final phonon branches
	  for(s3 = 0; s3 < numBranches; s3++){
	    omega3Plus = omega(iq3Plus,s3); //third phonon(+) ang. freq.
	    omega3Minus = omega(iq3Minus,s3); //third phonon(-) ang. freq.
	  
	    //Bose distribution for final phonon mode
	    n03Plus = bose(omega3Plus,T);
	    n03Minus = bose(omega3Minus,T);

	    //Calculate tetrahedron weight for plus and minus processes
	    tetWeightPlus = fillTetsWeights(omega3Plus-omega1,s2,iq2,tetra);
	    tetWeightMinus = fillTetsWeights(omega1-omega3Minus,s2,iq2,tetra);
	    
	    //Plus processes
	    if(tetWeightPlus > 0 && omega3Plus > 0){
	      //Increase processes counter
	      plusProcessCount++;
				
	      //Time reverse second phonon to get Vplus2
	      Vplus2 = Vm2(timeReverse(iq2,grid),s2,iq3Plus,s3);
	      
	      //Calculatate transition probability W+
	      Wplus = a*(n02-n03Plus)*Vplus2*tetWeightPlus \
		/(omega1*omega2*omega3Plus); //THz
	      
	      //Write plus process info to file
	      WplusFile << iq2 << " " << s2 << " " << iq3Plus << " " << s3 << " "\
			<< Wplus << "\n";
	    }
	    
	    //Minus processes
	    if(tetWeightMinus > 0 && omega3Minus > 0){
	      //Increase processes counter
	      minusProcessCount++;
		
	      Vminus2 = Vm2(iq2,s2,iq3Minus,s3);

	      //Calculatate transition probability W-
	      Wminus = a*(n02+n03Minus+1.0)*Vminus2*tetWeightMinus	\
		/(omega1*omega2*omega3Minus); //THz

	      //Write minus process info to disk
	      WminusFile << iq2 << " " << s2 << " " << iq3Minus << " " << s3 << " "\
			<< Wminus << "\n";
	    }
	  }//s3
	}//zero of second phonon
      }//s2
    }//iq2
  }//zero of first phonon
  WplusFile.close();
  WminusFile.close();

  //Write total number of plus and minus processes to disk
  string counterFileName = "WCounter.iq"+to_string(iq1)+".s"+to_string(s1);
  ofstream counterFile;
  counterFile.open(counterFileName, ios::trunc);
  counterFile << plusProcessCount << "\n";
  counterFile << minusProcessCount << "\n";
  counterFile.close();
}


