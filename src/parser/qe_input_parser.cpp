#include "qe_input_parser.h"
#include "constants.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>  // to declare istringstream
#include <algorithm> // to use .remove_if

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>

struct FileFormatNotRecognized : public std::exception {
	const char * what () const throw ()
    {
    	return "Error reading the file input parameter";
    }
};

//void setupmat(q,dyn,nat,at,bg,tau,itau_blk,nsc,alat, &
//	     &         dyn_blk,nat_blk,at_blk,bg_blk,tau_blk,omega_blk, &
//	     &         loto_2d, &
//	     &         epsil,zeu,frc,nr1,nr2,nr3,has_zstar,rws,nrws,na_ifc,f_of_q,fd)
//{
//
//	// compute the dynamical matrix (the analytic part only)
//
//	q_gen(nsc,qbid,at_blk,bg_blk,at,bg);
//
//	std::vector<double> qp[3];
//	int nb_blk;
//
//	for ( int iq=0; iq<nsc; iq++ ) {
//
//	     for ( int k=0; k<3; k++ ) {
//	        qp[k]= q[k] + qbid[k][iq];
//	     }
//
//	     std::fill(dyn_blk, dyn_blk+nat3*nat3, {0.,0.}); // from <algorithm>
//
//	     CALL frc_blk (dyn_blk,qp,tau_blk,nat_blk,              &
//	          &              nr1,nr2,nr3,frc,at_blk,bg_blk,rws,nrws,f_of_q,fd)
//	      IF (has_zstar .and. .not.na_ifc) &
//	           CALL rgd_blk(nr1,nr2,nr3,nat_blk,dyn_blk,qp,tau_blk,   &
//	                         epsil,zeu,bg_blk,omega_blk,celldm(1), loto_2d,+1.d0)
//	           ! LOTO 2D added celldm(1)=alat to passed arguments
//	     !
//
//		 for ( int iat=0; iat<numAtoms; iat++ ) {
//			iatBulk = itau_blk[iat];
//			for ( int jat=0; jat<numAtoms; jat++ ) {
//				jatBulk = itau_blk[jat];
//
//	            arg = twoPi * ( qp(1) * ( (tau(1,na)-tau_blk(1,na_blk)) -   &
//	                                (tau(1,nb)-tau_blk(1,nb_blk)) ) + &
//	                      qp(2) * ( (tau(2,na)-tau_blk(2,na_blk)) -   &
//	                                (tau(2,nb)-tau_blk(2,nb_blk)) ) + &
//	                      qp(3) * ( (tau(3,na)-tau_blk(3,na_blk)) -   &
//	                                (tau(3,nb)-tau_blk(3,nb_blk)) ) )
//	           cfac[jat] = { COS(arg)/nsc , SIN(arg)/nsc };
//			}
//
//			for ( int i=0; i<3; i++ ) {
//				for ( int j=0; j<3; j++ ) {
//					for ( int jat=0; jat<numAtoms; jat++ ) {
//						jatBulk = itau_blk[jat];
//		                dyn[i][j][iat][jat] = dyn[i][j][iat][jat] + cfac[jta] *
//		                     dyn_blk(i][j][iatBulk][jatBulk];
//					}
//				}
//			}
//		 } // iat
//	} // iq
//
//	  return
//}



void dyndiag(int numAtoms, std::vector<double> speciesMasses,
		std::vector<int> atomicSpecies,
		const Eigen::Tensor<std::complex<double>,4>& dyn,
		Eigen::VectorXd& w2, Eigen::Tensor<std::complex<double>,3>& z)
{
	// diagonalise the dynamical matrix
	// On input:  amass = masses, in amu
	// On output: w2 = energies, z = displacements

	// fill the two-indices dynamical matrix

  	int nat3 = 3 * numAtoms;
  	int iType, jType;
  	std::complex<double> cx;
  	Eigen::MatrixXcd dyn2Tmp(nat3, nat3);
  	Eigen::MatrixXcd dyn2(nat3, nat3);

	for (int iat = 0; iat<numAtoms; iat++) {
		for (int jat = 0; jat<numAtoms; jat++) {
			for (int ipol = 0; ipol<3; ipol++) {
				for (int jpol = 0; jpol<3; jpol++) {
					cx = dyn(ipol,jpol,iat,jat);
					dyn2Tmp(iat*3 + ipol, jat*3 + jpol) = cx;
				}
			}
		}
	}

    // impose hermiticity

//	dyn2 = dyn2 + dyn2.adjoint().eval();
//	dyn2 = dyn2Tmp + dyn2Tmp.adjoint();
	for ( int i=0; i<nat3; i++ ) {
		for ( int j=0; j<i-1; j++ ) {
			cx = ( dyn2(i,j) + std::conj(dyn2(j,i)) ) / 2.;
			dyn2(i,j) = cx;
			dyn2(j,i) = std::conj(cx);
		}
	}

    //  divide by the square root of masses

	for ( int iat=0; iat<numAtoms; iat++ ) {
		iType = atomicSpecies[iat];
		for ( int jat = 0; jat < numAtoms; jat++ ) {
			jType = atomicSpecies[jat];
			for ( int ipol = 0; ipol < 3; ipol++ ) {
				for ( int jpol = 0; jpol < 3; jpol++ ) {
					 dyn2(iat*3 + ipol, jat*3 + jpol) /=
							sqrt(speciesMasses[iType]*speciesMasses[jType]) / massRyToAmu;
				}
			}
		}
	}

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(dyn2);

    w2 = eigensolver.eigenvalues();

    Eigen::MatrixXcd zTemp = eigensolver.eigenvectors();

    //  displacements are eigenvectors divided by sqrt(amass)

    for ( int iband=0; iband<nat3; iband++ ) {
	    for ( int iat=0; iat<numAtoms; iat++ ) {
	    	iType = atomicSpecies[iat];
	    	for ( int ipol=0; ipol<3; ipol++ ) {
	    		z(ipol,iat,iband) = zTemp(iat*3 + ipol, iband) / sqrt(speciesMasses[iType] / massRyToAmu);
	    	}
	    }
    }
};

std::vector<std::string> split(const std::string& s, char delimiter)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);

   if ( delimiter == ' ' ) {
	   for (std::string s; tokenStream >> s; ) {
		   tokens.push_back(s);
	   }
   } else {
	   while (std::getline(tokenStream, token, delimiter)) {
		   token.erase(std::remove_if(token.begin(), token.end(), ::isspace),
				   token.end());
		   tokens.push_back(token);
	   }
   }

   return tokens;
}

void QEParser::parsePhHarmonic(std::string fileName) {
//  Here we read the dynamical matrix of interatomic force constants
//	in real space.
//	Since the file is typically small, we don't worry about memory management

	std::string line;
    std::vector<std::string> lineSplit;

// open input file
    std::ifstream infile(fileName);

//    this would read all content
//	std::vector<std::string> lines;
//	while (std::getline(infile, line)) {
//		lines.push_back(line);
//	}

//  First line contains ibrav, celldm and other variables

    std::getline(infile, line);
    lineSplit = split(line, ' ');

    int numElements = std::stoi(lineSplit[0]);
    int numAtoms = std::stoi(lineSplit[1]);
    int ibrav = std::stoi(lineSplit[2]);

    std::vector<double> celldm = {0.,0.,0.,0.,0.,0.};
    celldm[0] = std::stod(lineSplit[3]);
    celldm[1] = std::stod(lineSplit[4]);
    celldm[2] = std::stod(lineSplit[5]);
    celldm[3] = std::stod(lineSplit[6]);
    celldm[4] = std::stod(lineSplit[7]);
    celldm[5] = std::stod(lineSplit[8]);

	std::vector<std::vector<double>> unitCell(3, std::vector<double> (3, 0.));
    if ( ibrav == 0 ) {
//    	In this case, unitCell is written in the file, in angstroms
    	for ( int i=0; i<3; i++ ) {
			std::getline(infile, line);
			lineSplit = split(line, ' ');
			for ( int j=0; j<3; j++) {
				unitCell[i][j] = std::stod(lineSplit[j]) / distanceRyToAng;
			}
    	};
    };


//  Next, we read the atomic species
    std::vector<std::string> speciesNames;
    std::vector<double> speciesMasses;
    for ( int i=0; i<numElements; i++ ) {
		std::getline(infile, line);
		lineSplit = split(line, '\'');
		speciesNames.push_back(lineSplit[1]);
		speciesMasses.push_back(std::stod(lineSplit[2]) / massRyToAmu );
    };


    //  we read the atomic positions
    std::vector<std::vector<double>> atomicPositions(numAtoms,
    		std::vector<double> (3,0.));
    std::vector<int> atomicSpecies(numAtoms, 0);
    for ( int i=0; i<numAtoms; i++ ) {
		std::getline(infile, line);
		lineSplit = split(line, ' ');
		atomicSpecies[i] = std::stoi(lineSplit[1]);
		atomicPositions[i][0] = std::stod(lineSplit[2]);
		atomicPositions[i][1] = std::stod(lineSplit[3]);
		atomicPositions[i][2] = std::stod(lineSplit[4]);
    }

//  Read if hasDielectric
	std::getline(infile, line);
	line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
	bool hasDielectric;
	if ( line == "T" ) {
		hasDielectric = true;
	} else {
		hasDielectric = false;
	}

//	if there are the dielectric info, we can read dielectric matrix
//	and the Born charges
	if ( hasDielectric ) {
		std::vector<std::vector<double>> dielectricMatrix(3, std::vector<double> (3, 0.));
	    for ( int i=0; i<3; i++) {
			std::getline(infile, line);
			lineSplit = split(line, ' ');
			for ( int j=0; j<3; j++) {
				dielectricMatrix[i][j] = std::stod(lineSplit[j]);
			}
	    }

		std::vector<std::vector<std::vector<double>>> bornCharges(numAtoms,
				std::vector<std::vector<double>> (3,
						std::vector<double>(3,0.)) );
	    for ( int iAtom=0; iAtom < numAtoms; iAtom++ ) {
	    	std::getline(infile, line);
	    	for ( int i=0; i<3; i++ ) {
		    	std::getline(infile, line);
				lineSplit = split(line, ' ');
				for ( int j=0; j<3; j++ ) {
					bornCharges[iAtom][i][j] = std::stod(lineSplit[j]);
				}
	    	}
	    }
	}

//	Now we parse the coarse q grid
	std::getline(infile, line);
	lineSplit = split(line, ' ');
	std::vector<int> qCoarseGrid(3,0);
	qCoarseGrid[0] = std::stoi(lineSplit[0]);
	qCoarseGrid[1] = std::stoi(lineSplit[1]);
	qCoarseGrid[2] = std::stoi(lineSplit[2]);

//	dimensions of forceConstants
//	(qCoarseGrid[0], qCoarseGrid[1], qCoarseGrid[2], 3, 3, nat, nat)



	std::cout << qCoarseGrid[0];
	std::cout << qCoarseGrid[1];
	std::cout << qCoarseGrid[2];
	std::cout << numAtoms << "\n";

//  This would work:
//	double * forceConstants = new double[qCoarseGrid[0]*qCoarseGrid[1]*qCoarseGrid[2]*3*3*numAtoms*numAtoms];
//	*(forceConstants+0) = 1.2;
//	std::cout << *(forceConstants+0) << '\n';
//	delete [] forceConstants;

	Eigen::Tensor<double, 7> forceConstants(qCoarseGrid[0], qCoarseGrid[1],
		qCoarseGrid[2], 3, 3, numAtoms, numAtoms);

	int m1Test;
	int m2Test;
	int m3Test;
	double x;

	for ( int ic=0; ic<3; ic++ ) {
		for ( int jc=0; jc<3; jc++ ) {
			for ( int iat=0; iat<numAtoms; iat++ ) {
				for ( int jat=0; jat<numAtoms; jat++ ) {
					// a line containing ic, jc, iat, jat
					std::getline(infile, line);

					for ( int r3=0; r3<qCoarseGrid[2]; r3++ ) {
						for ( int r2=0; r2<qCoarseGrid[1]; r2++ ) {
							for ( int r1=0; r1<qCoarseGrid[0]; r1++ ) {
								std::getline(infile, line);
								istringstream iss(line);
								iss >> m1Test >> m2Test >> m3Test >> x;
								forceConstants(r1, r2, r3, ic, jc, iat, jat) = x;
							}
						}
					}
				}
			}
		}
	}

    infile.close();

    Eigen::MatrixXcd m(qCoarseGrid[2],2);



//	Since I'm here, let's try to diagonalize phonons
    // dyn(3,3,nat,nat)
    // z(3*nat,3*nat), w2(3*nat,nq)
//    set a q-point in cartesian coordinates

//    dyn = {0.,0.};





//
//
//	  ALLOCATE ( frc(nr1,nr2,nr3,3,3,nat,nat) )
//	  frc(:,:,:,:,:,:,:) = 0.d0
//	  DO i=1,3
//	     DO j=1,3
//	        DO na=1,nat
//	           DO nb=1,nat
//	              IF (ionode) READ (1,*) ibid, jbid, nabid, nbbid
//	              CALL mp_bcast(ibid,ionode_id, world_comm)
//	              CALL mp_bcast(jbid,ionode_id, world_comm)
//	              CALL mp_bcast(nabid,ionode_id, world_comm)
//	              CALL mp_bcast(nbbid,ionode_id, world_comm)
//	              IF(i .NE.ibid  .OR. j .NE.jbid .OR.                   &
//	                 na.NE.nabid .OR. nb.NE.nbbid)                      &
//	                 CALL errore  ('readfc','error in reading',1)
//	              IF (ionode) READ (1,*) (((m1bid, m2bid, m3bid,        &
//	                          frc(m1,m2,m3,i,j,na,nb),                  &
//	                           m1=1,nr1),m2=1,nr2),m3=1,nr3)
//
//	              CALL mp_bcast(frc(:,:,:,i,j,na,nb),ionode_id, world_comm)
//	           END DO
//	        END DO
//	     END DO
//	  END DO
//





//    std::cout << numElements << " , " << numAtoms <<  "!!!!\n";

	return;
};


