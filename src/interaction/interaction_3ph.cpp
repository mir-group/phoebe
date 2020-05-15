#include "ph_interaction_3ph.h"

// default constructor
Interaction3Ph::Interaction3Ph(Crystal & crystal_,
		long & numTriplets_,
		Eigen::Tensor<double,4> & ifc3Tensor_,
		Eigen::Tensor<double,3> & cellPositions_,
		Eigen::Tensor<long,2> & displacedAtoms_) :
		crystal(crystal_),
		numTriplets(numTriplets_),
		ifc3Tensor(ifc3Tensor_),
		cellPositions(cellPositions_),
		displacedAtoms(displacedAtoms_) {
}

// copy constructor
Interaction3Ph(const Interaction3Ph & that) :
	crystal(that.crystal),
	numTriplets(that.numTriplets),
	ifc3Tensor(that.ifc3Tensor),
	cellPositions(that.cellPositions),
	displacedAtoms(that.displacedAtoms) {
}

// assignment operator
Interaction3Ph & operator=(const Interaction3Ph & that) {
	if ( this != &that ) {
		crystal = that.crystal_;
		numTriplets = that.numTriplets_;
		ifc3Tensor = that.ifc3Tensor_;
		cellPositions = that.cellPositions_;
		displacedAtoms = that.displacedAtoms_;
	}
	return *this;
}

std::tuple<Eigen::Tensor<std::complex<double>,3>,
		Eigen::Tensor<std::complex<double>,3>>
		Interaction3Ph::getCouplingSquared(
		State & state1, State & state2, State & state3) {
	return calcCouplingSquared(state1, state2, state3);
}

std::tuple<Eigen::Tensor<std::complex<double>,3>,
		Eigen::Tensor<std::complex<double>,3>>
		Interaction3Ph::calcCouplingSquared(State & state1, State & state2,
				State & state3Plus, State & state3Mins ) {

	Eigen::Vector3d cell2Pos, cell3Pos;

	auto ev1 = state1.getEigenvectors(); // size (3,numAtoms,numBands)
	auto ev2 = state2.getEigenvectors();
	auto ev3Plus = state3Plus.getEigenvectors();
	auto ev3Mins = state3Mins.getEigenvectors();

	//phonon branches:
	// we allow the number of bands to be different in each direction
	long nb1 = ev1.dimension(2); // <- numBands
	long nb2 = ev2.dimension(2);
	long nb3Plus = ev3Plus.dimension(2);
	long nb3Mins = ev3Mins.dimension(2);

	//Cartesian phonon wave vectors: q1,q2,q3
	auto q1 = state1.getPoint().getCoords("cartesian");
	auto q2 = state2.getPoint().getCoords("cartesian");
	auto q3Plus = state3Plus.getPoint().getCoords("cartesian");
	auto q3Mins = state3Mins.getPoint().getCoords("cartesian");

	Eigen::Tensor<std::complex<double>,3> vPlus(nb1,nb2,nb3);
	Eigen::Tensor<std::complex<double>,3> vMins(nb1,nb2,nb3);
	vPlus.setZero();
	vMins.setZero();

	for ( long it=0; it<numTriplets; it++ ) { // sum over all triplets

		Eigen::Tensor<std::complex<double>,3> v0Plus(nb1,nb2,nb3);
		Eigen::Tensor<std::complex<double>,3> v0Mins(nb1,nb2,nb3);
		v0Plus.setZero();
		v0Mins.setZero();

		for ( int ib1=0; ib1<nb1; ib1++ ) {
			double omega1 = ...;
			if ( omega1 < energyCutoff ) continue;

			for ( int ib2=0; ib2<nb2; ib2++ ) {
				double omega2 = ...;
				if ( omega2 < energyCutoff ) continue;

				for ( int ic3 : {0,1,2} ) {
					for ( int ic2 : {0,1,2} ) {
						for ( int ic1 : {0,1,2} ) {
							for ( int ib3=0; ib3<nb3Plus; ib3++ ) {
								double omega3Plus = ;
								if ( omega3Plus < energyCutoff ) continue;
								double freqsPlus = omega1 * omega2 * omega3Plus;

								v0Plus(ib1,ib2,ib3) +=
										ifc3Tensor(it,ic1,ic2,ic3)
										* ev1(ic1,displacedAtoms(it,0),ib1)
										* ev2(ic2,displacedAtoms(it,1),ib2)
								* std::conj(ev3(ic3,displacedAtoms(it,2),ib3));
							}

							for ( int ib3=0; ib3<nb3Mins; ib3++ ) {
								double omega3Mins = ;
								if ( omega3Mins < energyCutoff ) continue;
								double freqsMins = omega1 * omega2 * omega3Mins;

								v0Mins(ib1,ib2,ib3) +=
										ifc3Tensor(it,ic1,ic2,ic3)
										* ev1(ic1,displacedAtoms(it,0),ib11)
								* std::conj(ev2(ic2,displacedAtoms(it,1),ib2))
								* std::conj(ev3(ic3,displacedAtoms(it,2),ib3));
							}
						}
					}
				}
			}
		}


		for ( int ic : {0,1,2} ) {
			cell2Pos(ic) = cellPositions(it,0,ic);
			cell3Pos(ic) = cellPositions(it,1,ic);
		}

		// As a convention, the first primitive cell in the triplet is
		// restricted to the origin, so the phase for that cell is unity.

		double arg = complexI *( q2.dot(cell2Pos) - q3.dot(cell3Pos) );
		std::complex phasePlus = exp( complexI * arg );
		arg = complexI * ( q2.dot(cell2Pos) + q3.dot(cell3Pos) );
		std::complex phaseMins = exp(-complexI * arg );

		for ( int ib1=0; ib1<nb1; ib1++ ) {
			for ( int ib2=0; ib2<nb2; ib2++ ) {
				for ( int ib3=0; ib3<nb3Plus; ib3++ ) {
					// case +
					vPlus(ib1,ib2,ib3) += v0Plus(ib1,ib2,ib3) * phasePlus;
				}
				for ( int ib3=0; ib3<nb3Mins; ib3++ ) {
					// case -
					vMins(ib1,ib2,ib3) += v0Mins(ib1,ib2,ib3) * phaseMins;
				}
			}
		}
	}

	Eigen::Tensor<double,3> couplingPlus(nb1,nb2,nb3);
	Eigen::Tensor<double,3> couplingMins(nb1,nb2,nb3);
	for ( int ib1=0; ib1<nb1; ib1++ ) {
		for ( int ib2=0; ib2<nb2; ib2++ ) {
			for ( int ib3=0; ib3<nb3Plus; ib3++ ) {
				// case +
				couplingPlus(ib1,ib2,ib3) += vPlus(ib1,ib2,ib3) *
						std::conj(vPlus(ib1,ib2,ib3));
			}
			for ( int ib3=0; ib3<nb3Mins; ib3++ ) {
				// case -
				couplingMins(ib1,ib2,ib3) += vMins(ib1,ib2,ib3) *
						std::conj(vMins(ib1,ib2,ib3));
			}
		}
	}
	return {couplingPlus, couplingMins};
}





//// Function to calculate the full set of V_minus processes for a given IBZ mode
//void PhInteraction3Ph::calculateAllVminus(const int grid[3], const PhononMode &mode,\
//		const Eigen::MatrixXd &qFBZ,	\
//		const Eigen::Tensor<complex<double>,3> &ev,const int numTriplets, \
//		const Eigen::Tensor<double,4> &ifc3Tensor, \
//		const Eigen::Tensor<double,3> &cellPositions,	\
//		const Eigen::Tensor<int,2> &displacedAtoms, const CrystalInfo &crysInfo){
//
//	int iq1,iq2,iq3,ib,jb,s1,s2,s3,idim,iat;
//	int i1x,i1y,i1z,i2x,i2y,i2z,i3x,i3y,i3z;
//
//	Eigen::Tensor<complex <double>,3> ev1(3,crysInfo.numAtoms,crysInfo.numBranches);
//	Eigen::Tensor<complex <double>,3> ev2(3,crysInfo.numAtoms,crysInfo.numBranches);
//	Eigen::Tensor<complex <double>,3> ev3(3,crysInfo.numAtoms,crysInfo.numBranches);
//
//	// Edge lengths of BZ
//	int nx = grid[0];
//	int ny = grid[1];
//	int nz = grid[2];
//	int gridSize = nx*ny*nz;
//
//	PhononTriplet interactingPhonons;
//	PhInteraction3Ph phInt;
//
//	int numAtoms = crysInfo.numAtoms;
//	int numBranches = crysInfo.numBranches;
//	Eigen::Tensor<double,4> Vm2(gridSize,numBranches,gridSize,numBranches);
//
//	// Grab irred phonon mode info:
//	iq1 = mode.iq; //index of wave vector in the full BZ
//	s1 = mode.s; //branch
//
//	//Demux 1st phonon wave vector
//	//!!WARNING: For testing purposes using ShengBTE ordering!!!
//	i1x = iq1%nx;
//	i1y = (iq1/nx)%ny;
//	i1z = iq1/nx/ny;
//
//	interactingPhonons.s1 = s1;
//	interactingPhonons.iq1 = iq1;
//	for(idim = 0; idim < 3; idim++){
//		for(iat = 0; iat < numAtoms; iat++){
//			for(ib = 0; ib < numBranches; ib++){
//				ev1(idim,iat,ib) = ev(iq1,idim+3*iat,ib);
//			}
//		}
//	}
//	interactingPhonons.ev1 = ev1;
//
//	int count = 0;
//
//	// Loop over all 2nd phonon wave vectors
//	for(i2z = 0; i2z < nz; i2z++){
//		for(i2y = 0; i2y < ny; i2y++){
//			for(i2x = 0; i2x < nx; i2x++){
//				//Muxed index of 2nd phonon wave vector
//				//!!WARNING: For testing purposes using ShengBTE ordering!!!
//				iq2 = (i2z*ny + i2y)*nx + i2x;
//
//				interactingPhonons.iq2 = iq2;
//
//				for(idim = 0; idim < 3; idim++){
//					for(iat = 0; iat < numAtoms; iat++){
//						for(ib = 0; ib < numBranches; ib++){
//							ev2(idim,iat,ib) = ev(iq2,idim+3*iat,ib);
//						}
//					}
//				}
//				interactingPhonons.ev2 = ev2;
//
//				// Third phonon wave vector (Umklapped, if needed)
//				i3x = (i1x - i2x + nx)%nx;
//				i3y = (i1y - i2y + ny)%ny;
//				i3z = (i1z - i2z + nz)%nz;
//				//!!WARNING: For testing purposes using ShengBTE ordering!!!
//				iq3 = (i3z*ny + i3y)*nx + i3x;
//
//				interactingPhonons.iq3 = iq3;
//				for(idim = 0; idim < 3; idim++){
//					for(iat = 0; iat < numAtoms; iat++){
//						for(ib = 0; ib < numBranches; ib++){
//							ev3(idim,iat,ib) = ev(iq3,idim+3*iat,ib);
//						}
//					}
//				}
//				interactingPhonons.ev3 = ev3;
//
//				// Sum over 2nd phonon branches
//				for(ib = 0; ib < numBranches; ib++){
//					interactingPhonons.s2 = ib;
//					// Sum over 3rd phonon branches
//					for(jb = 0; jb < numBranches; jb++){
//						interactingPhonons.s3 = jb;
//
//						// Call calculateSingleV
//						Vm2(iq2,ib,iq3,jb) = phInt.calculateSingleV(interactingPhonons, qFBZ, numTriplets, ifc3Tensor, \
//								cellPositions, displacedAtoms, crysInfo, '-');
//
//						//cout << iq2+1 << " " << ib+1 << " " << iq3+1 << " " << jb+1 << " " << Vm2(iq2,ib,iq3,jb) << "\n";
//					}
//				}
//			}
//		}
//		cout << ++count*2304.0/18432.0*100 << " % done.\n";
//	}
//
//
//	// Write to disk
//	string fileName = "Vm2.iq"+to_string(iq1)+".s"+to_string(s1);
//	ofstream outFile;
//	//outFile.open(fileName, ios::out | ios::trunc | ios::binary);
//	outFile.open(fileName,ios::trunc);
//	for(iq2 = 0; iq2 < gridSize; iq2++){
//		for(iq3 = 0; iq3 < gridSize; iq3++){
//			for(ib = 0; ib < numBranches; ib++){
//				for(jb = 0; jb < numBranches; jb++){
//					outFile << Vm2(iq2,ib,iq3,jb) << "\n";
//				}
//			}
//		}
//	}
//	outFile.close();
//}
//
////Transition probabilities for a given irreducible phonon mode
//void PhInteraction3Ph::calculateAllW(const double T,const int grid[3], const PhononMode &mode,\
//		const Eigen::MatrixXi &indexMesh,const CrystalInfo &crysInfo,\
//		const Eigen::MatrixXd omega,const TetraData tetra){
//
//	int iq1,iq2,iq3,iq3Plus,iq3Minus,s1,s2,s3,iDim,plusProcessCount,minusProcessCount;
//	int numBranches = crysInfo.numBranches;
//	int nq = grid[0]*grid[1]*grid[2];
//
//	double omega1,omega2,omega3Plus,omega3Minus,n01,n02,n03Plus,n03Minus,Vplus2,Vminus2,\
//	Wplus,Wminus,tetWeightPlus,tetWeightMinus;
//	const double a = M_PI*hbar/4.0*5.60626442*1.0e30;
//
//	Eigen::Vector3i q1,q2,q3Plus,q3Minus;
//
//	// Grab irred phonon mode info:
//	iq1 = mode.iq; //index of wave vector in the full BZ
//	s1 = mode.s; //branch
//	omega1 = omega(iq1,s1); //irred mode energy
//	q1 = indexMesh.row(iq1);
//
//	//Read full set of V-(q2,s2;q3,s3) for given mode from file
//	Eigen::Tensor<double,4> Vm2(nq,numBranches,nq,numBranches);
//	string fileName = "Vm2.iq"+to_string(iq1)+".s"+to_string(s1);
//	ifstream inFile;
//	inFile.open(fileName);
//	for(iq2 = 0; iq2 < nq; iq2++){
//		for(iq3 = 0; iq3 < nq; iq3++){
//			for(s2 = 0; s2 < numBranches; s2++){
//				for(s3 = 0; s3 < numBranches; s3++){
//					inFile >> Vm2(iq2,s2,iq3,s3);
//				}
//			}
//		}
//	}
//	inFile.close();
//
//	//Open W+, W-, and process counter files
//	string WplusFileName = "Wp2.iq"+to_string(iq1)+".s"+to_string(s1);
//	string WminusFileName = "Wm2.iq"+to_string(iq1)+".s"+to_string(s1);
//	ofstream WplusFile, WminusFile;
//	WplusFile.open(WplusFileName, ios::trunc);
//	WminusFile.open(WminusFileName, ios::trunc);
//
//	plusProcessCount = 0;
//	minusProcessCount = 0;
//	if(omega1 > 0){ //skip zero energy phonon
//		//Bose distribution for first phonon mode
//		n01 = bose(omega1,T);
//
//		//Sum over second phonon wave vector in full BZ
//		for(iq2 = 0; iq2 < nq; iq2++){
//			q2 = indexMesh.row(iq2);
//
//			//Get final phonon wave vector location
//			//modulo reciprocal lattice vector
//			for(iDim = 0; iDim < 3; iDim++){
//				//plus process
//				q3Plus(iDim) = (q1(iDim)+q2(iDim))%grid[iDim];
//				//minus process
//				q3Minus(iDim) = (q1(iDim)-q2(iDim)+grid[iDim])%grid[iDim];
//			}
//			//!!WARNING: For testing purposes using ShengBTE ordering!!!
//			iq3Plus = (q3Plus(2)*grid[1] + q3Plus(1))*grid[0] + q3Plus(0);
//			iq3Minus = (q3Minus(2)*grid[1] + q3Minus(1))*grid[0] + q3Minus(0);
//
//			//Sum over second phonon branches
//			for(s2 = 0; s2 < numBranches; s2++){
//				omega2 = omega(iq2,s2); //second phonon ang. freq.
//
//				if(omega2 > 0){ //skip zero energy phonon
//
//					//Bose distribution for second phonon mode
//					n02 = bose(omega2,T);
//
//					//Sum over final phonon branches
//					for(s3 = 0; s3 < numBranches; s3++){
//						omega3Plus = omega(iq3Plus,s3); //third phonon(+) ang. freq.
//						omega3Minus = omega(iq3Minus,s3); //third phonon(-) ang. freq.
//
//						//Bose distribution for final phonon mode
//						n03Plus = bose(omega3Plus,T);
//						n03Minus = bose(omega3Minus,T);
//
//						//Calculate tetrahedron weight for plus and minus processes
//						tetWeightPlus = fillTetsWeights(omega3Plus-omega1,s2,iq2,tetra);
//						tetWeightMinus = fillTetsWeights(omega1-omega3Minus,s2,iq2,tetra);
//
//						//Plus processes
//						if(tetWeightPlus > 0 && omega3Plus > 0){
//							//Increase processes counter
//							plusProcessCount++;
//
//							//Time reverse second phonon to get Vplus2
//							Vplus2 = Vm2(timeReverse(iq2,grid),s2,iq3Plus,s3);
//
//							//Calculatate transition probability W+
//							Wplus = a*(n02-n03Plus)*Vplus2*tetWeightPlus \
//									/(omega1*omega2*omega3Plus); //THz
//
//							//Write plus process info to file
//							WplusFile << iq2 << " " << s2 << " " << iq3Plus << " " << s3 << " "\
//									<< Wplus << "\n";
//						}
//
//						//Minus processes
//						if(tetWeightMinus > 0 && omega3Minus > 0){
//							//Increase processes counter
//							minusProcessCount++;
//
//							Vminus2 = Vm2(iq2,s2,iq3Minus,s3);
//
//							//Calculatate transition probability W-
//							Wminus = a*(n02+n03Minus+1.0)*Vminus2*tetWeightMinus	\
//									/(omega1*omega2*omega3Minus); //THz
//
//							//Write minus process info to disk
//							WminusFile << iq2 << " " << s2 << " " << iq3Minus << " " << s3 << " "\
//									<< Wminus << "\n";
//						}
//					}//s3
//				}//zero of second phonon
//			}//s2
//		}//iq2
//	}//zero of first phonon
//	WplusFile.close();
//	WminusFile.close();
//
//	//Write total number of plus and minus processes to disk
//	string counterFileName = "WCounter.iq"+to_string(iq1)+".s"+to_string(s1);
//	ofstream counterFile;
//	counterFile.open(counterFileName, ios::trunc);
//	counterFile << plusProcessCount << "\n";
//	counterFile << minusProcessCount << "\n";
//	counterFile.close();
//}
