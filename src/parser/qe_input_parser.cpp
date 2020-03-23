#include <math.h> // round()
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>  // to declare istringstream
#include <algorithm> // to use .remove_if
#include <stdlib.h> // abs()
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>

#include "constants.h"
#include "exceptions.h"
#include "qe_input_parser.h"
#include "phononH0.h"

double calcVolume(const Eigen::Matrix3d& directUnitCell, const double alat)
{
	Eigen::Vector3d a1 = directUnitCell.col(0);
	Eigen::Vector3d a2 = directUnitCell.col(1);
	Eigen::Vector3d a3 = directUnitCell.col(2);
	double volume;
	volume = abs( a1.dot(( a2.cross(a3) )) );
	volume+= abs( a2.dot(( a3.cross(a1) )) );
	volume+= abs( a3.dot(( a1.cross(a2) )) );
	volume *= alat * alat * alat / 3.;
	return volume;
}

Eigen::MatrixXd calcReciprocalCell(const Eigen::Matrix3d& directUnitCell)
{
	Eigen::Matrix3d reciprocalCell = directUnitCell.inverse().transpose();
	return reciprocalCell;
}

void latgen(const int ibrav, Eigen::VectorXd& celldm, Eigen::Matrix3d& unitCell)
{
	//  !     sets up the crystallographic vectors a1, a2, and a3.
	//  !
	//  !     ibrav is the structure index:
	//  !       1  cubic P (sc)                8  orthorhombic P
	//  !       2  cubic F (fcc)               9  1-face (C) centered orthorhombic
	//  !       3  cubic I (bcc)              10  all face centered orthorhombic
	//  !       4  hexagonal and trigonal P   11  body centered orthorhombic
	//  !       5  trigonal R, 3-fold axis c  12  monoclinic P (unique axis: c)
	//  !       6  tetragonal P (st)          13  one face (base) centered monoclinic
	//  !       7  tetragonal I (bct)         14  triclinic P
	//  !     Also accepted:
	//  !       0  "free" structure          -12  monoclinic P (unique axis: b)
	//  !      -3  cubic bcc with a more symmetric choice of axis
	//  !      -5  trigonal R, threefold axis along (111)
	//  !      -9  alternate description for base centered orthorhombic
	//  !     -13  one face (base) centered monoclinic (unique axis: b)
	//  !      91  1-face (A) centered orthorombic
	//  !
	//  !     celldm are parameters which fix the shape of the unit cell
	//  !     volumeUnitCell is the unit-cell volume
	//  !
	//  !     NOTA BENE: all axis sets are right-handed
	//  !     Boxes for US PPs do not work properly with left-handed axis

	const double sr2 = 1.414213562373, sr3 = 1.732050807569;

	//  user-supplied lattice vectors

	Eigen::Vector3d a1, a2, a3;

	a1 = unitCell.col(0);
	a2 = unitCell.col(1);
	a3 = unitCell.col(2);

	if ( ibrav == 0 ) {
		if ( sqrt( a1.transpose()*a1 ) == 0 || sqrt( a2.transpose()*a2 ) == 0
				|| sqrt( a3.transpose()*a3 ) == 0 ) {
			Error e("wrong at for ibrav=0", 1);
		}
		if ( celldm(0) != 0. ) {
			// ... input at are in units of alat => convert them to a.u.
			unitCell *= celldm(0);
		} else {
			// ... input at are in atomic units: define celldm(1) from a1
			celldm(0) = sqrt( a1.transpose() * a1 );
		}
	} else {
		a1.setZero();
		a2.setZero();
		a3.setZero();
	}

	if ( celldm(0) <= 0. ) {
		Error e("wrong celldm(1)", 1 );
	}

	//  index of bravais lattice supplied

	if ( ibrav == 1 ) { // simple cubic lattice
		a1(0) = celldm(0);
		a2(1) = celldm(0);
		a3(2) = celldm(0);
	} else if (ibrav == 2) { //     fcc lattice
		double term = celldm(0) / 2.;
		a1(0) =-term;
		a1(2) = term;
		a2(1) = term;
		a2(2) = term;
		a3(0) =-term;
		a3(1) = term;
	} else if (abs(ibrav) == 3) { // bcc lattice
		double term = celldm(0) / 2.;
		for ( int ir=0; ir<3; ir++ ) {
			a1(ir) = term;
			a2(ir) = term;
			a3(ir) = term;
		} if ( ibrav < 0 ) {
			a1(0) = -a1(0);
			a2(1) = -a2(1);
			a3(2) = -a3(2);
		} else {
			a2(0) = -a2(0);
			a3(0) = -a3(0);
			a3(1) = -a3(1);
		}
	} else if ( ibrav == 4 ) {// hexagonal lattice
		if ( celldm(2) <= 0. ) {
			Error e("wrong celldm(2)", ibrav);
		}
		double cbya  = celldm(2);
		a1(1) = celldm(0);
		a2(1) =-celldm(0) / 2.;
		a2(2) = celldm(0) * sr3 / 2.;
		a3(3) = celldm(0) * cbya;

	} else if (abs(ibrav) == 5) { // trigonal lattice
		if ( celldm(3) <= -0.5 || celldm(3) >= 1. ) {
			Error e("wrong celldm(4)", abs(ibrav));
		}

		double term1 = sqrt(1. + 2. * celldm(3) );
		double term2 = sqrt(1. - celldm(3) );

		if ( ibrav == 5 ) { // threefold axis along c (001)
			a2(1) = sr2 * celldm(0) * term2 / sr3;
			a2(2) = celldm(0) * term1 / sr3;
			a1(0) = celldm(0) * term2 / sr2;
			a1(1) =-a1(0) / sr3;
			a1(2) = a2(2);
			a3(0) =-a1(0);
			a3(1) = a1(1);
			a3(2) = a2(2);
		} else if ( ibrav == -5 ) { // threefold axis along (111)
			// Notice that in the cubic limit (alpha=90, celldm(4)=0, term1=term2=1)
			// does not yield the x,y,z axis, but an equivalent rotated triplet:
			//   a/3 (-1,2,2), a/3 (2,-1,2), a/3 (2,2,-1)
			// If you prefer the x,y,z axis as cubic limit, you should modify the
			// definitions of a1(1) and a1(2) as follows:'
			// a1(1) = celldm(1)*(term1+2.0_dp*term2)/3.0_dp
			// a1(2) = celldm(1)*(term1-term2)/3.0_dp
			// (info by G. Pizzi and A. Cepellotti)
			a1(0) = celldm(0) * ( term1 - 2. * term2 ) / 3.;
			a1(1) = celldm(0) * ( term1 + term2 ) / 3.;
			a1(2) = a1(1);
			a2(0) = a1(2);
			a2(1) = a1(0);
			a2(2) = a1(1);
			a3(0) = a1(1);
			a3(1) = a1(2);
			a3(2) = a1(0);
		}
	} else if (ibrav == 6) { // tetragonal lattice
		if ( celldm(2) <= 0. ) {
			Error e("wrong celldm(3)", 6);
		}
		double cbya = celldm(2);
		a1(0) = celldm(0);
		a2(1) = celldm(0);
		a3(2) = celldm(0) * cbya;

	} else if (ibrav == 7) { // body centered tetragonal lattice
		if ( celldm(2) <= 0. ) {
			Error e("wrong celldm(3)", 7);
		}
		double cbya = celldm(2);
		a2(0) = celldm(0) / 2.;
		a2(1) = a2(0);
		a2(2) = cbya * celldm(0) / 2.;
		a1(0) = a2(0);
		a1(1) = - a2(0);
		a1(2) = a2(2);
		a3(0) = - a2(0);
		a3(1) = - a2(0);
		a3(2) = a2(2);
	} else if ( ibrav == 8 ) { // Simple orthorhombic lattice
		if ( celldm(1) <= 0. ) {
			Error e("wrong celldm(2)", ibrav);
		}
		if ( celldm(2) <= 0. )
		{
			Error e("wrong celldm(3)", ibrav);
		}
		a1(0) = celldm(0);
		a2(1) = celldm(0) * celldm(1);
		a3(2) = celldm(0) * celldm(2);
	} else if ( abs(ibrav) == 9) { // One face (base) centered orthorhombic lattice  (C type)
		if ( celldm(1) <= 0. ) {
			Error e("wrong celldm(2)", abs(ibrav));
		}
		if ( celldm(2) <= 0. ) {
			Error e("wrong celldm(3)", abs(ibrav));
		}
		if ( ibrav == 9 ) {// old PWscf description
			a1(0) = 0.5 * celldm(0);
			a1(1) = a1(0) * celldm(1);
			a2(0) = - a1(0);
			a2(1) = a1(1);
		} else {// alternate description
			a1(0) =  0.5 * celldm(0);
			a1(1) = -a1(0) * celldm(1);
			a2(0) =  a1(0);
			a2(1) = -a1(1);
		}
		a3(2) = celldm(0) * celldm(2);
	} else if ( ibrav == 91 ) { // One face(base)centered orthorhombic lattice (A type)
		if ( celldm(1) <= 0. ) {
			Error e("wrong celldm(2)", ibrav);
		}
		if ( celldm(2) <= 0. ) {
			Error e("wrong celldm(3)", ibrav);
		}
		a1(0) = celldm(0);
		a2(1) = celldm(0) * celldm(1) * 0.5;
		a2(2) = - celldm(0) * celldm(2) * 0.5;
		a3(1) = a2(1);
		a3(2) = - a2(2);
	} else if (ibrav == 10) {// All face centered orthorhombic lattice
		if ( celldm(1) <= 0. ) {
			Error e("wrong celldm(2)", ibrav);
		}
		if ( celldm(2) <= 0. ) {
			Error e("wrong celldm(3)", ibrav);
		}
		a2(0) = 0.5 * celldm(0);
		a2(1) = a2(0) * celldm(1);
		a1(0) = a2(0);
		a1(2) = a2(0) * celldm(2);
		a3(1) = a2(0) * celldm(1);
		a3(2) = a1(2);
	} else if (ibrav == 11) { // Body centered orthorhombic lattice
		if ( celldm(1) <= 0. ) {
			Error e("wrong celldm(2)", ibrav);
		}
		if ( celldm(2) <= 0. ) {
			Error e("wrong celldm(3)", ibrav);
		}
		a1(0) = 0.5 * celldm(0);
		a1(1) = a1(0) * celldm(1);
		a1(2) = a1(0) * celldm(2);
		a2(0) = - a1(0);
		a2(1) = a1(1);
		a2(2) = a1(2);
		a3(0) = - a1(0);
		a3(1) = - a1(1);
		a3(2) = a1(2);
	} else if (ibrav == 12) { // Simple monoclinic lattice, unique (i.e. orthogonal to a) axis: c
		if ( celldm(1) <= 0. ) {
			Error e("wrong celldm(2)", ibrav);
		}
		if ( celldm(2) <= 0. ) {
			Error e("wrong celldm(3)", ibrav);
		}
		if ( abs(celldm(3)) >= 1. ) {
			Error e("wrong celldm(4)", ibrav);
		}
		double sen = sqrt( 1. - celldm(3)*celldm(3) );
		a1(0) = celldm(0);
		a2(0) = celldm(0) * celldm(1) * celldm(3);
		a2(1) = celldm(0) * celldm(1) * sen;
		a3(2) = celldm(0) * celldm(2);
	} else if ( ibrav == - 12 ) { // Simple monoclinic lattice, unique axis: b (more common)
		if ( celldm(1) <= 0. ) {
			Error e("wrong celldm(2)",-ibrav);
		}
		if ( celldm(2) <= 0. ) {
			Error e("wrong celldm(3)",-ibrav);
		}
		if ( abs(celldm(4))>=1. ) {
			Error e("wrong celldm(5)",-ibrav);
		}
		double sen = sqrt( 1. - celldm(4)*celldm(4) );
		a1(0) = celldm(0);
		a2(1) = celldm(0) * celldm(1);
		a3(0) = celldm(0) * celldm(2) * celldm(4);
		a3(2) = celldm(0) * celldm(2) * sen;
	} else if ( ibrav == 13 ) { // One face centered monoclinic lattice unique axis c
		if ( celldm(1) <= 0. ) {
			Error e("wrong celldm(2)", ibrav);
		}
		if ( celldm(2) <= 0. ) {
			Error e("wrong celldm(3)", ibrav);
		}
		if ( abs(celldm(3)) >= 1. ) {
			Error e("wrong celldm(4)", ibrav);
		}
		double sen = sqrt( 1. - celldm(4)*celldm(4) );
		a1(0) = 0.5 * celldm(0);
		a1(2) =-a1(0) * celldm(2);
		a2(0) = celldm(0) * celldm(1) * celldm(2);
		a2(1) = celldm(0) * celldm(1) * sen;
		a3(0) = a1(0);
		a3(2) =-a1(2);
	} else if ( ibrav == -13 ) { // One face centered monoclinic lattice unique axis b
		if ( celldm(1) <= 0. ) {
			Error e("wrong celldm(2)", -ibrav);
		}
		if ( celldm(2) <= 0. ) {
			Error e("wrong celldm(3)", -ibrav);
		}
		if ( abs(celldm(4)) >= 1. ) {
			Error e("wrong celldm(5)", -ibrav);
		}
		double sen = sqrt( 1. - celldm(4)*celldm(4) );
		a1(0) = 0.5 * celldm(0);
		a1(1) =-a1(0) * celldm(1);
		a2(0) = a1(0);
		a2(1) =-a1(1);
		a3(0) = celldm(0) * celldm(2) * celldm(4);
		a3(2) = celldm(0) * celldm(2) * sen;
	} else if (ibrav == 14) { // Triclinic lattice
		if ( celldm(1) <= 0. ) {
			Error e("wrong celldm(2)", ibrav);
		}
		if ( celldm(2) <= 0. ) {
			Error e("wrong celldm(3)", ibrav);
		}
		if ( abs(celldm(3)) >= 1. ) {
			Error e("wrong celldm(4)", ibrav);
		}
		if ( abs(celldm(4)) >= 1. ) {
			Error e("wrong celldm(5)", ibrav);
		}
		if ( abs(celldm(5)) >= 1. ) {
			Error e("wrong celldm(6)", ibrav);
		}
		double singam = sqrt( 1. - celldm(5)*celldm(5) );
		double term = ( 1. + 2. * celldm(3)*celldm(4)*celldm(5)
				- celldm(3)*celldm(3) - celldm(4)*celldm(4) - celldm(5)*celldm(5));
		if ( term < 0. )
		{
			Error e("celldm does not make sense, check your data", ibrav);
		}
		term = sqrt( term / ( 1. - celldm(5)*celldm(5) ) );
		a1(0) = celldm(0);
		a2(0) = celldm(0) * celldm(1) * celldm(5);
		a2(1) = celldm(0) * celldm(1) * singam;
		a3(0) = celldm(0) * celldm(2) * celldm(4);
		a3(1) = celldm(0) * celldm(2) * (celldm(3)-celldm(4)*celldm(5))/singam;
		a3(2) = celldm(0) * celldm(2) * term;

	} else {
		Error e("nonexistent bravais lattice", ibrav);
	}

	if ( ibrav != 0 ) {
		unitCell.col(0) = a1;
		unitCell.col(1) = a2;
		unitCell.col(2) = a3;
	}
}

std::vector<std::string> split(const std::string& s, char delimiter) {
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

	//  First line contains ibrav, celldm and other variables

	std::getline(infile, line);
	lineSplit = split(line, ' ');

	int numElements = std::stoi(lineSplit[0]);
	int numAtoms = std::stoi(lineSplit[1]);
	int ibrav = std::stoi(lineSplit[2]);

	Eigen::VectorXd celldm(6);
	celldm(0) = std::stod(lineSplit[3]);
	celldm(1) = std::stod(lineSplit[4]);
	celldm(2) = std::stod(lineSplit[5]);
	celldm(3) = std::stod(lineSplit[6]);
	celldm(4) = std::stod(lineSplit[7]);
	celldm(5) = std::stod(lineSplit[8]);

	Eigen::Matrix3d directUnitCell(3,3);
	if ( ibrav == 0 ) {
		// In this case, unitCell is written in the file, in angstroms
		for ( int i=0; i<3; i++ ) {
			std::getline(infile, line);
			lineSplit = split(line, ' ');
			for ( int j=0; j<3; j++) {
				directUnitCell(i,j) = std::stod(lineSplit[j]); // / distanceRyToAng;
			}
		};
	};

	// generate the unit cell vectors (also for ibrav != 0)
	latgen(ibrav, celldm, directUnitCell);

	double alat = celldm(0);
	directUnitCell /= alat; // bring unit cell in units of the lattice parameter

	//  Next, we read the atomic species
	std::vector<std::string> speciesNames;
	Eigen::VectorXd speciesMasses(numElements);
	for ( int i=0; i<numElements; i++ ) {
		std::getline(infile, line);
		lineSplit = split(line, '\'');
		speciesNames.push_back(lineSplit[1]);
		speciesMasses(i) = std::stod(lineSplit[2]); // in rydbergs
	};

	//  we read the atomic positions
	Eigen::MatrixXd atomicPositions(numAtoms,3);
	Eigen::VectorXi atomicSpecies(numAtoms);
	for ( int i=0; i<numAtoms; i++ ) {
		std::getline(infile, line);
		lineSplit = split(line, ' ');
		atomicSpecies(i) = std::stoi(lineSplit[1]) - 1;
		atomicPositions(i,0) = std::stod(lineSplit[2]);
		atomicPositions(i,1) = std::stod(lineSplit[3]);
		atomicPositions(i,2) = std::stod(lineSplit[4]);
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
	Eigen::MatrixXd dielectricMatrix(3,3);
	dielectricMatrix.setZero();
	Eigen::Tensor<double,3> bornCharges(numAtoms, 3, 3);
	bornCharges.setZero();

	if ( hasDielectric ) {
		for ( int i=0; i<3; i++) {
			std::getline(infile, line);
			lineSplit = split(line, ' ');
			for ( int j=0; j<3; j++) {
				dielectricMatrix(i,j) = std::stod(lineSplit[j]);
			}
		}

		for ( int iAtom=0; iAtom < numAtoms; iAtom++ ) {
			std::getline(infile, line);
			for ( int i=0; i<3; i++ ) {
				std::getline(infile, line);
				lineSplit = split(line, ' ');
				for ( int j=0; j<3; j++ ) {
					bornCharges(iAtom,i,j) = std::stod(lineSplit[j]);
				}
			}
		}
	}

	//	Now we parse the coarse q grid
	std::getline(infile, line);
	lineSplit = split(line, ' ');
	Eigen::VectorXi qCoarseGrid(3);
	qCoarseGrid(0) = std::stoi(lineSplit[0]);
	qCoarseGrid(1) = std::stoi(lineSplit[1]);
	qCoarseGrid(2) = std::stoi(lineSplit[2]);

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
								std::istringstream iss(line);
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

	// Now we do postprocessing

	double volumeUnitCell = calcVolume(directUnitCell, alat);
	Eigen::Matrix3d reciprocalUnitCell = calcReciprocalCell(directUnitCell);

	if ( qCoarseGrid(0) <= 0 || qCoarseGrid(1) <= 0 || qCoarseGrid(2) <= 0 ) {
		Error e("qCoarseGrid smaller than zero", 1);
	}

	//	Now, let's try to diagonalize some points, and start debugging at q=0


	PhononH0 dynamicalMatrix(directUnitCell,
			reciprocalUnitCell,
			alat,
			volumeUnitCell,
			atomicSpecies,
			speciesMasses,
			atomicPositions,
			dielectricMatrix,
			bornCharges,
			qCoarseGrid,
			forceConstants);

	dynamicalMatrix.setAcousticSumRule("crystal");

	Eigen::VectorXd omega(numAtoms*3);
	Eigen::Tensor<std::complex<double>,3> z(3,numAtoms,numAtoms*3);
	Eigen::VectorXd q(3);
	q << 0., 0., 0.;

	dynamicalMatrix.diagonalize(q, omega, z);

	std::cout << "Finito\n";
	std::cout << omega.transpose() * ryToCmm1 << std::endl;

	return;
};



