#include "gtest/gtest.h"
#include "points.h"
#include "state.h"
#include "ifc3_parser.h"
#include "qe_input_parser.h"
#include "ph_scattering.h"
#include "bandstructure.h"
#include <fstream>

//TEST (Interaction3Ph, Coupling3Ph) {
//	// in this test, we compute the 3ph interaction |V|^2 for three arbitrary
//	// points, and compare the results with the coupling computed by ShengBTE
//
//	Context context;
//	context.setPhD2FileName("interaction3ph/QEspresso.fc");
//	context.setPhD3FileName("interaction3ph/ShengBTEForceConstants3rd");
//	Eigen::VectorXd temperatures(1);
//	temperatures << 300.;
//	context.setTemperatures(temperatures);
//	context.setSmearingMethod(2); // tetrahedron smearing
//
//	QEParser qeParser;
//	auto [crystal,phononH0] = qeParser.parsePhHarmonic(context);
//
//	IFC3Parser ifc3Parser;
//	auto coupling3Ph = ifc3Parser.parseFromShengBTE(context, crystal);
//
//	// set up the mesh of qpoints of shengbte
//	Eigen::Vector3i qMesh;
//	qMesh << 2, 2, 2;
//	FullPoints points(crystal, qMesh);
//	bool withVelocities = true;
//	bool withEigenvectors = true;
//	FullBandStructure<FullPoints> bandStructure = phononH0.populate(points,
//			withVelocities, withEigenvectors);
//
//	long iq1 = 2;
//	long iq2 = 7;
//
//	auto q1 = bandStructure.getPoint(iq1);
//	auto q2 = bandStructure.getPoint(iq2);
//
//	auto states1 = bandStructure.getState(iq1);
//
//	long iq2Inv = bandStructure.getPoints().getIndexInverted(iq2);
//
//	auto states2 = bandStructure.getState(iq2);
//	auto states2Plus = bandStructure.getState(iq2Inv);
//
//	auto q3Plus = q1 + q2;
//	auto q3Mins = q1 - q2;
//
//	auto states3Plus = bandStructure.getState(q3Plus);
//	auto states3Mins = bandStructure.getState(q3Mins);
//
//	auto [couplingPlus, couplingMins] = coupling3Ph.getCouplingSquared(
//			states1, states2Plus, states2,
//			states3Plus, states3Mins);
//	EXPECT_EQ ( 0., 0.);
//}


// this is a test with all Phoebe components
//int main() {
//	Context context;
//	context.setPhD2FileName("interaction3ph/QEspresso.fc");
//	context.setPhD3FileName("interaction3ph/ShengBTEForceConstants3rd");
//	Eigen::VectorXd temperatures(1);
//	temperatures << 300.;
//	context.setTemperatures(temperatures);
//
//	QEParser qeParser;
//	auto [crystal,phononH0] = qeParser.parsePhHarmonic(context);
//
//	IFC3Parser ifc3Parser;
//	auto coupling3Ph = ifc3Parser.parseFromShengBTE(context, crystal);
//
//	// set up the mesh of qpoints of shengbte
//	Eigen::Vector3i qMesh;
//	qMesh << 2, 2, 2;
//	FullPoints points(crystal, qMesh);
//	bool withVelocities = true;
//	bool withEigenvectors = true;
//	FullBandStructure<FullPoints> bandStructure = phononH0.populate(points,
//			withVelocities, withEigenvectors);
//
//	long iq1 = 1;
//	long iq2 = 7;
//
//	auto q1 = bandStructure.getPoint(iq1);
//	auto q2 = bandStructure.getPoint(iq2);
//
//	std::cout << q1.getCoords().transpose() << " : 1\n";
//	std::cout << q2.getCoords().transpose() << " : 2 \n";
//
//	auto states1 = bandStructure.getState(iq1);
//
//	long iq2Inv = bandStructure.getPoints().getIndexInverted(iq2);
//
//	auto states2 = bandStructure.getState(iq2);
//	auto states2Plus = bandStructure.getState(iq2Inv);
//
//	std::cout << bandStructure.getPoint(iq2Inv).getCoords().transpose() << " :2- \n";
//
////	auto q3Plus = q1 + q2;
////	auto q3Mins = q1 - q2;
//
//	auto q3Plus = bandStructure.getPoint(6);
//	auto q3Mins = bandStructure.getPoint(6);
//
//	std::cout << q3Plus.getCoords().transpose() << "\n";
//	std::cout << q3Mins.getCoords().transpose() << "\n";
//
//	auto states3Plus = bandStructure.getState(6);
//	auto states3Mins = bandStructure.getState(6);
//
//	auto [couplingPlus, couplingMins] = coupling3Ph.getCouplingSquared(
//			states1, states2Plus, states2,
//			states3Plus, states3Mins);
//
//	auto en1 = states1.getEnergies();
//	auto en2 = states2.getEnergies();
//	auto en3P = states3Plus.getEnergies();
//	auto en3M = states3Mins.getEnergies();
//
//	int numBands = en1.size();
//
//	std::cout << en1.transpose() * ryToCmm1 << "\n";
//
//	Eigen::Tensor<std::complex<double>,3> ev1;
//	states1.getEigenvectors(ev1);
//	for ( int k : {0,1,2,3,4,5} ) {
//	for ( int i : {0,1,2} ) {
//		for ( int j : {0,1} ) {
//			std::cout << i << " " << j << " " << ev1(i,j,k) << "\n";
//		}
//	}
//	std::cout << "\n";
//	}
//	std::cout << "\n";
//
//	for ( int i=0; i<numBands; i++ ) {
//		for ( int j=0; j<numBands; j++ ) {
//			for ( int k=0; k<numBands; k++ ) {
//				std::cout << i+1 << " " << j+1 << " " << k+1 << " "
//						<< couplingPlus(i,j,k) //* en1(i) * en2(j) * en3P(k)
//							* energyRyToEv * energyRyToEv
//							/ pow(distanceBohrToAng,6)
//							/ pow(massRyToAmu,3)
//						<< " "
//						<< couplingMins(i,j,k) //* en1(i) * en2(j) * en3M(k)
//							* energyRyToEv * energyRyToEv
//							/ pow(distanceBohrToAng,6)
//							/ pow(massRyToAmu,3)
//						<< "\n";
//			}
//		}
//		std::cout << "\n";
//	}
//
//	auto q0 = bandStructure.getPoint(0);
//	std::cout << q0.getCoords().transpose() << "\n";
//	auto states0 = bandStructure.getState(0);
//	std::cout << states0.getEnergies().transpose() * ryToCmm1 << "\n";
//	Eigen::MatrixXcd ev0;
//	states1.getEigenvectors(ev0);
//	for ( int i : {0,1,2,3,4,5} ) {
//		for ( int j : {0,1,2,3,4,5} ) {
//			std::cout << ev0(i,j) * sqrt(25598.822984527003) << "\n";
//		}
//	}
//	std::cout << "\n";
//
//	std::cout << crystal.getSpeciesMasses().transpose() << "\n";
//}


TEST (Interaction3Ph, Coupling3Ph) {
	// here we test the 3-ph matrix element (squared)
	// against results produced by ShengBTE.
	// in particular, we compute the couplings at q1=0, q2=0, q3=0

	// Note: Eigenvectors are complex: they can get a phase
	//       the phase doesn't change the coupling
	// Note: Eigenvectors can be degenerate! Hence, it's not straightforward to
	//       compare the coupling with ShengBTE or other codes, since the
	//       coupling at fixed mode indices might be different. However, the
	//       SUM of the coupling matrix elements within a degenerate subspace
	//       is the correct invariant quantity to be compared.
	// In this test, specifically, we use the three-fold degenerate modes at
	// gamma, and we compare the sum of the interaction within the optical
	// manifold

	Context context;
	context.setPhD2FileName("test/interaction3ph/QEspresso.fc");
	context.setPhD3FileName("test/interaction3ph/ShengBTEForceConstants3rd");

	QEParser qeParser;
	auto [crystal,phononH0] = qeParser.parsePhHarmonic(context);
	phononH0.setAcousticSumRule("simple");

	IFC3Parser ifc3Parser;
	auto coupling3Ph = ifc3Parser.parseFromShengBTE(context, crystal);

	//const double amu=1.66053904e-27; //Kg

	//------Set up cubic silicon-------//
	//Grid points along the 3 lattice vectors
//	long grid[3] = {20,20,20};
	//Total number of wave vectors
//	long nq = 20*20*20;
	//Number of atoms
	long numAtoms = crystal.getNumAtoms();
	//Number of bands
	long numBands = 3*numAtoms;
	//Number of types of atoms
//	long numTypes = crystal.getNumSpecies();
	//Types of atoms
//	Eigen::VectorXi types(numAtoms);
//	types << 0,0;
	//Masses of the types of atoms
//	Eigen::VectorXd masses(numTypes);
//	masses << 28.085509989600006; //amu
	//fill CrystalInfo
//	CrystalInfo crysInfo;
//	crysInfo.numAtoms = numAtoms;
//	crysInfo.numBands = numBands;
//	crysInfo.types = types;
//	crysInfo.masses = masses;
	//------------end setup-----------//

	// Read ifc3 from file
//	std::string ifc3FileName = "./interaction3ph/FORCE_CONSTANTS_3RD";
//	long numTriplets;
//	Eigen::Tensor<double,4> ifc3Tensor;
//	Eigen::Tensor<double,3> cellPositions;
//	Eigen::Tensor<long,2> displacedAtoms;

	// Read full BZ wave vectors (Cartesian) from file
//	std::string qFileName = "./interaction3ph/qpoints_full_cart";
//	Eigen::MatrixXd q(nq,3);
//	std::ifstream qfile(qFileName);
//	for ( long iq = 0; iq < nq; iq++ ) {
//		for ( long ip = 0; ip < 3; ip++ ) {
//			double x;
//			qfile >> x;
//			q(iq,ip) = x * distanceBohrToAng;
//		}
//	}
//	qfile.close();

	// Read full BZ eigenvectors from file
	// Note that this is read as outputted by ShengBTE
	// and has to be reshaped to conform to PHOEBE's format
	// of the eigenvector.
//	std::string evFileName = "./interaction3ph/evecs_full";
//	Eigen::Tensor<complex<double>,3> ev(nq,numBands,numBands);
//	ifstream evfile(evFileName);
//	double re,im;
//	complex<double> z;

	auto atomicMasses = crystal.getAtomicMasses();

//	for(long iq = 0; iq < nq; iq++){
//		for(long ib = 0; ib < numBands; ib++){
//			auto [idim,ia] = decompress2Indeces(ib,3,numAtoms);
//			for(long jb = 0; jb < numBands; jb++){
//				evfile >> re >> im;
//				z = std::complex<double>(re,im);
//				ev(iq,ib,jb) = z / sqrt(atomicMasses(ia));
//			}
//		}
//	}
//	evfile.close();

	// Create indexMesh
	// For testing purposes this is in ShengBTE ordering
//	Eigen::MatrixXi indexMesh(nq,3);
//	long count = 0;
//	for ( long k = 0; k < grid[2]; k++ ) {
//		for ( long j = 0; j < grid[1]; j++ ) {
//			for ( long i = 0; i < grid[0]; i++ ) {
//				indexMesh(count,0) = i;
//				indexMesh(count,1) = j;
//				indexMesh(count,2) = k;
//				count++;
//			}
//		}
//	}

	// Form a triplet to test vertex calculator
//	long iq1 = 0;
//	long iq2 = 0;
//	long iq3 = 0;
	Eigen::VectorXd energies(numBands);
	energies.setConstant(1.);

	Eigen::Vector3d q1, q2, q3;
//	q1 = q.row(iq1);
//	q2 = q.row(iq2);
//	q3 = q.row(iq3);
	q1.setZero();
	q2.setZero();
	q3.setZero();

	// Reshape the eigenvectors read from file
//	Eigen::Tensor<complex<double>,3> ev1(3,numAtoms,numBands);
//	Eigen::Tensor<complex<double>,3> ev2(3,numAtoms,numBands);
//	Eigen::Tensor<complex<double>,3> ev3(3,numAtoms,numBands);
//
//	for ( long idim = 0; idim < 3; idim++ ) {
//		for ( long iat = 0; iat < numAtoms; iat++ ) {
//			for ( long ib = 0; ib < numBands; ib++ ) {
//				ev1(idim,iat,ib) = ev(iq1,idim+3*iat,ib);
//				ev2(idim,iat,ib) = ev(iq2,idim+3*iat,ib);
//				ev3(idim,iat,ib) = ev(iq3,idim+3*iat,ib);
//			}
//		}
//	}
//	Eigen::MatrixXcd evm1(numBands,numBands),
//			evm2(numBands,numBands), evm3(numBands,numBands);
//	for ( int i=0; i<numBands; i++ ) {
//		for ( int j=0; j<numBands; j++ ) {
//			auto [idim,iat] = decompress2Indeces(i,3,numAtoms);
//			evm1(i,j) = ev1(idim,iat,j);
//			evm2(i,j) = ev2(idim,iat,j);
//			evm3(i,j) = ev3(idim,iat,j);
//		}
//	}

	auto [energies1,ev1] = phononH0.diagonalizeFromCoords(q1);
	auto [energies2,ev2] = phononH0.diagonalizeFromCoords(q2);
	auto [energies3,ev3] = phononH0.diagonalizeFromCoords(q3);

	Eigen::MatrixXcd evm1(numBands,numBands), evm2(numBands,numBands),
			evm3(numBands,numBands);

	// note: I want eigenvectors to be normalized with the masses
	evm1 = ev1 / sqrt(atomicMasses(0));
	evm2 = ev2 / sqrt(atomicMasses(0));
	evm3 = ev3 / sqrt(atomicMasses(0));

	DetachedState s1(q1, energies, numAtoms, numBands, evm1, nullptr);
	DetachedState s2(q2, energies, numAtoms, numBands, evm2, nullptr);
	DetachedState s3(q3, energies, numAtoms, numBands, evm3, nullptr);

	auto [couplingPlus,couplingMins] = coupling3Ph.getCouplingSquared(
														s1, s2, s3, s3);

	// we load reference data

	Eigen::Tensor<double,3> referenceCoupling(numBands,numBands,numBands);
	referenceCoupling.setZero();
	std::ifstream tfile("./test/interaction3ph/referenceCmins000");
	int i_, j_, k_;
	double x;
    if (tfile.is_open()) {
		for ( long i=0; i<numBands; i++ ) {
			for ( long j=0; j<numBands; j++ ) {
				for ( long k=0; k<numBands; k++ ) {
					tfile >> i_ >> j_ >> k_ >> x;
					referenceCoupling(i,j,k) = x;
				}
			}
		}
    }
	tfile.close();

	// now we compare the difference

	//	for ( int i=0; i<numBands; i++ ) {
	//		for ( int j=0; j<numBands; j++ ) {
	//			for ( int k=0; k<numBands; k++ ) {
	//				std::cout << i+1 << " " << j+1 << " " << k+1 << " "
	//						<< couplingPlus(i,j,k)
	//						<< " "
	//						<< couplingMins(i,j,k)
	//						<< " "
	//						<< referenceCoupling(i,j,k)
	//						<< "\n";
	//			}
	//		}
	//		std::cout << "\n";
	//	}

	double x1,x2;
	x1 = 0.;
	x2 = 0.;
	for ( int i=3; i<numBands; i++ ) {
		for ( int j=3; j<numBands; j++ ) {
			for ( int k=3; k<numBands; k++ ) {
				x2 += couplingPlus(i,j,k);
				x1 += referenceCoupling(i,j,k);
			}
		}
	}
	double relativeError = abs((x1-x2)/x1);
//	std::cout << x1 << " " << x2 << "\n";
//	std::cout << "Relative error on degenerate space: " << relativeError << "\n";

	ASSERT_NEAR(relativeError, 0., 1.0e-4);

//	if ( relativeError > 1.0e-4 ) {
//		Error e("Coupling 3ph is different than the reference values");
//	}

	// Note: to compare results with coupling from ShengBTE, multiply
	// the ShengBTE results by this factor:
	//	  double conversion = 1. / pow(13.60569300974785278,2)
	//	    * pow(0.52917721067000001,6)
	//	    * pow(0.00109715981930014,3);
}
