#include "gtest/gtest.h"
#include "points.h"
#include "state.h"
#include "ifc3_parser.h"
#include "qe_input_parser.h"
#include "ph_scattering.h"
#include "bandstructure.h"
#include <fstream>

TEST (Interaction3Ph, Coupling3Ph000) {
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

	DetachedState s1(q1, energies, numAtoms, numBands, ev1, nullptr);
	DetachedState s2(q2, energies, numAtoms, numBands, ev2, nullptr);
	DetachedState s3(q3, energies, numAtoms, numBands, ev3, nullptr);

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

//		for ( int i=0; i<numBands; i++ ) {
//			for ( int j=0; j<numBands; j++ ) {
//				for ( int k=0; k<numBands; k++ ) {
//					std::cout << i+1 << " " << j+1 << " " << k+1 << " "
//							<< couplingPlus(i,j,k)
//							<< " "
//							<< couplingMins(i,j,k)
//							<< " "
//							<< referenceCoupling(i,j,k)
//							<< "\n";
//				}
//			}
//			std::cout << "\n";
//		}

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

TEST (Interaction3Ph, Coupling3Ph210) {
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
	context.setPhD2FileName("./test/interaction3ph/QEspresso.fc");
	context.setPhD3FileName("./test/interaction3ph/ShengBTEForceConstants3rd");

	QEParser qeParser;
	auto [crystal,phononH0] = qeParser.parsePhHarmonic(context);
	phononH0.setAcousticSumRule("simple");

	IFC3Parser ifc3Parser;
	auto coupling3Ph = ifc3Parser.parseFromShengBTE(context, crystal);

	//------Set up cubic silicon-------//
	//Total number of wave vectors
	long nq = 20*20*20;
	//Number of atoms
	long numAtoms = crystal.getNumAtoms();
	//Number of bands
	long numBands = 3*numAtoms;

	auto atomicMasses = crystal.getAtomicMasses();

	//------------end setup-----------//

	Eigen::Vector3d q1, q2, q3;
	Eigen::MatrixXcd evm1(numBands,numBands),
			evm2(numBands,numBands), evm3(numBands,numBands);

	// Form a triplet to test vertex calculator
	long iq1 = 10;
	long iq2 = 210;
	long iq3 = 200;

	if ( false ) {
		// left here in case we need debugging in the future
		// if we use eigenvectors from file, the results are equal to the
		// reference values (even for degenerate points)

		// Read full BZ wave vectors (Cartesian) from file
		Eigen::MatrixXd q(nq,3);
		std::ifstream qfile("./test/interaction3ph/qpoints_full_cart");
		for ( long iq = 0; iq < nq; iq++ ) {
			for ( long ip = 0; ip < 3; ip++ ) {
				double x;
				qfile >> x;
				q(iq,ip) = x * distanceBohrToAng / 10.;
				// convert from nm to bohr
			}
		}
		qfile.close();

		// Read full BZ eigenvectors from file
		// Note that this is read as outputted by ShengBTE
		// and has to be reshaped to conform to PHOEBE's format
		// of the eigenvector.
		std::ifstream evfile("./test/interaction3ph/evecs_full");
		Eigen::Tensor<complex<double>,3> ev(nq,numBands,numBands);
		for(long iq = 0; iq < nq; iq++){
			double re,im;
			for(long ib = 0; ib < numBands; ib++){
				auto [idim,ia] = decompress2Indeces(ib,3,numAtoms);
				for(long jb = 0; jb < numBands; jb++){
					evfile >> re >> im;
					ev(iq,ib,jb) = {re / sqrt(atomicMasses(ia)),
							im / sqrt(atomicMasses(ia))};
				}
			}
		}
		evfile.close();

		q1 = q.row(iq1);
		q2 = q.row(iq2);
		q3 = q.row(iq3);

		// Reshape the eigenvectors read from file
		Eigen::Tensor<complex<double>,3> ev1(3,numAtoms,numBands);
		Eigen::Tensor<complex<double>,3> ev2(3,numAtoms,numBands);
		Eigen::Tensor<complex<double>,3> ev3(3,numAtoms,numBands);

		for ( long idim = 0; idim < 3; idim++ ) {
			for ( long iat = 0; iat < numAtoms; iat++ ) {
				for ( long ib = 0; ib < numBands; ib++ ) {
					ev1(idim,iat,ib) = ev(iq1,idim+3*iat,ib);
					ev2(idim,iat,ib) = ev(iq2,idim+3*iat,ib);
					ev3(idim,iat,ib) = ev(iq3,idim+3*iat,ib);
				}
			}
		}
		for ( int i=0; i<numBands; i++ ) {
			for ( int j=0; j<numBands; j++ ) {
				auto [iat,idim] = decompress2Indeces(i,numAtoms,3);
				evm1(i,j) = ev1(idim,iat,j);
				evm2(i,j) = ev2(idim,iat,j);
				evm3(i,j) = ev3(idim,iat,j);
			}
		}

	} else {

		Eigen::Vector3i qMesh;
		qMesh << 20, 20, 20;
		FullPoints points(crystal, qMesh);
		auto p1 = points.getPoint(iq1);
		auto p2 = points.getPoint(iq2);
		auto p3 = points.getPoint(iq3);
		auto [energies1,evv1] = phononH0.diagonalize(p1);
		auto [energies2,evv2] = phononH0.diagonalize(p2);
		auto [energies3,evv3] = phononH0.diagonalize(p3);

		for ( int i=0; i<numBands; i++ ) {
			for ( int j=0; j<numBands; j++ ) {
				auto [iat,idim] = decompress2Indeces(i,numAtoms,3);
				evm1(i,j) = evv1(idim,iat,j);
				evm2(i,j) = evv2(idim,iat,j);
				evm3(i,j) = evv3(idim,iat,j);
			}
		}
		q1 = p1.getCoords(Points::cartesianCoords);
		q2 = p2.getCoords(Points::cartesianCoords);
		q3 = p3.getCoords(Points::cartesianCoords);
	}

	// note: the reference was generated without the normalization by energies
	// so we set them to one.
	Eigen::VectorXd energies(numBands);
	energies.setConstant(1.);

	DetachedState s1(q1, energies, numAtoms, numBands, evm1, nullptr);
	DetachedState s2(q2, energies, numAtoms, numBands, evm2, nullptr);
	DetachedState s3(q3, energies, numAtoms, numBands, evm3, nullptr);

	auto [couplingPlus,couplingMins] = coupling3Ph.getCouplingSquared(
														s1, s2, s3, s3);

	// we load reference data

	Eigen::Tensor<double,3> referenceCoupling(numBands,numBands,numBands);
	referenceCoupling.setZero();
	std::ifstream tfile("./test/interaction3ph/referenceCMins10-210-200");
	int i_, j_, k_;
    if (tfile.is_open()) {
		double x1, x2;
		for ( long i=0; i<numBands; i++ ) {
			for ( long j=0; j<numBands; j++ ) {
				for ( long k=0; k<numBands; k++ ) {
					tfile >> i_ >> j_ >> k_ >> x1 >> x2;
					referenceCoupling(i,j,k) = x1/ pow(13.60569300974785278,2)
    	    * pow(0.52917721067000001,6)
    	    * pow(0.00109715981930014,3);
				}
			}
		}
    }
	tfile.close();

	// now we compare the difference

	double x1, x2, x3;
	x1 = 0.;
	x2 = 0.;
	x3 = 0.;
	for ( int i=0; i<2; i++ ) {
		for ( int j=0; j<numBands; j++ ) {
			for ( int k=0; k<numBands; k++ ) {
				x1 += couplingPlus(i,j,k);
				x2 += couplingMins(i,j,k);
				x3 += referenceCoupling(i,j,k);
			}
		}
	}

	double relativeErrorP = abs((x1-x3)/x1);
	double relativeErrorM = abs((x1-x3)/x1);

	ASSERT_NEAR(relativeErrorP, 0., 1.0e-3);
	ASSERT_NEAR(relativeErrorM, 0., 1.0e-3);

	// now we test the same configuration, but we use State<FullPoints>
	// rather than DetachedState, to generate the coupling

	Eigen::Vector3i qMesh;
	qMesh << 20, 20, 20;
	FullPoints points(crystal, qMesh);
	bool withVelocities = true;
	bool withEigenvectors = true;
	FullBandStructure<FullPoints> bandStructure = phononH0.populate(points,
			withVelocities, withEigenvectors);

	auto states1 = bandStructure.getState(iq1);
	auto states2 = bandStructure.getState(iq2);
	auto states3Plus = bandStructure.getState(iq3);
	auto states3Mins = bandStructure.getState(iq3);

	auto p1 = states1.getPoint();
	auto p2 = states2.getPoint();
	auto p3PlusTest = p1 + p2;
	auto p3MinsTest = p1 - p2;
	auto p3Plus = states3Plus.getPoint();
	auto p3Mins = states3Mins.getPoint();

	// check that the sum of Point works
	ASSERT_EQ((p3PlusTest.getCoords(Points::cartesianCoords)-p3Plus.getCoords(Points::cartesianCoords)).norm(), 0.);
	ASSERT_EQ((p3MinsTest.getCoords(Points::cartesianCoords)-p3Mins.getCoords(Points::cartesianCoords)).norm(), 0.);

	auto [couplingPlus2,couplingMins2] = coupling3Ph.getCouplingSquared(
			states1, states2, states3Plus, states3Mins);

	auto en1 = states1.getEnergies();
	auto en2 = states2.getEnergies();
	auto en3Plus = states3Plus.getEnergies();
	auto en3Mins = states3Mins.getEnergies();

	x1 = 0.;
	x2 = 0.;
	for ( int i=0; i<2; i++ ) {
		for ( int j=0; j<numBands; j++ ) {
			for ( int k=0; k<numBands; k++ ) {
				x1 += pow(couplingPlus(i,j,k) -
						couplingPlus2(i,j,k),2);
				x2 += pow(couplingMins(i,j,k) -
						couplingMins2(i,j,k),2);
			}
		}
	}
	ASSERT_NEAR(x1, 0., 1.0e-40);
	ASSERT_NEAR(x2, 0., 1.0e-40);
}
