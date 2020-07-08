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
	context.setPhD2FileName("../test/interaction3ph/QEspresso.fc");
	context.setPhD3FileName("../test/interaction3ph/ShengBTEForceConstants3rd");
	context.setSumRuleD2("simple");

	QEParser qeParser;
	auto [crystal,phononH0] = qeParser.parsePhHarmonic(context);

	IFC3Parser ifc3Parser;
	auto coupling3Ph = ifc3Parser.parseFromShengBTE(context, crystal);

	//Number of atoms
	long numAtoms = crystal.getNumAtoms();
	//Number of bands
	long numBands = 3*numAtoms;

	auto atomicMasses = crystal.getAtomicMasses();

	// Form a triplet to test vertex calculator
	Eigen::VectorXd energies(numBands);
	energies.setConstant(1.);

	Eigen::Vector3d q1, q2, q3;
	q1.setZero();
	q2.setZero();
	q3.setZero();

	auto [energies1,ev1] = phononH0.diagonalizeFromCoords(q1);
	auto [energies2,ev2] = phononH0.diagonalizeFromCoords(q2);
	auto [energies3,ev3] = phononH0.diagonalizeFromCoords(q3);

	DetachedState s1(q1, energies, numBands, numBands, ev1, nullptr);
	DetachedState s2(q2, energies, numBands, numBands, ev2, nullptr);
	DetachedState s3(q3, energies, numBands, numBands, ev3, nullptr);

	auto [couplingPlus,couplingMins] = coupling3Ph.getCouplingSquared(
            <#initializer #>, <#initializer #>, <#initializer #>,
            <#initializer #>, <#initializer #>, <#initializer #>);

	// we load reference data

	Eigen::Tensor<double,3> referenceCoupling(numBands,numBands,numBands);
	referenceCoupling.setZero();
	std::ifstream tfile("../test/interaction3ph/referenceCmins000");
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

	ASSERT_NEAR(relativeError, 0., 1.0e-4);

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
	context.setPhD2FileName("../test/interaction3ph/QEspresso.fc");
	context.setPhD3FileName("../test/interaction3ph/ShengBTEForceConstants3rd");
	context.setSumRuleD2("simple");

	QEParser qeParser;
	auto [crystal,phononH0] = qeParser.parsePhHarmonic(context);

	IFC3Parser ifc3Parser;
	auto coupling3Ph = ifc3Parser.parseFromShengBTE(context, crystal);

	//Number of atoms
	long numAtoms = crystal.getNumAtoms();
	//Number of bands
	long numBands = 3*numAtoms;

	auto atomicMasses = crystal.getAtomicMasses();

	// Form a triplet to test vertex calculator
	long iq1 = 10;
	long iq2 = 210;
	long iq3 = 200;

	Eigen::Vector3i qMesh;
	qMesh << 20, 20, 20;
	FullPoints points(crystal, qMesh);
	auto p1 = points.getPoint(iq1);
	auto p2 = points.getPoint(iq2);
	auto p3 = points.getPoint(iq3);
	auto [energies1,evm1] = phononH0.diagonalize(p1);
	auto [energies2,evm2] = phononH0.diagonalize(p2);
	auto [energies3,evm3] = phononH0.diagonalize(p3);
	auto q1 = p1.getCoords(Points::cartesianCoords);
	auto q2 = p2.getCoords(Points::cartesianCoords);
	auto q3 = p3.getCoords(Points::cartesianCoords);

	// note: the reference was generated without the normalization by energies
	// so we set them to one.
	Eigen::VectorXd energies(numBands);
	energies.setConstant(1.);

	DetachedState s1(q1, energies, numBands, numBands, evm1, nullptr);
	DetachedState s2(q2, energies, numBands, numBands, evm2, nullptr);
	DetachedState s3(q3, energies, numBands, numBands, evm3, nullptr);

	auto [couplingPlus,couplingMins] = coupling3Ph.getCouplingSquared(
            <#initializer #>, <#initializer #>, <#initializer #>,
            <#initializer #>, <#initializer #>, <#initializer #>);
	// we load reference data

	Eigen::Tensor<double,3> referenceCoupling(numBands,numBands,numBands);
	referenceCoupling.setZero();
	std::ifstream tfile("../test/interaction3ph/referenceCMins10-210-200");
	int i_, j_, k_;
    if (tfile.is_open()) {
		double x1, x2;
		for ( long i=0; i<numBands; i++ ) {
			for ( long j=0; j<numBands; j++ ) {
				for ( long k=0; k<numBands; k++ ) {
					tfile >> i_ >> j_ >> k_ >> x1 >> x2;
					referenceCoupling(i,j,k) = x1 / pow(13.60569300974785278,2)
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

	bool withVelocities = true;
	bool withEigenvectors = true;
	FullBandStructure bandStructure = phononH0.populate(points,
			withVelocities, withEigenvectors);

	auto states1 = bandStructure.getState(iq1);
	auto states2 = bandStructure.getState(iq2);
	auto states3Plus = bandStructure.getState(iq3);
	auto states3Mins = bandStructure.getState(iq3);

	p1 = states1.getPoint();
	p2 = states2.getPoint();
	auto p3PlusTest = p1 + p2;
	auto p3MinsTest = p1 - p2;
	auto p3Plus = states3Plus.getPoint();
	auto p3Mins = states3Mins.getPoint();

	// check that the sum of Point works
	ASSERT_EQ((p3PlusTest.getCoords(Points::cartesianCoords)-p3Plus.getCoords(Points::cartesianCoords)).norm(), 0.);
	ASSERT_EQ((p3MinsTest.getCoords(Points::cartesianCoords)-p3Mins.getCoords(Points::cartesianCoords)).norm(), 0.);

	auto [couplingPlus2,couplingMins2] = coupling3Ph.getCouplingSquared(
            <#initializer #>, <#initializer #>, <#initializer #>,
            <#initializer #>, <#initializer #>, <#initializer #>);

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
