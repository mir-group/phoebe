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
	auto tup = qeParser.parsePhHarmonic(context);
 auto crystal = std::get<0>(tup);
 auto phononH0 = std::get<1>(tup);

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

  auto tup1 = phononH0.diagonalizeFromCoords(q1);
  auto energies1 = std::get<0>(tup1);
  auto ev1 = std::get<1>(tup1);
  auto tup2 = phononH0.diagonalizeFromCoords(q2);
  auto energies2 = std::get<0>(tup2);
  auto ev2 = std::get<1>(tup2);
  auto tup3 = phononH0.diagonalizeFromCoords(q3);
  auto energies3 = std::get<0>(tup3);
  auto ev3 = std::get<1>(tup3);

	DetachedState s1(q1, energies, numBands, numBands, ev1, nullptr);
	DetachedState s2(q2, energies, numBands, numBands, ev2, nullptr);
	DetachedState s3(q3, energies, numBands, numBands, ev3, nullptr);

  int nb1 = energies.size();
  int nb2 = energies.size();

  std::vector<Eigen::Vector3d> q1s_e(1);
  Eigen::Vector3d q2_e;
  std::vector<Eigen::MatrixXcd> ev1s_e(1);
  Eigen::MatrixXcd ev2_e;
  std::vector<Eigen::MatrixXcd> ev3Pluss_e(1);
  std::vector<Eigen::MatrixXcd> ev3Minss_e(1);
  std::vector<int> nb1s_e(1);
  std::vector<int> nb3Pluss_e(1);
  std::vector<int> nb3Minss_e(1);

  q1s_e[0] = s1.getCoords(Points::cartesianCoords);
  nb1s_e[0] = nb1;
  nb3Pluss_e[0] = energies.size();
  nb3Minss_e[0] = energies.size();
  s1.getEigenvectors(ev1s_e[0]);
  s3.getEigenvectors(ev3Pluss_e[0]);
  s3.getEigenvectors(ev3Minss_e[0]);
  s2.getEigenvectors(ev2_e);
  q2_e = q2;

  auto tup8 = coupling3Ph.getCouplingsSquared(q1s_e, q2_e, ev1s_e, ev2_e,
           ev3Pluss_e, ev3Minss_e, nb1s_e, nb2, nb3Pluss_e, nb3Minss_e);
  auto couplingPlus = std::get<0>(tup8)[0];
  auto couplingMins = std::get<1>(tup8)[0];

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
	auto tup = qeParser.parsePhHarmonic(context);
 auto crystal = std::get<0>(tup);
 auto phononH0 = std::get<1>(tup);

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
	auto tup1 = phononH0.diagonalize(p1);
 auto energies1 = std::get<0>(tup1);
 auto evm1 = std::get<1>(tup1);
	auto tup2 = phononH0.diagonalize(p2);
 auto energies2 = std::get<0>(tup2);
 auto evm2 = std::get<1>(tup2);
	auto tup3 = phononH0.diagonalize(p3);
 auto energies3 = std::get<0>(tup3);
 auto evm3 = std::get<1>(tup3);
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

  int nb1 = energies.size();
  int nb2 = energies.size();

  std::vector<Eigen::Vector3d> q1s_e(1);
  Eigen::Vector3d q2_e;
  std::vector<Eigen::MatrixXcd> ev1s_e(1);
  Eigen::MatrixXcd ev2_e;
  std::vector<Eigen::MatrixXcd> ev3Pluss_e(1);
  std::vector<Eigen::MatrixXcd> ev3Minss_e(1);
  std::vector<int> nb1s_e(1);
  std::vector<int> nb3Pluss_e(1);
  std::vector<int> nb3Minss_e(1);

  q1s_e[0] = s1.getCoords(Points::cartesianCoords);
  nb1s_e[0] = nb1;
  nb3Pluss_e[0] = energies.size();
  nb3Minss_e[0] = energies.size();
  s1.getEigenvectors(ev1s_e[0]);
  s3.getEigenvectors(ev3Pluss_e[0]);
  s3.getEigenvectors(ev3Minss_e[0]);
  s2.getEigenvectors(ev2_e);
  q2_e = q2;

  auto tup7 = coupling3Ph.getCouplingsSquared(q1s_e, q2_e, ev1s_e, ev2_e,
      ev3Pluss_e, ev3Minss_e, nb1s_e, nb2, nb3Pluss_e, nb3Minss_e);
  auto couplingPlus = std::get<0>(tup7)[0];
  auto couplingMins = std::get<1>(tup7)[0];

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

	auto en1 = states1.getEnergies();
	auto en2 = states2.getEnergies();
	auto en3Plus = states3Plus.getEnergies();
	auto en3Mins = states3Mins.getEnergies();

  nb1 = en1.size();
  nb2 = en2.size();
  q2 = states2.getCoords(Points::cartesianCoords);

  q1s_e[0] = s1.getCoords(Points::cartesianCoords);
  nb1s_e[0] = nb1;
  nb3Pluss_e[0] = en3Plus.size();
  nb3Minss_e[0] = en3Mins.size();
  states1.getEigenvectors(ev1s_e[0]);
  states3Plus.getEigenvectors(ev3Pluss_e[0]);
  states3Mins.getEigenvectors(ev3Minss_e[0]);
  states2.getEigenvectors(ev2_e);
  q2_e = q2;
  auto tup6 = coupling3Ph.getCouplingsSquared(q1s_e, q2_e, ev1s_e, ev2_e,
      ev3Pluss_e, ev3Minss_e, nb1s_e, nb2, nb3Pluss_e, nb3Minss_e);
  auto couplingPlus2 = std::get<0>(tup6)[0];
  auto couplingMins2 = std::get<1>(tup6)[0];


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
