#include "gtest/gtest.h"
#include "points.h"
#include "qe_input_parser.h"

TEST (PhononH0, Velocity) {
//int main() {
	Context context;
	context.setPhD2FileName("../test/interaction3ph/QEspresso.fc");

	QEParser qeParser;
	auto [crystal,phononH0] = qeParser.parsePhHarmonic(context);
	phononH0.setAcousticSumRule("simple");

	// now, let's create a fine mesh

	Eigen::Vector3i qMesh;
	qMesh << 40, 40, 40;
	FullPoints points(crystal, qMesh);

	// pick a point close to gamma and get energies/velocities
	long ik = 1;
	auto p = points.getPoint(ik);
	auto [omega,z] = phononH0.diagonalize(p);
	auto v = phononH0.diagonalizeVelocity(p);

	// take out the group velocity
	long numBands = omega.size();
	Eigen::MatrixXd groupV(3,numBands);
	for ( int ib=0; ib<numBands; ib++ ) {
		for ( int i : {0,1,2} ) {
			groupV(i,ib) = v(ib,ib,i).real();
		}
	}

	// select acoustic modes close to gamma
	auto v0 = groupV.col(0);
	auto v1 = groupV.col(1);
	auto v2 = groupV.col(2);
	auto q = p.getCoords(Points::cartesianCoords);

	// for these three acoustic modes, check velocity is parallel to wavevector
	ASSERT_NEAR(v0.dot(q)/q.norm()/v0.norm(), 1., 0.04);
	ASSERT_NEAR(v1.dot(q)/q.norm()/v1.norm(), 1., 0.04);
	ASSERT_NEAR(v2.dot(q)/q.norm()/v2.norm(), 1., 0.04);

	// for silicon, the velocity is around 2200 m/s
	ASSERT_NEAR(abs(v0.minCoeff())*velocityRyToSi, 2200., 100.);

	// for another sanity check
	// we can also verify that, for acoustic phonons in silicon close to gamma,
	// the velocity is approximately (omega/q)we can approximate the velocity

	double err0 = abs( omega(0)/q.norm() - v0.norm() ) / v0.norm();
	double err1 = abs( omega(1)/q.norm() - v1.norm() ) / v1.norm();
	double err2 = abs( omega(2)/q.norm() - v2.norm() ) / v2.norm();
	// we allow a 4% error, (anisotropies...)
	ASSERT_NEAR(err0, 0., 0.04);
	ASSERT_NEAR(err1, 0., 0.04);
	ASSERT_NEAR(err2, 0., 0.04);
};
