//#include "gtest/gtest.h"
#include "points.h"
#include "state.h"
#include "ifc3_parser.h"
#include "qe_input_parser.h"
#include "ph_scattering.h"
#include "bandstructure.h"

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

int main() {
	Context context;
	context.setPhD2FileName("interaction3ph/QEspresso.fc");
	context.setPhD3FileName("interaction3ph/ShengBTEForceConstants3rd");
	Eigen::VectorXd temperatures(1);
	temperatures << 300.;
	context.setTemperatures(temperatures);
	context.setSmearingMethod(2); // tetrahedron smearing

	QEParser qeParser;
	auto [crystal,phononH0] = qeParser.parsePhHarmonic(context);

	IFC3Parser ifc3Parser;
	auto coupling3Ph = ifc3Parser.parseFromShengBTE(context, crystal);

	// set up the mesh of qpoints of shengbte
	Eigen::Vector3i qMesh;
	qMesh << 2, 2, 2;
	FullPoints points(crystal, qMesh);
	bool withVelocities = true;
	bool withEigenvectors = true;
	FullBandStructure<FullPoints> bandStructure = phononH0.populate(points,
			withVelocities, withEigenvectors);

	long iq1 = 2;
	long iq2 = 7;

	auto q1 = bandStructure.getPoint(iq1);
	auto q2 = bandStructure.getPoint(iq2);

	auto states1 = bandStructure.getState(iq1);

	long iq2Inv = bandStructure.getPoints().getIndexInverted(iq2);

	auto states2 = bandStructure.getState(iq2);
	auto states2Plus = bandStructure.getState(iq2Inv);

	auto q3Plus = q1 + q2;
	auto q3Mins = q1 - q2;

	auto states3Plus = bandStructure.getState(q3Plus);
	auto states3Mins = bandStructure.getState(q3Mins);

	auto [couplingPlus, couplingMins] = coupling3Ph.getCouplingSquared(
			states1, states2Plus, states2,
			states3Plus, states3Mins);
}
