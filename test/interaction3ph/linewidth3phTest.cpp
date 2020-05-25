#include "gtest/gtest.h"
#include "points.h"
#include "state.h"
#include "ifc3_parser.h"
#include "qe_input_parser.h"
#include "ph_scattering.h"
#include "bandstructure.h"

TEST (PhScatteringMatrix, Linewidth) {
	// in this test, we compute the 3ph interaction |V|^2 for three arbitrary
	// points, and compare the results with the coupling computed by ShengBTE

	Context context;
	context.setPhD2FileName("interaction3ph/QEspresso.fc");
	context.setPhD3FileName("itneraction3ph/ShengBTEForceConstants3rd");
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
	qMesh << 8, 8, 8;
	long numPoints = 8*8*8;

	FullPoints innerPoints(crystal, qMesh);

	Eigen::Vector3i outerMesh;
	outerMesh << 1, 1, 1;
	Eigen::Vector3d outerOffset;
	outerOffset << 1/8., 0., 0.;

	FullPoints outerPoints(crystal, outerMesh, outerOffset);

	bool withVelocities = true;
	bool withEigenvectors = true;
	FullBandStructure<FullPoints> outerBandStructure = phononH0.populate(
			outerPoints, withVelocities, withEigenvectors);
	FullBandStructure<FullPoints> innerBandStructure = phononH0.populate(
			innerPoints, withVelocities, withEigenvectors);

	// set the chemical potentials to zero, load temperatures
	StatisticsSweep statisticsSweep(context);

	// build/initialize the scattering matrix and the smearing
	PhScatteringMatrix scatteringMatrix(context, statisticsSweep,
			innerBandStructure, outerBandStructure, &coupling3Ph);
	// compute linewidths
	scatteringMatrix.setup();

	// load the linewidth
	VectorBTE sMatrixDiagonal = scatteringMatrix.diagonal();
	std::cout << sMatrixDiagonal.data(0,0) << "\t" <<
			sMatrixDiagonal.data(0,1) << "\t" <<
			sMatrixDiagonal.data(0,2) << "\t" <<
			sMatrixDiagonal.data(0,3) << "\t" <<
			sMatrixDiagonal.data(0,4) << "\t" <<
			sMatrixDiagonal.data(0,5);

	EXPECT_EQ ( sMatrixDiagonal.data(0,5), 0.);
}
