#include "gtest/gtest.h"
#include "points.h"
#include "qe_input_parser.h"

TEST (PointsTest, PointsHandling) {
	Eigen::Matrix3d directUnitCell;
	directUnitCell.row(0) << 0., 0., 5.1;
	directUnitCell.row(1) << 0., 5.1, 5.1;
	directUnitCell.row(2) << -5.1, 5.1, 0.;
	Eigen::MatrixXd atomicPositions(2,3);
	atomicPositions.row(0) << 0., 0., 0.;
	atomicPositions.row(1) << 2.55, 2.55, 2.55;
	Eigen::VectorXi atomicSpecies(2);
	atomicSpecies << 0, 0;
	std::vector<std::string> speciesNames;
	speciesNames.push_back("Si");
	Eigen::VectorXd speciesMasses(1);
	speciesMasses(0) = 28.086;
	long dimensionality = 3;

	Crystal crystal(directUnitCell, atomicPositions, atomicSpecies,
			speciesNames, speciesMasses, dimensionality);

	Eigen::Vector3i mesh;
	mesh << 4,4,4;
	FullPoints points(crystal, mesh);

	auto [mesh_,offset_] = points.getMesh();

	EXPECT_EQ ((mesh-mesh_).norm(),0.);

	auto p1 = points.getPoint(4);

	// find the index of the inverted point
	long i4 = points.getIndex(-p1.getCoords(Points::cartesianCoords));
//	long i4 = points.getIndexInverted(4);

	auto p2 = points.getPoint(i4);
	auto p3 = p1 + p2;
	EXPECT_EQ ( p3.getCoords(Points::cartesianCoords).norm(),0.);

	//

	mesh << 2,2,2;
	points = FullPoints(crystal, mesh);
	long iq = 7;
	p1 = points.getPoint(iq);
//	long iqr = points.getIndexInverted(iq);
	long iqr = points.getIndex(-p1.getCoords(Points::cartesianCoords));
	p2 = points.getPoint(iqr);
	p3 = p1 + p2;
	EXPECT_EQ ( p3.getCoords().norm(),0.);


}
