#include "gtest/gtest.h"
#include "points.h"
#include "qe_input_parser.h"

TEST (PointsTest, PointsHandling) {
  Eigen::Matrix3d directUnitCell;
  directUnitCell(0,0) = 0.;
  directUnitCell(0,1) = 0.;
  directUnitCell(0,2) = 5.1;;
  directUnitCell(1,0) = 0.;
  directUnitCell(1,1) = 5.1;
  directUnitCell(1,2) = 5.1;
  directUnitCell(2,0) = -5.1;
  directUnitCell(2,1) = 5.1;
  directUnitCell(2,2) = 0.;
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
}
