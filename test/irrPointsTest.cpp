#include "full_points.h"
#include "gtest/gtest.h"

TEST(IrrPointsTest, Symmetries) {
  Eigen::Matrix3d directUnitCell;
  directUnitCell.col(0) << -5.1, 0., 5.1;
  directUnitCell.col(1) << 0., 5.1, 5.1;
  directUnitCell.col(2) << -5.1, 5.1, 0.;
  Eigen::MatrixXd atomicPositions(2, 3);
  atomicPositions.row(0) << 0., 0., 0.;
  atomicPositions.row(1) << 2.55, 2.55, 2.55;
  Eigen::VectorXi atomicSpecies(2);
  atomicSpecies << 0, 0;
  std::vector<std::string> speciesNames;
  speciesNames.push_back("Si");
  Eigen::VectorXd speciesMasses(1);
  speciesMasses(0) = 28.086;
  long dimensionality = 3;

  Context context;
  context.setUseSymmetries(true);

  Crystal crystal(context, directUnitCell, atomicPositions, atomicSpecies,
                  speciesNames, speciesMasses, dimensionality);

  ASSERT_EQ(crystal.getNumSymmetries(),48);

  //-----------------

  Eigen::Vector3i mesh;
  mesh << 4, 4, 4;
  FullPoints fullPoints(crystal, mesh);
  fullPoints.setIrreduciblePoints();

  //-----------------------------------

  // I hard code that I expect 8 irreducible points
  long numIrrPoints = fullPoints.irrPointsIterator().size();
  ASSERT_EQ(numIrrPoints,8);

  int numFullPoints = 0;
  for (long ik : fullPoints.irrPointsIterator()) {
    numFullPoints += fullPoints.getRotationsStar(ik).size();
  }
  ASSERT_EQ(numFullPoints,fullPoints.getNumPoints());

  // here I check the symmetries matrices
  // loop over irreducible points, unfold the star, and check that we can
  // reconstruct the whole list of points of the full grid
  int counter = 0;
  std::vector<long> allIndeces;
  for (long ikIrr : fullPoints.irrPointsIterator()) {
    auto kIrr = fullPoints.getPointCoords(ikIrr, Points::cartesianCoords);

    long ikIrrAsRed = fullPoints.asIrreducibleIndex(ikIrr);
    ASSERT_EQ(ikIrrAsRed, counter);
    counter++;

    auto rots = fullPoints.getRotationsStar(ikIrr);
    for ( auto s : rots ) {
      auto kRedCart = s * kIrr; // in cartesian coordinates
      auto kRedCrys = fullPoints.cartesianToCrystal(kRedCart);
      long oldIndex = fullPoints.getIndex(kRedCrys); // getIndex needs crystal coords
      allIndeces.push_back(oldIndex);

      auto t = fullPoints.getRotationToIrreducible(kRedCart, Points::cartesianCoords);
      long ik2 = std::get<0>(t);
      ASSERT_EQ(ik2,ikIrr);

      Eigen::Matrix3d rot = std::get<1>(t);
      // this rotation should map the point to the irreducible one
      Eigen::Vector3d kIrr2 = rot * kRedCart;

      ASSERT_NEAR((rot.inverse()-s).squaredNorm(),0.,0.0001);

      Eigen::Vector3d x1 = fullPoints.cartesianToCrystal(kIrr2).transpose();
      Eigen::Vector3d x2 = fullPoints.cartesianToCrystal(kIrr).transpose();
      ASSERT_NEAR((x1-x2).squaredNorm(),0.,0.0001);
    }
  }
  int uniqueCount = std::unique(allIndeces.begin(), allIndeces.end()) - allIndeces.begin();
  ASSERT_EQ(mesh.prod(),uniqueCount);
}
