#include "crystal.h"
#include "qe_input_parser.h"
#include "gtest/gtest.h"

TEST(Crystal, Test1) {

  // test the crystal class for silicon
  Eigen::Matrix3d directUnitCell;
  directUnitCell.row(0) << -5.1, 0., 5.1;
  directUnitCell.row(1) << 0., 5.1, 5.1;
  directUnitCell.row(2) << -5.1, 5.1, 0.;
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

  // set up the crystal object
  Crystal crystal(context, directUnitCell, atomicPositions, atomicSpecies,
                  speciesNames, speciesMasses, dimensionality);

  EXPECT_EQ(crystal.getNumAtoms(), 2);
  EXPECT_EQ(crystal.getNumSpecies(), 1);

  int numSym = crystal.getNumSymmetries();
  EXPECT_EQ(numSym, 48);


}
