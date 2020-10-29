#include "gtest/gtest.h"
#include "points.h"
#include "qe_input_parser.h"

TEST (PhononH0, Velocity) {
//int main() {
  Context context;
  context.setPhD2FileName("../test/interaction3ph/QEspresso.fc");
  context.setSumRuleD2("simple");

  QEParser qeParser;
  auto tup = qeParser.parsePhHarmonic(context);
  auto crystal = std::get<0>(tup);
  auto phononH0 = std::get<1>(tup);

  // now, let's create a fine mesh

  Eigen::Vector3i qMesh;
  qMesh << 40, 40, 40;
  FullPoints points(crystal, qMesh);

  // pick a point close to gamma and get energies/velocities
  long ik = 1;
  auto p = points.getPoint(ik);
  auto tup1 = phononH0.diagonalize(p);
  auto omega = std::get<0>(tup1);
  auto z = std::get<1>(tup1);
  auto v = phononH0.diagonalizeVelocity(p);

  // take out the group velocity
  long numBands = omega.size();
  Eigen::MatrixXd groupV(3, numBands);
  for (int ib = 0; ib < numBands; ib++) {
    for (int i : {0, 1, 2}) {
      groupV(i, ib) = v(ib, ib, i).real();
    }
  }

  // select acoustic modes close to gamma
  auto v0 = groupV.col(0);
  auto v1 = groupV.col(1);
  auto v2 = groupV.col(2);
  auto q = p.getCoords(Points::cartesianCoords);

  // for these three acoustic modes, check velocity is parallel to wavevector
  ASSERT_NEAR(v0.dot(q) / q.norm() / v0.norm(), 1., 0.04);
  ASSERT_NEAR(v1.dot(q) / q.norm() / v1.norm(), 1., 0.04);
  ASSERT_NEAR(v2.dot(q) / q.norm() / v2.norm(), 1., 0.04);

  // for silicon, the velocity is around 2200 m/s
  ASSERT_NEAR(abs(v0.minCoeff()) * velocityRyToSi, 2200., 100.);

  // for another sanity check
  // we can also verify that, for acoustic phonons in silicon close to gamma,
  // the velocity is approximately (omega/q)we can approximate the velocity

  double err0 = abs(omega(0) / q.norm() - v0.norm()) / v0.norm();
  double err1 = abs(omega(1) / q.norm() - v1.norm()) / v1.norm();
  double err2 = abs(omega(2) / q.norm() - v2.norm()) / v2.norm();
  // we allow a 4% error, (anisotropies...)
  ASSERT_NEAR(err0, 0., 0.04);
  ASSERT_NEAR(err1, 0., 0.04);
  ASSERT_NEAR(err2, 0., 0.04);
};


TEST (WannierH0, Velocity) {
//int main() {
  Context context;
  context.setElectronH0Name("../test/data/si_tb.dat");

  Eigen::MatrixXd atomicPositions(2,3);
  atomicPositions.row(0) << 0., 0., 0.;
  atomicPositions.row(1) << 1.34940, 1.34940, 1.34940;

  Eigen::VectorXi atomicSpecies(2);
  atomicSpecies(0) = 0;
  atomicSpecies(1) = 0;

  std::vector<std::string> speciesNames;
  speciesNames.push_back("Si");

  context.setInputAtomicPositions(atomicPositions);
  context.setInputAtomicSpecies(atomicSpecies);
  context.setInputSpeciesNames(speciesNames);

  auto tup = QEParser::parseElHarmonicWannier(context);
  auto crystal = std::get<0>(tup);
  auto electronH0 = std::get<1>(tup);

  Eigen::Vector3d k1;
  k1.setZero();
  auto v1 = electronH0.diagonalizeVelocityFromCoords(k1);

  Eigen::Vector3d k2;
  k2.setZero();
  double deltaK = 0.0005;

  k2(0) += deltaK;
  auto v2 = electronH0.diagonalizeVelocityFromCoords(k2);

  auto tup1 = electronH0.diagonalizeFromCoords(k1);
  auto ens1 = std::get<0>(tup1);
  auto tup2 = electronH0.diagonalizeFromCoords(k2);
  auto ens2 = std::get<0>(tup2);

  int ib = 3;

  double mass1 = 2*(ens2(ib)-ens1(ib)) / deltaK / deltaK;
  double mass2 = 2*v2(ib,ib,0).real() / deltaK;

  // for these three acoustic modes, check velocity is parallel to wavevector
  ASSERT_NEAR(abs((mass1-mass2)/mass1), 0., 0.02);

};
