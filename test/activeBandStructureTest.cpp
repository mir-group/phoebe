#include "active_bandstructure.h"
#include "ifc3_parser.h"
#include "ph_scattering.h"
#include "points.h"
#include "qe_input_parser.h"
#include <fstream>
#include <gtest/gtest.h>

/* This test checks that the active bandstructure construction
 * for a phonon Hamiltonian is generated the same way regardless 
 * of if it was made via buildOnTheFly or buildAsPostProcessing.
 *
 * We use the phonon case, because the two functions should behave
 * the same when they are done using the same chemical potential, 
 * and for phonons, of course, the chemical potential is zero. 
 */
TEST(ActiveBandStructureTest,BandStructureStorage) {

  // set up a phononH0
  Context context;
  context.setPhD2FileName("../test/interaction3ph/QEspresso.fc");
  context.setupFromInput("../test/data/activeBands.in");
  QEParser qeParser;
  auto tup = qeParser.parsePhHarmonic(context);
  auto crystal = std::get<0>(tup);
  auto phH0 = std::get<1>(tup);

  // Number of atoms
  long numAtoms = crystal.getNumAtoms();
  // Number of bands
  long numBands = 3 * numAtoms;

  //------------end setup-----------//

  // generate an active bandstructure for phonons
  Eigen::Vector3i qMesh;
  qMesh << 10, 10, 10;
  FullPoints points(crystal, qMesh);

  bool withVelocities = true;
  bool withEigenvectors = true;

  // set up the window for filtering active points
  Eigen::VectorXd temperatures = context.getTemperatures();
  double temperatureMin = temperatures.minCoeff();
  double temperatureMax = temperatures.maxCoeff();
  Particle particle = phH0.getParticle();
  Window window(context, particle, temperatureMin, temperatureMax);

  // set up empty bandstructures
  ActivePoints bogusPoints(points, Eigen::VectorXi::Zero(1));
  ActiveBandStructure absOTF(particle, bogusPoints);
  ActiveBandStructure absAPP(particle, bogusPoints);

  // create two bandstructures, built with different methods
  absOTF.buildOnTheFly(window, points, phH0, withEigenvectors, withVelocities);
  absAPP.buildAsPostprocessing(context, points, phH0, withEigenvectors, withVelocities);

  // TEST check that they selected the same number of states
  ASSERT_EQ(absOTF.getNumPoints(),absAPP.getNumPoints()); 

  // pick a wavevector to check 
  long ik = 7;
  auto ikIndex = WavevectorIndex(ik);

  // grab the points associated with this wavevector
  Point pointOTF = absOTF.getPoint(ik);
  Point pointAPP = absAPP.getPoint(ik);

  // generate values directly from the hamiltonian
  auto tup1 = phH0.diagonalize(pointOTF);
  auto ensT = std::get<0>(tup1);
  auto eigvecsT = std::get<1>(tup1);
  auto velsT = phH0.diagonalizeVelocity(pointOTF);

  // get the values stored in each active bandstructure
  auto ensOTF = absOTF.getEnergies(ikIndex);
  auto ensAPP = absAPP.getEnergies(ikIndex);

  //Eigen::Tensor<std::complex<double>, 3> eigvecs = bandStructure.getPhEigenvectors(ikIndex);
  //auto vels = bandStructure.getVelocities(ikIndex);

  // now we check the difference between each bandstructure and the H0 values
  // TEST check OTF built bandstructure
  double x1 = (ensOTF - ensT).norm();
  ASSERT_EQ(x1, 0.);
 
  // TEST check APP built bandstructure
  x1 = (ensAPP - ensT).norm();
  ASSERT_EQ(x1, 0.);
} 
