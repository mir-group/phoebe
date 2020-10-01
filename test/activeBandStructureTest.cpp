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

  // TEST check that the number of bands are the same
  long sumBands = 0; 
  for (long point = 0; point < absOTF.getNumPoints(); point++) { 
    auto ikTemp = WavevectorIndex(point); 
    sumBands += absOTF.getNumBands(ikTemp) - absAPP.getNumBands(ikTemp);
  }
  ASSERT_EQ(sumBands, 0); 

  // pick a wavevector to check energies, eigenvectors, vels
  long ik = 7;
  auto ikIndex = WavevectorIndex(ik);

  // grab the points associated with this wavevector
  // TODO might be nice to check that they reference the same 
  // point 
  Point pointOTF = absOTF.getPoint(ik);
  //Point pointAPP = absAPP.getPoint(ik);

  // get the number of bands at this kpoint
  double nbOTF = absOTF.getNumBands(ikIndex);
  double nbAPP = absAPP.getNumBands(ikIndex);

  // generate values directly from the hamiltonian
  auto tup1 = phH0.diagonalize(pointOTF);
  auto ensT = std::get<0>(tup1);
  auto eigvecsT = std::get<1>(tup1);
  auto velsT = phH0.diagonalizeVelocity(pointOTF);

  // get the values stored in each active bandstructure
  auto ensOTF = absOTF.getEnergies(ikIndex);
  auto ensAPP = absAPP.getEnergies(ikIndex);

  Eigen::Tensor<std::complex<double>, 3> eigvecsOTF = absOTF.getPhEigenvectors(ikIndex);
  Eigen::Tensor<std::complex<double>, 3> eigvecsAPP = absAPP.getPhEigenvectors(ikIndex); 

  auto velsOTF = absOTF.getVelocities(ikIndex);
  auto velsAPP = absAPP.getVelocities(ikIndex);

  // TEST check OTF built bandstructure -----------------------

  // check the energies 
  double otfEns = (ensT - ensOTF).norm();
  ASSERT_EQ(otfEns, 0.);

  // check the velocities 
  std::complex<double> otfVels = complexZero;
  for (long ib1 = 0; ib1 < nbOTF; ib1++) {
    for (long ib2 = 0; ib2 < nbOTF; ib2++) {
      for (long ic = 0; ic < 3; ic++) {
        otfVels += pow(velsT(ib1, ib2, ic) - velsOTF(ib1, ib2, ic), 2);
      }
    }
  }
  ASSERT_EQ(otfVels, complexZero);
 
  // check the eigenvectors 
  std::complex<double> otfEigs = complexZero;
  for ( long i = 0; i<nbOTF; i++ ) {
        auto tup = decompress2Indeces(i,numAtoms,3);
        auto iat = std::get<0>(tup);
        auto ic = std::get<1>(tup);
        for ( long j = 0; j<nbOTF; j++ ) {
          otfEigs += pow(eigvecsT(i,j) - eigvecsOTF(ic, iat, j), 2);
    }
  }
  ASSERT_EQ(otfEigs, complexZero);

  // TEST check APP built bandstructure -----------------------

  // check the energies 
  double appEns = (ensT - ensAPP).norm();
  ASSERT_EQ(appEns, 0.);

  // check the velocities 
  std::complex<double> appVels = complexZero;
  for (long ib1 = 0; ib1 < nbAPP; ib1++) {
    for (long ib2 = 0; ib2 < nbAPP; ib2++) {
      for (long ic = 0; ic < 3; ic++) {
        appVels += pow(velsT(ib1, ib2, ic) - velsAPP(ib1, ib2, ic), 2);
      }
    }
  }
  ASSERT_EQ(appVels, complexZero);

  // check the eigenvectors 
  std::complex<double> appEigs = complexZero;
  for ( long i = 0; i<nbAPP; i++ ) {
        auto tup = decompress2Indeces(i,numAtoms,3);
        auto iat = std::get<0>(tup);
        auto ic = std::get<1>(tup);
        for ( long j = 0; j<nbAPP; j++ ) {
          appEigs += pow(eigvecsT(i,j) - eigvecsAPP(ic, iat, j), 2);
    }
  }
  ASSERT_EQ(appEigs, complexZero);

} 
