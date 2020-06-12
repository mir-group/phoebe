#include "bandstructure.h"
#include "ifc3_parser.h"
#include "ph_scattering.h"
#include "points.h"
#include "qe_input_parser.h"
#include "state.h"
#include <fstream>
#include <gtest/gtest.h>

TEST(FullBandStructureTest, BandStructureStorage) {
  // in this test, we check whether we are storing data inside the
  // bandstructure object correctly and consistently

  Context context;
  context.setPhD2FileName("../test/interaction3ph/QEspresso.fc");

  QEParser qeParser;
  auto [crystal, phononH0] = qeParser.parsePhHarmonic(context);

  // Number of atoms
  long numAtoms = crystal.getNumAtoms();
  // Number of bands
  long numBands = 3 * numAtoms;

  //------------end setup-----------//

  Eigen::Vector3i qMesh;
  qMesh << 10, 10, 10;
  FullPoints points(crystal, qMesh);

  bool withVelocities = true;
  bool withEigenvectors = true;
  FullBandStructure bandStructure =
      phononH0.populate(points, withVelocities, withEigenvectors);

  long ik = 7;
  Point point = bandStructure.getPoint(ik);
  auto [ensT, eigvecsT] = phononH0.diagonalize(point);
  auto velsT = phononH0.diagonalizeVelocity(point);

  auto s = bandStructure.getState(ik);
  auto ens = s.getEnergies();
  Eigen::Tensor<std::complex<double>, 3> eigvecs(3, numAtoms, numBands);
  s.getEigenvectors(eigvecs);
  auto vels = s.getVelocities();

  // now we check the difference
  double x1 = (ens - ensT).norm();
  ASSERT_EQ(x1, 0.);

  std::complex<double> c1 = complexZero;
  for (long ib1 = 0; ib1 < numBands; ib1++) {
    for (long ib2 = 0; ib2 < numBands; ib2++) {
      for (long ic = 0; ic < 3; ic++) {
        c1 += pow(velsT(ib1, ib2, ic) - vels(ib1, ib2, ic), 2);
      }
    }
  }
  ASSERT_EQ(c1, complexZero);

  std::complex<double> c2 = complexZero;
  for (long ib = 0; ib < numBands; ib++) {
    for (long ia = 0; ia < numAtoms; ia++) {
      for (long ic = 0; ic < 3; ic++) {
        c2 += pow(eigvecsT(ib*numAtoms+ia, ic) - eigvecs(ib, ia, ic), 2);
      }
    }
  }
  ASSERT_EQ(c2, complexZero);

  // now we check that we get the same eigenvectors in the two different
  // shapes
  auto k = point.getCoords(Points::cartesianCoords);
  auto [ensC, eigvecsC] = phononH0.diagonalizeFromCoords(k);
  x1 = (ens - ensC).norm();
  ASSERT_EQ(x1, 0.);

  Eigen::MatrixXcd eigvecsT2;
  s.getEigenvectors(eigvecsT2);
  c2 = complexZero;
  for (long ib1 = 0; ib1 < numBands; ib1++) {
    for (long ib2 = 0; ib2 < numBands; ib2++) {
      c2 += pow(eigvecsT2(ib1, ib2) - eigvecsC(ib1, ib2), 2);
    }
  }
  ASSERT_EQ(c2, complexZero);

  // we check what happens if we set eigenvectors as a matrix
  Eigen::MatrixXcd eigvecsC2(numBands, numBands);
  bandStructure.setEigenvectors(point, eigvecsC);
  bandStructure.getState(ik).getEigenvectors(eigvecsC2);
  c2 = complexZero;
  for (long ib1 = 0; ib1 < numBands; ib1++) {
    for (long ib2 = 0; ib2 < numBands; ib2++) {
      c2 += pow(eigvecsC2(ib1, ib2) - eigvecsC(ib1, ib2), 2);
    }
  }
  ASSERT_EQ(c2, complexZero);

  // make a comparison between tensor and matrix
  c2 = complexZero;
  for (long iband = 0; iband < numBands; iband++) {
    for (long iat = 0; iat < numAtoms; iat++) {
      for (long ipol = 0; ipol < 3; ipol++) {
        auto ind = compress2Indeces(iat, ipol, numAtoms, 3);
        c2 += pow(eigvecs(ipol, iat, iband) - eigvecsC(ind, iband), 2);
      }
    }
  }
  ASSERT_EQ(c2, complexZero);

  // now we use DetachedState, and verify that results are the same

  DetachedState d(k, ens, numAtoms, numBands, eigvecsC, &vels);

  // check energy
  auto ensD = d.getEnergies();
  x1 = (ens - ensD).norm();
  ASSERT_EQ(x1, 0.);

  // check velocity
  auto velsD = d.getVelocities();
  c1 = complexZero;
  for (long ib1 = 0; ib1 < numBands; ib1++) {
    for (long ib2 = 0; ib2 < numBands; ib2++) {
      for (long ic = 0; ic < 3; ic++) {
        c1 += pow(velsD(ib1, ib2, ic) - vels(ib1, ib2, ic), 2);
      }
    }
  }
  ASSERT_EQ(c1, complexZero);

  Eigen::MatrixXcd eigvecsD;
  d.getEigenvectors(eigvecsD);
  c2 = complexZero;
  for (long ib1 = 0; ib1 < numBands; ib1++) {
    for (long ib2 = 0; ib2 < numBands; ib2++) {
      c2 += pow(eigvecsD(ib1, ib2) - eigvecsC(ib1, ib2), 2);
    }
  }
  ASSERT_EQ(c2, complexZero);

  c2 = complexZero;
  Eigen::Tensor<std::complex<double>, 3> eigvecsDT(3, numAtoms, numBands);
  d.getEigenvectors(eigvecsDT);
  for (long ib = 0; ib < numBands; ib++) {
    for (long ia = 0; ia < numAtoms; ia++) {
      for (long ic = 0; ic < 3; ic++) {
        c2 += pow(eigvecsDT(ib, ia, ic) - eigvecs(ib, ia, ic), 2);
      }
    }
  }
  ASSERT_EQ(c2, complexZero);
}
