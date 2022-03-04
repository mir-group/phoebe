#include "gtest/gtest.h"
#include "points.h"
#include "qe_input_parser.h"

/** Here I estimate the mass at the top of the valence band of silicon
 */
TEST (Kokkos, Wannier) {
  Context context;
  context.setElectronH0Name("../test/data/666_si_tb.dat");

  // setup crystal input
  Eigen::MatrixXd atomicPositions(2,3);
  atomicPositions.row(0) << 0., 0., 0.;
  atomicPositions.row(1) << 1.34940, 1.34940, 1.34940;
  Eigen::VectorXi atomicSpecies(2);
  atomicSpecies(0) = 0;
  atomicSpecies(1) = 0;
  std::vector<std::string> speciesNames;
  speciesNames.emplace_back("Si");
  context.setInputAtomicPositions(atomicPositions);
  context.setInputAtomicSpecies(atomicSpecies);
  context.setInputSpeciesNames(speciesNames);

  // read electron hamiltonian
  auto tup = QEParser::parseElHarmonicWannier(context);
  auto crystal = std::get<0>(tup);
  auto electronH0 = std::get<1>(tup);

  // gamma point
  Eigen::Vector3d k1 = {0.1,0.23,-0.17};

  // first, we diagonalize on the CPU
  auto tup1 = electronH0.diagonalizeFromCoordinates(k1);
  auto ens1 = std::get<0>(tup1);

  // now we add the code on Kokkos
  Eigen::VectorXd ens2(electronH0.getNumBands());
  {
    int numK = 1;

    // we need to copy the wavevectors to the GPU
    DoubleView2D q3Cs_d("q3", numK, 3);
    auto q3Cs_h = Kokkos::create_mirror_view(q3Cs_d);

#pragma omp parallel for
    for (int ik = 0; ik < numK; ik++) {
      for (int i = 0; i < 3; i++) {
        q3Cs_h(ik, i) = k1(i);
      }
    }
    Kokkos::deep_copy(q3Cs_d, q3Cs_h);

    auto t2 = electronH0.kokkosBatchedDiagonalizeFromCoordinates(q3Cs_d);
    DoubleView2D batchedEnergies = std::get<0>(t2);
    ComplexView3D batchedEigenvectors = std::get<1>(t2);

    // now we copy back to host
    auto tmpEnergies_h = Kokkos::create_mirror_view(batchedEnergies);
    Kokkos::deep_copy(tmpEnergies_h, batchedEnergies);
    for (int ib=0; ib<electronH0.getNumBands(); ++ib) {
      ens2(ib) = tmpEnergies_h(0, ib);
    }
  }

  // the two masses should be similar
  ASSERT_NEAR( (ens1-ens2).norm(), 0., 0.00001);
}
