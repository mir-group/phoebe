#include "elph_qe_to_phoebe_app.h"
#include "qe_input_parser.h"
#include "gtest/gtest.h"

/** In this test we convert the input from QE to Phoebe data for Wannier
 * interpolation of el-ph coupling.
 * We check that if we Fourier back-transform to the same q/k point grid of QE
 * we find the same result.
 * By "same" result, note that only the square modulus of the el-ph coupling
 * |g|^2 (phases in the back-Fourier transform are still random)
 * summed over degenerate states, can be compared.
 * Note also that this test doesn't check that |g|^2 is smooth.
 */
TEST(InteractionElPh, Test1) {

  // setup input file
  Context context;
  context.setPhD2FileName("../example/development_silicon/silicon.fc");
  context.setElectronH0Name("../example/development_silicon/si_tb.dat");
  context.setWannier90Prefix("../example/development_silicon/si");
  context.setQuantumEspressoPrefix("../example/development_silicon/silicon");

  // actually, we only need the crystal
  auto t1 = QEParser::parsePhHarmonic(context);
  auto crystal = std::get<0>(t1);
  auto phononH0 = std::get<1>(t1);

  ElPhQeToPhoebeApp app;

  std::string phoebePrefixQE = context.getQuantumEspressoPrefix();
  auto t0 = app.readQEPhoebeHeader(crystal, phoebePrefixQE);
  Eigen::Vector3i qMesh = std::get<0>(t0);
  Eigen::Vector3i kMesh = std::get<1>(t0);
  Eigen::MatrixXd kgridFull = std::get<2>(t0);
  Eigen::MatrixXd qgridFull = std::get<3>(t0);
  Eigen::MatrixXd energies = std::get<4>(t0);
  int numIrrQPoints = std::get<5>(t0);
  int numQEBands = std::get<6>(t0);
  int numElectrons = std::get<7>(t0);
  int numSpin = std::get<8>(t0);

  FullPoints kPoints(crystal, kMesh);
  FullPoints qPoints(crystal, qMesh);

  int numModes = 3 * crystal.getNumAtoms();

  app.postProcessingWannier(
      context, crystal, phononH0, kPoints, qPoints, numQEBands, numModes,
      numIrrQPoints, numElectrons, numSpin, energies, kgridFull, kMesh, qMesh,
      true);
}
