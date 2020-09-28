#include "elph_bloch_to_wannier_app.h"
#include "bandstructure.h"
#include "eigen.h"
#include "interaction_elph.h"
#include "io.h"
#include "qe_input_parser.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

void ElPhQeToPhoebeApp::run(Context &context) {
  (void)context;

  if (mpi->mpiHead()) {
    std::cout << "\nLaunching app ElPhBlochToWannier\n" << std::endl;
  }

  std::string phoebePrefixQE = context.getQuantumEspressoPrefix();
  std::string wannierPrefix = context.getWannier90Prefix();
  std::string interpolation = context.getElPhInterpolation();

  // read Hamiltonian of phonons and electrons

  // actually, we only need the crystal
  auto t1 = QEParser::parsePhHarmonic(context);
  auto crystal = std::get<0>(t1);
  auto phononH0 = std::get<1>(t1);

  int numWannier; // number of bands used for g, at the end of the program
  if (interpolation == "wannier") {
    auto t2 = QEParser::parseElHarmonicWannier(context, &crystal);
    auto electronH0 = std::get<1>(t2);
    numWannier = electronH0.getNumBands();
  } else {
    auto t2 = QEParser::parseElHarmonicFourier(context);
    auto electronH0 = std::get<1>(t2);
    numWannier = electronH0.getNumBands();
  }

  auto t0 = readQEPhoebeHeader(crystal, phoebePrefixQE);
  Eigen::Vector3i qMesh = std::get<0>(t0);
  Eigen::Vector3i kMesh = std::get<1>(t0);
  Eigen::MatrixXd kgridFull = std::get<2>(t0);
  Eigen::MatrixXd qgridFull = std::get<3>(t0);
  Eigen::MatrixXd energies = std::get<4>(t0);
  int numIrrQPoints = std::get<5>(t0);
  int numQEBands = std::get<6>(t0);
  int numElectrons = std::get<7>(t0);
  int numSpin = std::get<8>(t0);

  //----------------------------------------------------------------------------
  // read Wannier90 rotation matrices

  FullPoints kPoints(crystal, kMesh);
  FullPoints qPoints(crystal, qMesh);

  int numModes = 3 * crystal.getNumAtoms();

  int numBands; // number of bands in the QE calculation
  Eigen::Tensor<std::complex<double>, 3> uMatrices;
  if (interpolation == "wannier") {
    uMatrices = setupRotationMatrices(wannierPrefix, kPoints);
    numBands = uMatrices.dimension(0);            // number of entangled bands
    assert(numWannier == uMatrices.dimension(1)); // number of entangled bands
  } else { // epa: here we preserve the number of bands
    numBands = numWannier;
  }

  //----------------------------------------------------------------------------

  // read coupling from file
  auto t5 =
      readGFromQEFile(context, numModes, numBands, numWannier, kPoints, qPoints,
                      kgridFull, numIrrQPoints, numQEBands, energies);
  auto gFull = std::get<0>(t5);          // (nBands, nBands, nModes, numK, numQ)
  auto phEigenvectors = std::get<1>(t5); // (numModes, numModes, numQPoints)
  auto phEnergies = std::get<2>(t5);     // (numModes, numQPoints)

  //----------------------------------------------------------------------------

  // Find the lattice vectors for the Fourier transforms

  auto t3 = crystal.buildWignerSeitzVectors(kMesh);
  Eigen::MatrixXd elBravaisVectors = std::get<0>(t3);
  Eigen::VectorXd elDegeneracies = std::get<1>(t3);

  auto t4 = crystal.buildWignerSeitzVectors(qMesh);
  Eigen::MatrixXd phBravaisVectors = std::get<0>(t4);
  Eigen::VectorXd phDegeneracies = std::get<1>(t4);

  //----------------------------------------------------------------------------

  // Bloch to Wannier transformation of el-ph coupling
  if (interpolation == "wannier") {

    if (mpi->mpiHead()) {
      std::cout << "Start Wannier-transform of g" << std::endl;
    }
    Eigen::Tensor<std::complex<double>, 5> g_wannier =
        blochToWannier(elBravaisVectors, phBravaisVectors, gFull, uMatrices,
                       phEigenvectors, kPoints, qPoints, crystal, phononH0);
    if (mpi->mpiHead()) {
      std::cout << "Done Wannier-transform of g\n" << std::endl;
    }

    //----------------------------------------------------------------------------

    // Dump el-ph in Wannier representation to file

    if (mpi->mpiHead()) {
      std::cout << "Start writing g to file" << std::endl;
      std::string outFileName = phoebePrefixQE + ".phoebe.elph.dat";
      std::ofstream outfile(outFileName);
      if (not outfile.is_open()) {
        Error e("Output file couldn't be opened");
      }
      outfile << numElectrons << " " << numSpin << "\n";
      outfile << kMesh << "\n";
      outfile << qMesh << "\n";
      outfile << phBravaisVectors.rows() << " " << phBravaisVectors.cols()
              << "\n";
      outfile << phBravaisVectors << "\n";
      outfile << phDegeneracies << "\n";
      outfile << elBravaisVectors.rows() << " " << elBravaisVectors.cols()
              << "\n";
      outfile << elBravaisVectors << "\n";
      outfile << elDegeneracies << "\n";
      outfile << "\n";
      for (auto x : g_wannier.dimensions()) {
        outfile << x << " ";
      }
      outfile << g_wannier << "\n";
      std::cout << "Done writing g to file\n" << std::endl;
    }
  } else {
    epaPostProcessing(context, gFull, energies, phEnergies, kPoints, qPoints,
                      numElectrons, numSpin);
  }
}

void ElPhQeToPhoebeApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");
  throwErrorIfUnset(context.getPhD2FileName(), "PhD2FileName");
  throwErrorIfUnset(context.getQuantumEspressoPrefix(),
                    "QuantumEspressoPrefix");

  choices = {"wannier", "epa"};
  std::string x = context.getElPhInterpolation();
  throwErrorIfUnset(x, "elPhInterpolation");
  if (std::find(choices.begin(), choices.end(), x) == choices.end()) {
    Error e("The elPhInterpolation value has not been recognized.");
  }

  if (x == "wannier") {
    throwErrorIfUnset(context.getWannier90Prefix(), "Wannier90Prefix");
  } else {
    throwErrorIfUnset(context.getEpaDeltaEnergy(), "epaDeltaEnergy");
    throwErrorIfUnset(context.getEpaSmearingEnergy(), "epaSmearingEnergy");
  }
}

Eigen::Tensor<std::complex<double>, 5> ElPhQeToPhoebeApp::blochToWannier(
    const Eigen::MatrixXd &elBravaisVectors,
    const Eigen::MatrixXd &phBravaisVectors,
    Eigen::Tensor<std::complex<double>, 5> &gFull,
    const Eigen::Tensor<std::complex<double>, 3> &uMatrices,
    const Eigen::Tensor<std::complex<double>, 3> &phEigenvectors,
    FullPoints &kPoints, FullPoints &qPoints, Crystal &crystal,
    PhononH0 &phononH0) {

  int numBands = gFull.dimension(0); // # of entangled bands
  int numModes = gFull.dimension(2);
  int numKPoints = gFull.dimension(3);
  int numQPoints = gFull.dimension(4);
  int numElBravaisVectors = elBravaisVectors.cols();
  int numPhBravaisVectors = phBravaisVectors.cols();
  int numWannier = uMatrices.dimension(1);

  std::array<Eigen::Index, 5> zeros;
  for (auto &s : zeros) {
    s = 0;
  }

  bool usePolarCorrection;
  Eigen::Matrix3d epsilon = phononH0.getDielectricMatrix();
  if (epsilon.squaredNorm() > 1.0e-10) { // i.e. if epsilon wasn't computed
    if (crystal.getNumSpecies() > 1) {   // otherwise polar correction = 0
      usePolarCorrection = true;
    }
  }
  if (usePolarCorrection) {
    // we need to subtract the polar correction
    // this contribution will be reinstated during the interpolation
    auto volume = crystal.getVolumeUnitCell();
    auto reciprocalUnitCell = crystal.getReciprocalUnitCell();
    auto epsilon = phononH0.getDielectricMatrix();
    auto bornCharges = phononH0.getBornCharges();
    auto atomicPositions = crystal.getAtomicPositions();
    auto qCoarseMesh = phononH0.getCoarseGrid();

    for (long iq = 0; iq < numQPoints; iq++) {
      Eigen::Vector3d q = qPoints.getPointCoords(iq, Points::cartesianCoords);
      if (q.norm() > 1.0e-8) {

        Eigen::MatrixXcd ev3(numModes, numModes);
        for (int j = 0; j < numModes; j++) {
          for (int i = 0; i < numModes; i++) {
            ev3(i, j) = phEigenvectors(i, j, iq);
          }
        }

        for (long ik = 0; ik < numKPoints; ik++) {
          Eigen::Vector3d k =
              kPoints.getPointCoords(ik, Points::cartesianCoords);

          // Coordinates and index of k+q point
          Eigen::Vector3d kq = k + q;
          Eigen::Vector3d kqCrystal = kPoints.cartesianToCrystal(kq);
          long ikq = kPoints.getIndex(kqCrystal);

          // gather eigenvectors
          Eigen::MatrixXcd ev1(numBands, numWannier);
          Eigen::MatrixXcd ev2(numBands, numWannier);
          for (int j = 0; j < numBands; j++) {
            for (int i = 0; i < numBands; i++) {
              ev1(i, j) = uMatrices(i, j, ik);
              ev2(i, j) = uMatrices(i, j, ikq);
            }
          }
          ev1 = ev1.adjoint();
          ev2 = ev2.adjoint();

          auto v = InteractionElPhWan::getPolarCorrectionStatic(
              q, ev1, ev2, ev3, volume, reciprocalUnitCell, epsilon,
              bornCharges, atomicPositions, qCoarseMesh);
          for (int nu = 0; nu < numModes; nu++) {
            for (int j = 0; j < numBands; j++) {
              for (int i = 0; i < numBands; i++) {
                gFull(i, j, nu, ik, iq) -= v(i, j, nu);
              }
            }
          }
        }
      }
    }
  }

  if (mpi->mpiHead()) {
    std::cout << "Wannier rotation" << std::endl;
  }
  Eigen::Tensor<std::complex<double>, 5> gFullTmp(
      numWannier, numWannier, numModes, numKPoints, numQPoints);
  gFullTmp.setZero();
  for (long iq = 0; iq < numQPoints; iq++) {
    Eigen::Vector3d q = qPoints.getPointCoords(iq, Points::cartesianCoords);
    for (long ik = 0; ik < numKPoints; ik++) {
      Eigen::Vector3d k = kPoints.getPointCoords(ik, Points::cartesianCoords);

      // Coordinates and index of k+q point
      Eigen::Vector3d kq = k + q;
      Eigen::Vector3d kqCrystal = kPoints.cartesianToCrystal(kq);
      long ikq = kPoints.getIndex(kqCrystal);

      // First we transform from the Bloch to Wannier Gauge

      Eigen::MatrixXcd uK(numWannier, numBands);
      Eigen::MatrixXcd uKq(numWannier, numBands);
      for (int i = 0; i < numBands; i++) {
        for (int j = 0; j < numWannier; j++) {
          uK(i, j) = uMatrices(i, j, ik);
          uKq(i, j) = uMatrices(i, j, ikq);
        }
      }
      uKq = uKq.adjoint();

      Eigen::Tensor<std::complex<double>, 3> tmp(numWannier, numBands,
                                                 numModes);
      tmp.setZero();
      for (int nu = 0; nu < numModes; nu++) {
        for (int i = 0; i < numWannier; i++) {
          for (int j = 0; j < numBands; j++) {
            for (int l = 0; l < numBands; l++) {
              tmp(i, j, nu) += uKq(i, l) * gFull(l, j, nu, ik, iq);
            }
          }
        }
      }
      for (int nu = 0; nu < numModes; nu++) {
        for (int i = 0; i < numWannier; i++) {
          for (int j = 0; j < numWannier; j++) {
            for (int l = 0; l < numBands; l++) {
              gFullTmp(i, j, nu, ik, iq) += tmp(i, l, nu) * uK(l, j);
            }
          }
        }
      }
    } // ik
  }   // iq
  gFull.reshape(zeros);

  if (mpi->mpiHead()) {
    std::cout << "Electronic Fourier Transform" << std::endl;
  }
  // Fourier transform on the electronic coordinates
  Eigen::Tensor<std::complex<double>, 5> gMixed(
      numWannier, numWannier, numModes, numElBravaisVectors, numQPoints);
  gMixed.setZero();
  for (int iR = 0; iR < numElBravaisVectors; iR++) {
    for (long ik = 0; ik < numKPoints; ik++) {
      Eigen::Vector3d k = kPoints.getPointCoords(ik, Points::cartesianCoords);
      double arg = k.dot(elBravaisVectors.col(iR));
      std::complex<double> phase = exp(-complexI * arg) / double(numKPoints);
      for (int iq = 0; iq < numQPoints; iq++) {
        for (int i = 0; i < numWannier; i++) {
          for (int j = 0; j < numWannier; j++) {
            for (int nu = 0; nu < numModes; nu++) {
              gMixed(i, j, nu, iR, iq) += gFullTmp(i, j, nu, ik, iq) * phase;
            }
          }
        }
      }
    }
  } // iq
  gFullTmp.reshape(zeros);

  if (mpi->mpiHead()) {
    std::cout << "Phonon rotation" << std::endl;
  }
  Eigen::Tensor<std::complex<double>, 5> gWannierTmp(
      numWannier, numWannier, numModes, numElBravaisVectors, numQPoints);
  gWannierTmp.setZero();
  for (long iq = 0; iq < numQPoints; iq++) {

    Eigen::MatrixXcd uQ(numModes, numModes);
    for (int nu = 0; nu < numModes; nu++) {
      for (int nu2 = 0; nu2 < numModes; nu2++) {
        uQ(nu, nu2) = phEigenvectors(nu, nu2, iq);
      }
    }
    uQ = uQ.inverse();

    for (int nu = 0; nu < numModes; nu++) {
      for (int nu2 = 0; nu2 < numModes; nu2++) {
        for (int irE = 0; irE < numElBravaisVectors; irE++) {
          for (int i = 0; i < numWannier; i++) {
            for (int j = 0; j < numWannier; j++) {
              gWannierTmp(i, j, nu, irE, iq) +=
                  gMixed(i, j, nu2, irE, iq) * uQ(nu2, nu);
            }
          }
        }
      }
    }
  }
  gMixed.reshape(zeros);

  if (mpi->mpiHead()) {
    std::cout << "Phonon Fourier Transform" << std::endl;
  }
  Eigen::Tensor<std::complex<double>, 5> gWannier(numWannier, numWannier,
                                                  numModes, numElBravaisVectors,
                                                  numPhBravaisVectors);
  gWannier.setZero();
  for (int iq = 0; iq < numQPoints; iq++) {
    Eigen::Vector3d q = qPoints.getPointCoords(iq, Points::cartesianCoords);
    for (int irP = 0; irP < numPhBravaisVectors; irP++) {
      double arg = q.dot(phBravaisVectors.col(irP));
      std::complex<double> phase = exp(-complexI * arg) / double(numQPoints);
      for (int irE = 0; irE < numElBravaisVectors; irE++) {
        for (int i = 0; i < numWannier; i++) {
          for (int j = 0; j < numWannier; j++) {
            for (int nu = 0; nu < numModes; nu++) {
              gWannier(i, j, nu, irE, irP) +=
                  phase * gWannierTmp(i, j, nu, irE, iq);
            }
          }
        }
      }
    }
  }
  gWannierTmp.reshape(zeros);
  return gWannier;
}

Eigen::Tensor<std::complex<double>, 3>
ElPhQeToPhoebeApp::setupRotationMatrices(const std::string &wannierPrefix,
                                         FullPoints &fullPoints) {
  std::string line;

  if (wannierPrefix.empty()) {
    Error e("Must provide an input H0 file name");
  }

  std::string fileName = wannierPrefix + "_u.mat";

  // open input file
  std::ifstream infile(fileName);
  if (not infile.is_open()) {
    Error e("U-matrix file not found");
  }

  // Title line
  std::getline(infile, line);

  int numPoints, numWannier, tmpI;
  infile >> numPoints >> numWannier >> tmpI;

  assert(numPoints == fullPoints.getNumPoints());

  Eigen::Tensor<std::complex<double>, 3> uMatrix(numWannier, numWannier,
                                                 numPoints);
  uMatrix.setZero();

  for (int ik = 0; ik < numPoints; ik++) {
    // empty line
    std::getline(infile, line);

    Eigen::Vector3d thisK; // vector in crystal coords
    infile >> thisK(0) >> thisK(1) >> thisK(2);

    long ikk = fullPoints.getIndex(thisK);

    double re, im;
    for (int j = 0; j < numWannier; j++) {
      for (int i = 0; i < numWannier; i++) {
        infile >> re >> im;
        uMatrix(i, j, ikk) = {re, im};
      }
    }
  }
  infile.close();

  // ---------------------------------------------------------------------

  // Now we get the disentanglement matrix

  std::string fileNameDis = wannierPrefix + "_u_dis.mat";

  // open input file
  std::ifstream infileDis(fileNameDis);
  if (not infileDis.is_open()) {
    // if the disentanglement file is not found
    // we assume there's no disentanglement and quit the function
    return uMatrix;
  } // else, we parse the file

  // Title line
  std::getline(infileDis, line);

  int numPoints2, numWannier2, numBands;
  infileDis >> numPoints >> numWannier2 >> numBands;

  assert(numPoints2 == numPoints);
  (void)numPoints2;
  assert(numWannier2 == numWannier);
  assert(numBands >= numWannier);

  Eigen::Tensor<std::complex<double>, 3> uMatrixDis(numBands, numWannier,
                                                    numPoints);
  uMatrixDis.setZero();

  for (int ik = 0; ik < numPoints; ik++) {
    // empty line
    std::getline(infileDis, line);

    double x, y, z;
    infileDis >> x >> y >> z;
    Eigen::Vector3d thisK;
    thisK << x, y, z; // vector in crystal coords

    long ikk = fullPoints.getIndex(thisK);

    double re, im;
    for (int j = 0; j < numWannier; j++) {
      for (int i = 0; i < numBands; i++) {
        infileDis >> re >> im;
        uMatrixDis(i, j, ikk) = {re, im};
      }
    }
  }
  infileDis.close();

  // Now I multiply the two rotation matrices

  Eigen::Tensor<std::complex<double>, 3> u(numBands, numWannier, numPoints);
  u.setZero();
  for (int ivec = 0; ivec < numPoints; ivec++) {

    Eigen::MatrixXcd a(numBands, numWannier);
    for (int i = 0; i < numBands; i++) {
      for (int j = 0; j < numWannier; j++) {
        a(i, j) = uMatrixDis(i, j, ivec);
      }
    }

    Eigen::MatrixXcd b(numWannier, numWannier);
    for (int i = 0; i < numWannier; i++) {
      for (int j = 0; j < numWannier; j++) {
        b(i, j) = uMatrix(i, j, ivec);
      }
    }

    Eigen::MatrixXcd c = a * b; // matrix of shape (numBands, numWannier)

    for (int i = 0; i < numBands; i++) {
      for (int j = 0; j < numWannier; j++) {
        u(i, j, ivec) = c(i, j);
      }
    }
  }

  return u;
}

int ElPhQeToPhoebeApp::computeOffset(const Eigen::MatrixXd &energies,
                                     const std::string &wannierPrefix) {
  Eigen::VectorXd energiesQEAtZero = energies.col(0); // k = 0

  { // check the first point in Wannier90 is gamma
    std::string fileName = wannierPrefix + ".nnkp";
    std::ifstream infile(fileName);
    std::string line;
    for (int i = 0; i < 18; i++) {
      std::getline(infile, line); // skip the first 18 lines
    }
    double kx, ky, kz;
    infile >> kx >> ky >> kz;
    if (kx * kx + ky * ky + kz * kz > 1.0e-5) {
      Error e("Expecting first coarse grid kpoint in Wannier90 to be gamma");
    }
  }

  // read .eig file to get energies

  std::vector<double> energiesWannierAtZero;
  {
    std::string eigFileName = wannierPrefix + ".eig";
    std::ifstream eigfile(eigFileName);
    int ib, ik;
    double x;
    while (eigfile >> ib >> ik >> x) {
      if (ik > 1) {
        break;
      }
      x /= energyRyToEv;
      energiesWannierAtZero.push_back(x);
    }
  }

  int numBandsWannier = energiesWannierAtZero.size();
  int numFull = energiesQEAtZero.size();

  // we find the offset by comparing the energy differences
  // the offset which minimizes energy differences is the chosen one
  int possibleValues = numFull - numBandsWannier + 1;
  Eigen::VectorXd difference(possibleValues);
  difference.setZero();
  for (int i = 0; i < possibleValues; i++) {
    for (int ib = 0; ib < numBandsWannier; ib++) {
      difference(i) =
          pow(energiesQEAtZero(ib) - energiesWannierAtZero[ib + i], 2);
    }
  }

  // offset = index of min difference
  int offset = -1;
  for (int i = 0; i < possibleValues; i++) {
    if (difference(i) == difference.minCoeff()) {
      offset = i;
      break;
    }
  }

  if (offset == -1) {
    Error e("Bands offset not found");
  }

  return offset;
}

// read g, which is written to file on all k, q points
std::tuple<Eigen::Tensor<std::complex<double>, 5>,
           Eigen::Tensor<std::complex<double>, 3>, Eigen::MatrixXd>
ElPhQeToPhoebeApp::readGFromQEFile(Context &context, const int &numModes,
                                   const int &numBands, const int &numWannier,
                                   FullPoints &kPoints, FullPoints &qPoints,
                                   const Eigen::MatrixXd &kgridFull,
                                   const int &numIrrQPoints,
                                   const int &numQEBands,
                                   const Eigen::MatrixXd &energies) {

  if (mpi->mpiHead()) {
    std::cout << "Start reading from file" << std::endl;
  }

  std::string interpolation = context.getElPhInterpolation();
  int bandsOffset;
  if (interpolation == "wannier") {
    std::string wannierPrefix = context.getWannier90Prefix();
    bandsOffset = computeOffset(energies, wannierPrefix);
  } else {
    bandsOffset = 0;
  }

  long numKPoints = kPoints.getNumPoints();
  long numQPoints = qPoints.getNumPoints();
  Eigen::Tensor<std::complex<double>, 5> g_full(numBands, numBands, numModes,
                                                numKPoints, numQPoints);
  Eigen::Tensor<std::complex<double>, 3> phEigenvectors(numModes, numModes,
                                                        numQPoints);
  Eigen::MatrixXd phEnergies(numModes, numQPoints);

  g_full.setZero();
  phEigenvectors.setZero();
  phEnergies.setZero();

  Eigen::VectorXi ikMap(numKPoints);
  for (long ikOld = 0; ikOld < numKPoints; ikOld++) {
    Eigen::Vector3d kOld = kgridFull.col(ikOld);
    long ikNew = kPoints.getIndex(kOld);
    ikMap(ikOld) = ikNew;
  }

  std::string phoebePrefixQE = context.getQuantumEspressoPrefix();

  for (int iqIrr = 0; iqIrr < numIrrQPoints; iqIrr++) {
    std::stringstream ss;
    ss << std::setw(4) << std::setfill('0') << iqIrr + 1;
    std::string numString = ss.str();
    std::string fileName = phoebePrefixQE + ".phoebe." + numString + ".dat";
    std::ifstream infileQ(fileName);

    std::string line;

    int nqStar; // number of reducible q points in this file
    infileQ >> nqStar;
    std::vector<Eigen::Vector3d> qStar;
    for (int iq = 0; iq < nqStar; iq++) {
      Eigen::Vector3d thisQ; // in crystal coordinates
      infileQ >> thisQ(0) >> thisQ(1) >> thisQ(2);
      qStar.push_back(thisQ);
    }

    Eigen::VectorXd phononEnergies(numModes);
    for (int nu = 0; nu < numModes; nu++) {
      infileQ >> phononEnergies(nu);
    }

    Eigen::Tensor<std::complex<double>,3> phononEigenvectorsStar(numModes, numModes, nqStar);
    for (int iq = 0; iq < nqStar; iq++) {
      for (int j = 0; j < numModes; j++) {
        for (int i = 0; i < numModes; i++) {
          double re, im;
          infileQ >> re >> im;
          phononEigenvectorsStar(i, j, iq) = {re, im}; // j is the eig index
        }
      }
    }
    std::getline(infileQ, line); // empty line

    // read the g-coupling
    Eigen::Tensor<std::complex<double>, 5> thisG(numQEBands, numQEBands,
                                                 numModes, numKPoints, nqStar);
    thisG.setZero();
    for (int iq = 0; iq < nqStar; iq++) {
      for (int nu = 0; nu < numModes; nu++) {
        for (int ik = 0; ik < numKPoints; ik++) {
          for (int ib2 = 0; ib2 < numQEBands; ib2++) {
            for (int ib1 = 0; ib1 < numQEBands; ib1++) {
              double re, im;
              infileQ >> re >> im;
              thisG(ib1, ib2, nu, ik, iq) = {re, im};
            }
          }
        }
      }
    }
    infileQ.close();

    // reorder the q/k indices
    for (int iqStar = 0; iqStar < nqStar; iqStar++) {
      Eigen::Vector3d qVec = qStar[iqStar];
      long iqFull = qPoints.getIndex(qVec);

      for (int nu = 0; nu < numModes; nu++) {
        for (int ik = 0; ik < numKPoints; ik++) {
          for (int ib2 = 0; ib2 < numWannier; ib2++) {
            for (int ib1 = 0; ib1 < numWannier; ib1++) {
              g_full(ib1, ib2, nu, ikMap(ik), iqFull) =
                  thisG(bandsOffset + ib1, bandsOffset + ib2, nu, ik, iqStar);
            }
          }
        }
      }

      for (int j = 0; j < numModes; j++) {
        for (int i = 0; i < numModes; i++) {
          phEigenvectors(i, j, iqFull) = phononEigenvectorsStar(i, j, iqStar);
        }
      }

      for (int i = 0; i < numModes; i++) {
        phEnergies(i, iqFull) = phononEnergies(i);
      }
    }
  }
  if (mpi->mpiHead()) {
    std::cout << "Done reading g from file" << std::endl << std::endl;
  }
  return {g_full, phEigenvectors, phEnergies};
}

std::tuple<Eigen::Vector3i, Eigen::Vector3i, Eigen::MatrixXd, Eigen::MatrixXd,
           Eigen::MatrixXd, int, int, int, int>
ElPhQeToPhoebeApp::readQEPhoebeHeader(Crystal &crystal,
                                      const std::string &phoebePrefixQE) {
  std::string fileName = phoebePrefixQE + ".phoebe.0000.dat";
  std::ifstream infile(fileName);
  std::string line;
  if (not infile.is_open()) {
    Error e("H0 file not found");
  }
  std::getline(infile, line); // first line is a title
  int numQEBands;             // number of Kohn-Sham states
  double numElectrons;        // number of electrons (spin degeneracy included)
  int numSpin;                // should always be one, without support for spin
  infile >> numQEBands >> numElectrons >> numSpin;
  Eigen::Vector3i kMesh, qMesh;
  infile >> qMesh(0) >> qMesh(1) >> qMesh(2) >> kMesh(0) >> kMesh(1) >>
      kMesh(2);
  int bogusI;
  double bogusD;
  int numAtoms;
  infile >> bogusD >> numAtoms; // alat and nat
  assert(numAtoms == crystal.getNumAtoms());
  (void)crystal;
  // unit cell
  for (int i = 0; i < 9; i++) {
    infile >> bogusD;
  }
  // reciprocal unit cell
  for (int i = 0; i < 9; i++) {
    infile >> bogusD;
  }
  // ityp
  for (int i = 0; i < numAtoms; i++) {
    infile >> bogusI;
  }
  // positions
  for (int i = 0; i < 3 * numAtoms; i++) {
    infile >> bogusD;
  }

  int numQPoints, numIrrQPoints;
  infile >> numQPoints >> numIrrQPoints;

  Eigen::MatrixXd qgridFull(3, numQPoints);
  for (int iq = 0; iq < numQPoints; iq++) {
    infile >> qgridFull(0, iq) >> qgridFull(1, iq) >> qgridFull(2, iq);
  }

  int numKPoints;
  infile >> numKPoints;
  Eigen::MatrixXd kgridFull(3, numKPoints);
  for (int ik = 0; ik < numKPoints; ik++) {
    infile >> kgridFull(0, ik) >> kgridFull(1, ik) >> kgridFull(2, ik);
  }

  Eigen::MatrixXd energies(numQEBands, numKPoints);
  for (int ik = 0; ik < numKPoints; ik++) {
    for (int ib = 0; ib < numQEBands; ib++) {
      infile >> energies(ib, ik);
    }
  }

  return {qMesh,         kMesh,      kgridFull,    qgridFull, energies,
          numIrrQPoints, numQEBands, numElectrons, numSpin};
}

void ElPhQeToPhoebeApp::epaPostProcessing(
    Context &context, Eigen::Tensor<std::complex<double>, 5> gFull,
    Eigen::MatrixXd &elEnergies, Eigen::MatrixXd &phEnergies,
    FullPoints &kPoints, FullPoints &qPoints, const int &numElectrons,
    const int &numSpin) {
  // input
  double smearing = context.getEpaDeltaEnergy();
  double smearing2 = 2. * smearing * smearing;
  double deltaEnergy = context.getEpaDeltaEnergy();

  // prepare energy bins
  double minEnergy = elEnergies.minCoeff();
  double maxEnergy = elEnergies.maxCoeff();
  long numEpaEnergies = (maxEnergy - minEnergy) / deltaEnergy + 1;
  Eigen::VectorXd epaEnergies(numEpaEnergies);
  for (int i = 0; i < numEpaEnergies; i++) {
    epaEnergies[i] = i * deltaEnergy + minEnergy;
  }

  if (mpi->mpiHead()) {
    std::cout << "Building EPA with " << numEpaEnergies << " energy bins.";
  }

  int numModes = gFull.dimension(2);
  int numBands = gFull.dimension(0);

  Eigen::Tensor<double, 3> g2Epa(numModes, numEpaEnergies, numEpaEnergies);
  g2Epa.setZero();

  int numKPoints = gFull.dimension(3);
  int numQPoints = gFull.dimension(4);

  Eigen::Tensor<double, 3> gaussians(numEpaEnergies, numBands, numKPoints);
  for (int ib1 = 0; ib1 < numBands; ib1++) {
    for (int ik = 0; ik < numKPoints; ik++) {
      for (int i = 0; i < numEpaEnergies; i++) {
        double arg = pow(elEnergies(ib1, ik) - epaEnergies(i), 2) / smearing2;
        gaussians(i, ib1, ik) = exp(-arg);
      }
    }
  }

  LoopPrint loopPrint("Computing coupling EPA", "q-points", numQPoints);
  for (int iq = 0; iq < numQPoints; iq++) {
    loopPrint.update();
    Eigen::Vector3d q = qPoints.getPointCoords(iq, Points::cartesianCoords);
    for (int ik = 0; ik < numKPoints; ik++) {
      Eigen::Vector3d k = kPoints.getPointCoords(ik, Points::cartesianCoords);

      // Coordinates and index of k+q point
      Eigen::Vector3d kq = k + q;
      Eigen::Vector3d kqCrystal = kPoints.cartesianToCrystal(kq);
      int ikq = kPoints.getIndex(kqCrystal);

      for (int ib1 = 0; ib1 < numBands; ib1++) {
        for (int ib2 = 0; ib2 < numBands; ib2++) {

          for (int j = 0; j < numEpaEnergies; j++) {
            for (int i = 0; i < numEpaEnergies; i++) {

              double gaussian = gaussians(i, ib1, ik) * gaussians(j, ib2, ikq);

              for (int nu = 0; nu < numModes; nu++) {
                g2Epa(nu, i, j) += std::norm(gFull(ib1, ib2, nu, ik, iq)) *
                                   gaussian / 2. / phEnergies(nu, iq);
                // /2omega, because there is a difference between the
                // coupling <k+q| dV_q |k> from quantum espresso
                // and the coupling g to be used for transport calcs
              }
            }
          }
        }
      }
    }
  }
  loopPrint.close();

  Eigen::VectorXd phAvgEnergies(numModes); // phEnergies(numModes, numQPoints);
  for (int nu = 0; nu < numModes; nu++) {
    phAvgEnergies(nu) = phEnergies.row(nu).sum() / phEnergies.cols();
  }

  if (mpi->mpiHead()) {
    std::cout << "\nStart writing g to file" << std::endl;
    std::string phoebePrefixQE = context.getQuantumEspressoPrefix();
    std::string outFileName = phoebePrefixQE + ".phoebe.epa.dat";
    std::ofstream outfile(outFileName);
    if (not outfile.is_open()) {
      Error e("Output file couldn't be opened");
    }
    outfile << numElectrons << " " << numSpin << "\n";
    outfile << numEpaEnergies << "\n";
    outfile << epaEnergies.transpose() << "\n";
    outfile << phAvgEnergies.size() << "\n";
    outfile << phAvgEnergies.transpose() << "\n";
    outfile << "\n";
    for (auto x : g2Epa.dimensions()) {
      outfile << x << " ";
    }
    outfile << g2Epa << "\n";
    std::cout << "Done writing g to file\n" << std::endl;
  }
}
