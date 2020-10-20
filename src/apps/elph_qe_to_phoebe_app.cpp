#include "elph_qe_to_phoebe_app.h"
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

  // actually, we only need the crystal
  auto t1 = QEParser::parsePhHarmonic(context);
  auto crystal = std::get<0>(t1);
  auto phononH0 = std::get<1>(t1);

  std::string phoebePrefixQE = context.getQuantumEspressoPrefix();
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

  FullPoints kPoints(crystal, kMesh);
  FullPoints qPoints(crystal, qMesh);

  int numModes = 3 * crystal.getNumAtoms();

  std::string interpolation = context.getElPhInterpolation();

  if (interpolation == "wannier") {

    postProcessingWannier(context, crystal, phononH0, kPoints, qPoints,
                          numQEBands, numModes, numIrrQPoints, numElectrons,
                          numSpin, energies, kgridFull, kMesh, qMesh);

  } else { // EPA

    epaPostProcessing(context, energies, kPoints, qPoints, numElectrons,
                      numSpin, numModes, numIrrQPoints, numQEBands, energies,
                      kgridFull);
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

  if (mpi->mpiHead()) {
    std::cout << "Start Wannier-transform of g" << std::endl;
  }

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

  bool usePolarCorrection = false;
  Eigen::Matrix3d epsilon = phononH0.getDielectricMatrix();
  if (epsilon.squaredNorm() > 1.0e-10) { // i.e. if epsilon wasn't computed
    if (crystal.getNumSpecies() > 1) {   // otherwise polar correction = 0
      usePolarCorrection = true;
    }
  }

  if (usePolarCorrection) {
    std::cout << "Polar correction\n";
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
          // ev1 = ev1.adjoint();
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

      // u has size (numBands, numWannier, numKPoints)
      Eigen::MatrixXcd uK(numBands, numWannier);
      Eigen::MatrixXcd uKq(numBands, numWannier);
      for (int i = 0; i < numBands; i++) {
        for (int j = 0; j < numWannier; j++) {
          uK(i, j) = uMatrices(i, j, ik);
          uKq(i, j) = uMatrices(i, j, ikq);
        }
      }
      Eigen::MatrixXcd uKDagger(numWannier, numBands);
      uKDagger = uK.adjoint();

      Eigen::Tensor<std::complex<double>, 3> tmp(numWannier, numBands,
                                                 numModes);
      tmp.setZero();
      for (int nu = 0; nu < numModes; nu++) {
        for (int i = 0; i < numWannier; i++) {
          for (int j = 0; j < numBands; j++) {
            for (int l = 0; l < numBands; l++) {
              // ukq has size(numWannier,numBands)
              // gFull has size numBands,numBands,...
              tmp(i, j, nu) += uKq(l, i) * gFull(l, j, nu, ik, iq);
            }
          }
        }
      }
      for (int nu = 0; nu < numModes; nu++) {
        for (int i = 0; i < numWannier; i++) {
          for (int j = 0; j < numWannier; j++) {
            for (int l = 0; l < numBands; l++) {
              gFullTmp(i, j, nu, ik, iq) += tmp(i, l, nu) * uKDagger(j, l);
            }
          }
        }
      }
    } // ik
  }   // iq

  //  gFull.reshape(zeros);

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
    auto uQM1 = uQ.inverse();
    // this isn't equal to the adjoint, due to mass renormalization

    for (int nu = 0; nu < numModes; nu++) {
      for (int nu2 = 0; nu2 < numModes; nu2++) {
        for (int irE = 0; irE < numElBravaisVectors; irE++) {
          for (int i = 0; i < numWannier; i++) {
            for (int j = 0; j < numWannier; j++) {
              gWannierTmp(i, j, nu, irE, iq) +=
                  gMixed(i, j, nu2, irE, iq) * uQM1(nu2, nu);
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
                                                  numModes, numPhBravaisVectors,
                                                  numElBravaisVectors);
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
              gWannier(j, i, nu, irP, irE) +=
                  phase * gWannierTmp(i, j, nu, irE, iq);
            }
          }
        }
      }
    }
  }
  gWannierTmp.reshape(zeros);

  if (mpi->mpiHead()) {
    std::cout << "Done Wannier-transform of g\n" << std::endl;
  }

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

    for (int i = 0; i < numBands; i++) {
      for (int j = 0; j < numWannier; j++) {
        for (int k = 0; k < numWannier; k++) {
          u(i, j, ivec) += b(k, j) * a(i, k);
        }
      }
    }
  } // ivec
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
    std::cout << "Start reading el-ph coupling from file" << std::endl;
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

    for (int iq = 0; iq < nqStar; iq++) {
      Eigen::Vector3d thisQ; // in same as above, in cartesian coordinates
      infileQ >> thisQ(0) >> thisQ(1) >> thisQ(2);
    }

    Eigen::VectorXd phononEnergies(numModes);
    for (int nu = 0; nu < numModes; nu++) {
      infileQ >> phononEnergies(nu);
    }

    Eigen::Tensor<std::complex<double>, 3> phononEigenvectorsStar(
        numModes, numModes, nqStar);
    for (int iq = 0; iq < nqStar; iq++) {
      for (int j = 0; j < numModes; j++) {
        for (int i = 0; i < numModes; i++) {
          // Note, in Fortran I was writing:
          // do jj = 1,nmodes
          //   do k = 1,nat
          //     do i_cart = 1,3
          // This has to be aligned with what done by PhononH0
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
    std::cout << "Done reading el-ph coupling from file\n" << std::endl;
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


void ElPhQeToPhoebeApp::epaPostProcessing(Context &context, Eigen::MatrixXd &elEnergies,
                       FullPoints &kPoints, FullPoints &qPoints,
                       const int &numElectrons, const int &numSpin,
                       const int &numModes, const int &numIrrQPoints,
                       const int &numQEBands, const Eigen::MatrixXd &energies,
                       const Eigen::MatrixXd &kgridFull) {

  if (mpi->mpiHead()) {
    std::cout << "Starting EPA post-processing\n" << std::endl;
  }

  auto t2 = QEParser::parseElHarmonicFourier(context);
  auto electronH0 = std::get<1>(t2);
  int numBands = electronH0.getNumBands();

  // read coupling from file
  auto t5 =
      readGFromQEFile(context, numModes, numBands, numBands, kPoints, qPoints,
                      kgridFull, numIrrQPoints, numQEBands, energies);
  auto gFull = std::get<0>(t5);          // (nBands, nBands, nModes, numK, numQ)
  auto phEigenvectors = std::get<1>(t5); // (numModes, numModes, numQPoints)
  auto phEnergies = std::get<2>(t5);     // (numModes, numQPoints)

  assert(numBands == gFull.dimension(0));
  assert(numModes == gFull.dimension(2));

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

  Eigen::Tensor<double, 3> g2Epa(numModes, numEpaEnergies, numEpaEnergies);
  g2Epa.setZero();

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

void ElPhQeToPhoebeApp::testElectronicTransform(
    Points &kPoints, const std::string &wannierPrefix,
    const Eigen::MatrixXd &elBravaisVectors,
    const Eigen::Tensor<std::complex<double>, 3> &uMatrices,
    const Eigen::VectorXd &elDegeneracies, ElectronH0Wannier &electronH0) {
  /** This is a simple test:
   * 1) Fourier Transform the electronic Hamiltonian to Wannier representation
   *    Here I use the U matrices from file
   * 2) FT back to Bloch representation, using the U matrices from ElectronH0
   *    on the original grid of kpoints
   * If everything works, I expect to find the same electronic energies
   * Phases of rotation matrices in the back-FT will be random.
   */

  int numBands = uMatrices.dimension(0);
  int numWannier = uMatrices.dimension(1);
  assert(numBands >= numWannier);

  Eigen::MatrixXd blochEnergies(numBands, kPoints.getNumPoints());
  blochEnergies.setZero();

  auto t = kPoints.getMesh();
  auto kMesh = std::get<0>(t);

  // I try the FFT of the energies
  for (int ik = 0; ik < kPoints.getNumPoints(); ik++) {
    auto kCryst = kPoints.getPointCoords(ik);
    kCryst(0) *= kMesh(0);
    kCryst(1) *= kMesh(1);
    kCryst(2) *= kMesh(2);

    int ikOld =
        kCryst[0] * kMesh(2) * kMesh(1) + kCryst[1] * kMesh(2) + kCryst[2];
    {
      std::string eigFileName = wannierPrefix + ".eig";
      std::ifstream eigfile(eigFileName);
      int ib, ikk;
      double x;
      while (eigfile >> ib >> ikk >> x) {
        if (ikk - 1 == ikOld) {
          // Note: this causes some warnings from Eigen
          blochEnergies(ib - 1, ik) = x;
        }
      }
    }
  }

  //----------------------------------------------------------------------------
  // Now FT to Wannier representation
  Eigen::Tensor<std::complex<double>, 3> h0R(elBravaisVectors.cols(),
                                             numWannier, numWannier);
  h0R.setZero();
  for (int ik1 = 0; ik1 < kPoints.getNumPoints(); ik1++) {
    auto k1C = kPoints.getPointCoords(ik1, Points::cartesianCoords);

    // u has size (numBands, numWannier, numKPoints)
    Eigen::MatrixXcd uK(numBands, numWannier);
    for (int i = 0; i < numBands; i++) {
      for (int j = 0; j < numWannier; j++) {
        uK(i, j) = uMatrices(i, j, ik1);
      }
    }

    // Eq. 26 of Giustino PRB 2007. Note that the U are inverted
    Eigen::MatrixXcd h0K1(numBands, numBands);
    for (int ib = 0; ib < numBands; ib++) {
      h0K1(ib, ib) = {blochEnergies(ib, ik1), 0};
    }
    Eigen::MatrixXcd h0K(numWannier, numWannier);
    h0K = uK.transpose() * h0K1 * uK.adjoint().transpose();

    for (int iR = 0; iR < elBravaisVectors.cols(); iR++) {
      Eigen::Vector3d R = elBravaisVectors.col(iR);
      double arg = k1C.dot(R);
      std::complex<double> phase =
          exp(-complexI * arg) / double(kPoints.getNumPoints());
      for (int m = 0; m < numWannier; m++) {
        for (int n = 0; n < numWannier; n++) {
          h0R(iR, m, n) += phase * h0K(m, n);
        }
      }
    }
  }

  //  --------------------------------------------------------------------------
  // FFT back

  for (int ik = 0; ik < kPoints.getNumPoints(); ik++) {
    // get U
    auto k1C = kPoints.getPointCoords(ik, Points::cartesianCoords);
    auto t = electronH0.diagonalizeFromCoords(k1C);
    auto en = std::get<0>(t);
    auto u = std::get<1>(t);

    Eigen::MatrixXcd h0K(numWannier, numWannier);
    h0K.setZero();
    for (long iR = 0; iR < elBravaisVectors.cols(); iR++) {
      Eigen::Vector3d R = elBravaisVectors.col(iR);
      double arg = k1C.dot(R);
      std::complex<double> phase = exp(complexI * arg) / elDegeneracies(iR);
      for (long m = 0; m < numWannier; m++) {
        for (long n = 0; n < numWannier; n++) {
          h0K(m, n) += phase * h0R(iR, m, n);
        }
      }
    }

    h0K = u.adjoint() * h0K * u;

    for (int ib = 0; ib < numWannier; ib++) {
      assert(abs(h0K(ib, ib).real() - blochEnergies(ib, ik)) < 1.0e-4);
    }
  }
}

void ElPhQeToPhoebeApp::testPhononTransform(
    Crystal &crystal, PhononH0 &phononH0, Points &qPoints,
    const Eigen::Tensor<std::complex<double>, 3> &phEigenvectors,
    const Eigen::MatrixXd &phBravaisVectors,
    const Eigen::VectorXd &phDegeneracies, const Eigen::MatrixXd &phEnergies) {
  /** Like the test above, we
   * 1) FT to Wannier representation.
   *    Since these are force constants, they should be real.
   * 2) FT back to Bloch space and check that we find the same results.
   *
   * We also verify that the eigenvectors are normalized by masses
   * Note that the test works fine for non-polar systems.
   */

  int numPhBands = phononH0.getNumBands();

  // Bloch To Wannier transform

  auto atomicPositions = crystal.getAtomicPositions();
  int numAtoms = atomicPositions.rows();
  auto atomicMasses = crystal.getAtomicMasses();

  // test mass normalization as expected
  for (int iq = 0; iq < qPoints.getNumPoints(); iq++) {
    Eigen::MatrixXcd norm(numPhBands, numPhBands);
    norm.setZero();
    for (int ib1 = 0; ib1 < numPhBands; ib1++) {
      for (int ib2 = 0; ib2 < numPhBands; ib2++) {
        for (int k1 = 0; k1 < numAtoms; k1++) {
          for (int iCart : {0, 1, 2}) {
            int i = compress2Indeces(k1, iCart, numAtoms, 3);
            norm(ib1, ib2) +=
                phEigenvectors(i, ib1, iq) * sqrt(atomicMasses(k1)) *
                phEigenvectors(i, ib2, iq) * sqrt(atomicMasses(k1));
          }
        }
      }
      norm(ib1) -= 1.; // It should be an identity matrix
    }
    assert(norm.sum() < 1.0e-6);
  }

  // FT to Wannier representation

  Eigen::Tensor<std::complex<double>, 5> h0R(
      numAtoms * numAtoms * phBravaisVectors.size(), numAtoms, numAtoms, 3, 3);
  h0R.setZero();

  for (int iq = 0; iq < qPoints.getNumPoints(); iq++) {
    auto qC = qPoints.getPointCoords(iq, Points::cartesianCoords);
    qC = qPoints.bzToWs(qC, Points::cartesianCoords);

    // u has size (numBands, numWannier, numKPoints)
    Eigen::MatrixXcd uK(numPhBands, numPhBands);
    for (int k1 = 0; k1 < numAtoms; k1++) {
      for (int iCart : {0, 1, 2}) {
        int i = compress2Indeces(k1, iCart, numAtoms, 3);
        for (int j = 0; j < numPhBands; j++) {
          uK(i, j) = phEigenvectors(i, j, iq) * sqrt(atomicMasses(k1));
        }
      }
    }
    assert(uK.inverse() == uK.adjoint()); // check is unitary matrix

    // build dynamical matrix
    Eigen::MatrixXcd h0K(numPhBands, numPhBands);
    h0K.setZero();
    for (int ib = 0; ib < numPhBands; ib++) {
      h0K(ib, ib) = {phEnergies(ib, iq) * phEnergies(ib, iq), 0.};
    }
    h0K = uK * h0K * uK.adjoint();
    // if here multiply by mass, we get the QE results

    for (long iR = 0; iR < phBravaisVectors.cols(); iR++) {
      Eigen::Vector3d R0 = phBravaisVectors.col(iR);
      for (int k1 = 0; k1 < numAtoms; k1++) {
        for (int k2 = 0; k2 < numAtoms; k2++) {
          Eigen::Vector3d R = R0; // - atomicPositions.col(k1)
          //+ atomicPositions.col(k2);
          double arg = qC.dot(R);
          std::complex<double> phase =
              exp(-complexI * arg) / double(qPoints.getNumPoints());
          for (int iCart : {0, 1, 2}) {
            for (int jCart : {0, 1, 2}) {
              int m = compress2Indeces(k1, iCart, numAtoms, 3);
              int n = compress2Indeces(k2, jCart, numAtoms, 3);
              h0R(iR, k1, k2, iCart, jCart) += phase * h0K(m, n);
            }
          }
        }
      }
    }
  }

  // check that h0R, the force constants, are real
  {
    double realSum = 0.;
    double imagSum = 0.;
    for (long iR0 = 0; iR0 < phBravaisVectors.cols(); iR0++) {
      for (int k1 = 0; k1 < numAtoms; k1++) {
        for (int k2 = 0; k2 < numAtoms; k2++) {
          for (int i : {0, 1, 2}) {
            for (int j : {0, 1, 2}) {
              double x = std::real(h0R(iR0, k1, k2, i, j));
              realSum += pow(x, 2);
              imagSum += pow(std::imag(h0R(iR0, k1, k2, i, j)), 2);
              // set to zero the imag part to clean noise
              // this is also what QE does
              h0R(iR0, k1, k2, i, j) = {x, 0.};
            }
          }
        }
      }
    }
    // I want the imag part to be much smaller than the real
    assert(imagSum * pow(10, 6) < realSum);
  }

  //--------------------------------------------------------------------------
  // FFT back

  for (int iq = 0; iq < qPoints.getNumPoints(); iq++) {
    // get U
    auto qC = qPoints.getPointCoords(iq, Points::cartesianCoords);
    auto t = phononH0.diagonalizeFromCoords(qC, false);
    auto en = std::get<0>(t);
    auto u = std::get<0>(t);

    Eigen::MatrixXcd hWK(numPhBands, numPhBands);
    hWK.setZero();
    for (long iR = 0; iR < phBravaisVectors.cols(); iR++) {
      Eigen::Vector3d R0 = phBravaisVectors.col(iR);
      for (int k1 = 0; k1 < numAtoms; k1++) {
        for (int k2 = 0; k2 < numAtoms; k2++) {
          Eigen::Vector3d R = R0; // - atomicPositions.col(k1)
          //+ atomicPositions.col(k2);
          double arg = qC.dot(R);
          std::complex<double> phase = exp(complexI * arg) / phDegeneracies(iR);
          for (int iCart : {0, 1, 2}) {
            for (int jCart : {0, 1, 2}) {
              int m = compress2Indeces(k1, iCart, numAtoms, 3);
              int n = compress2Indeces(k2, jCart, numAtoms, 3);
              hWK(m, n) += phase * h0R(iR, k1, k2, iCart, jCart);
            }
          }
        }
      }
    }

    // diagonalize it, using the matrices from phononH0
    auto dq = u.adjoint() * hWK * u;
    (void)dq;
    // check I found again the same eigenvalues
    for (int ib = 0; ib < numPhBands; ib++) {
      assert(abs(std::sqrt(dq(ib, ib).real()) - phEnergies(ib, iq)) < 1.0e-6);
    }
  }
}

void ElPhQeToPhoebeApp::testBackTransform(
    Context &context, PhononH0 &phononH0, Points &kPoints, Points &qPoints,
    ElectronH0Wannier &electronH0, Crystal &crystal,
    Eigen::Tensor<std::complex<double>, 5> &gFull) {
  /** This is the important test of el-ph Wannier interpolation
   * We compute the band structure
   * Read the el-ph interaction from file
   * Check that the el-ph coupling, interpolated on the same initial grid,
   * is the same of the el-ph coupling read from QE.
   */
  bool withVelocities = true;
  bool withEigenvectors = true;
  FullBandStructure bandStructure =
      electronH0.populate(kPoints, withVelocities, withEigenvectors);
  int numKPoints = kPoints.getNumPoints();
  int numModes = phononH0.getNumBands();

  // needed by ::parse()
  context.setEpwFileName(context.getQuantumEspressoPrefix() + ".phoebe.elph.dat");
  auto couplingElPh = InteractionElPhWan::parse(context, crystal, &phononH0);

  for (long ik1 = 0; ik1 < numKPoints; ik1++) {
    Eigen::Vector3d k1C = kPoints.getPointCoords(ik1, Points::cartesianCoords);
    for (long ik2 = 0; ik2 < numKPoints; ik2++) {
      Eigen::Vector3d k2C =
          kPoints.getPointCoords(ik2, Points::cartesianCoords);

      std::vector<Eigen::Vector3d> k2Cs;
      k2Cs.push_back(k2C);

      Eigen::Vector3d q3C = k2C - k1C;
      Eigen::Vector3d q3Cryst = qPoints.cartesianToCrystal(q3C);
      int iq3 = qPoints.getIndex(q3Cryst);
      std::vector<Eigen::Vector3d> q3Cs;
      q3Cs.push_back(q3C);

      auto ik1Index = WavevectorIndex(ik1);
      auto ik2Index = WavevectorIndex(ik2);

      Eigen::MatrixXcd eigvec1 = bandStructure.getEigenvectors(ik1Index);
      Eigen::MatrixXcd eigvec2 = bandStructure.getEigenvectors(ik2Index);
      std::vector<Eigen::MatrixXcd> eigvecs2;
      eigvecs2.push_back(eigvec2);

      auto t = phononH0.diagonalizeFromCoords(q3C);
      auto eigvec3 = std::get<1>(t);
      std::vector<Eigen::MatrixXcd> eigvecs3;
      eigvecs3.push_back(eigvec3);

      couplingElPh.calcCouplingSquared(eigvec1, eigvecs2, eigvecs3, k1C, k2Cs,
                                       q3Cs);
      auto coupling2 = couplingElPh.getCouplingSquared(0);

      double sum1 = 0.;
      double sum2 = 0.;
      for (int ib1 = 0; ib1 < 4; ib1++) {
        for (int ib2 = 0; ib2 < 4; ib2++) {
          for (int ib3 = 0; ib3 < numModes; ib3++) {
            sum1 += std::norm(gFull(ib2, ib1, ib3, ik1, iq3));
            sum2 += coupling2(ib1, ib2, ib3);
          }
        }
      }
      // note that I change the meaning of the indeces
      assert(abs((sum1 - sum2) / sum1) < 0.0001);
    }
  }
}

void ElPhQeToPhoebeApp::postProcessingWannier(
    Context &context, Crystal &crystal, PhononH0 &phononH0, FullPoints &kPoints,
    FullPoints &qPoints, int numQEBands, int numModes, int numIrrQPoints,
    int numElectrons, int numSpin, const Eigen::MatrixXd &energies,
    const Eigen::MatrixXd &kgridFull, const Eigen::Vector3i &kMesh,
    const Eigen::Vector3i &qMesh, bool runTests) {
  if (mpi->mpiHead()) {
    std::cout << "Starting Wannier post-processing\n" << std::endl;
  }

  std::string wannierPrefix = context.getWannier90Prefix();

  auto t2 = QEParser::parseElHarmonicWannier(context, &crystal);
  auto electronH0 = std::get<1>(t2);
  int numWannier = electronH0.getNumBands();

  //----------------------------------------------------------------------------
  // read Wannier90 rotation matrices

  Eigen::Tensor<std::complex<double>, 3> uMatrices;
  // uMatrices has size (numBands, numWannier, numKPoints)
  uMatrices = setupRotationMatrices(wannierPrefix, kPoints);
  int numBands = uMatrices.dimension(0);        // number of entangled bands
  assert(numWannier == uMatrices.dimension(1)); // number of entangled bands

  //----------------------------------------------------------------------------

  // read coupling from file
  auto t5 =
      readGFromQEFile(context, numModes, numBands, numWannier, kPoints, qPoints,
                      kgridFull, numIrrQPoints, numQEBands, energies);
  auto gFull = std::get<0>(t5);          // (nBands,nBands,nModes,numK,numQ)
  auto phEigenvectors = std::get<1>(t5); // (numModes,numModes,numQPoints)
  auto phEnergies = std::get<2>(t5);     // (numModes,numQPoints)

  //----------------------------------------------------------------------------

  // Find the lattice vectors for the Fourier transforms

  auto t3 = crystal.buildWignerSeitzVectors(kMesh);
  Eigen::MatrixXd elBravaisVectors = std::get<0>(t3);
  Eigen::VectorXd elDegeneracies = std::get<1>(t3);

  auto t4 = crystal.buildWignerSeitzVectors(qMesh);
  Eigen::MatrixXd phBravaisVectors = std::get<0>(t4);
  Eigen::VectorXd phDegeneracies = std::get<1>(t4);

  Eigen::Tensor<std::complex<double>, 5> gWannier =
      blochToWannier(elBravaisVectors, phBravaisVectors, gFull, uMatrices,
                     phEigenvectors, kPoints, qPoints, crystal, phononH0);

  //--------------------------------------------------------------------------

  // Dump el-ph in Wannier representation to file

  if (mpi->mpiHead()) {
    std::cout << "Start writing g to file" << std::endl;
    std::string phoebePrefixQE = context.getQuantumEspressoPrefix();
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
    for (auto x : gWannier.dimensions()) {
      outfile << x << " ";
    }
    outfile << "\n";

    outfile << std::setprecision(16);
    int numPhBands = 3 * crystal.getNumAtoms();
    for (int i5 = 0; i5 < elDegeneracies.size(); i5++) {
      for (int i4 = 0; i4 < phDegeneracies.size(); i4++) {
        for (int i3 = 0; i3 < numPhBands; i3++) {
          for (int i2 = 0; i2 < numWannier; i2++) {
            for (int i1 = 0; i1 < numWannier; i1++) {
              outfile << std::setw(22) << gWannier(i1, i2, i3, i4, i5).real()
                      << " " << std::setw(22)
                      << gWannier(i1, i2, i3, i4, i5).imag() << "\n";
            }
          }
        }
      }
    }
    std::cout << "Done writing g to file\n" << std::endl;
  }

  if (runTests) {
    testElectronicTransform(kPoints, wannierPrefix, elBravaisVectors, uMatrices,
                            elDegeneracies, electronH0);

    testPhononTransform(crystal, phononH0, qPoints, phEigenvectors,
                        phBravaisVectors, phDegeneracies, phEnergies);

    testBackTransform(context, phononH0, kPoints, qPoints, electronH0, crystal,
                      gFull);
  }
}
