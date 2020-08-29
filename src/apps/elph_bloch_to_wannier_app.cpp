#include "elph_bloch_to_wannier_app.h"
#include "bandstructure.h"
#include "eigen.h"
#include "electron_h0_wannier.h"
#include "phonon_h0.h"
#include "qe_input_parser.h"
#include <fstream>
#include <string>

void ElPhBlochToWannierApp::run(Context &context) {
  (void)context;

  auto t1 = QEParser::parsePhHarmonic(context);
  auto crystal = std::get<0>(t1);
  auto phononH0 = std::get<1>(t1);

  auto t2 = QEParser::parseElHarmonicWannier(context, &crystal);
  auto crystalEl = std::get<0>(t2);
  auto electronH0 = std::get<1>(t2);

  bool withVelocities = false;
  bool withEigenvectors = true;
  // we need the eigenvectors to perform the Bloch to Wannier rotation
  Eigen::Vector3i qMesh = phononH0.getCoarseGrid();
  FullPoints qPoints(crystal, qMesh);
  FullBandStructure phBandStructure =
      phononH0.populate(qPoints, withVelocities, withEigenvectors);
  Eigen::Vector3i kMesh;
  kMesh << 1, 1, 1;
  FullPoints kPoints(crystal, kMesh);
  FullBandStructure elBandStructure =
      electronH0.populate(kPoints, withVelocities, withEigenvectors);

  // read input file, and the Bloch representation of the el-ph coupling

  // Unfold the symmetries

  int numBands = 0; //????????????????????????????????????????
  int numModes = phBandStructure.getNumBands();
  int numKPoints = elBandStructure.getNumPoints();
  int numQPoints = phBandStructure.getNumPoints();
  Eigen::Tensor<std::complex<double>, 5> g_full(numBands, numBands, numModes,
                                                numKPoints, numQPoints);
  // g_full should come from the unsymmetrization

  // Note: I should check the disentanglement of Wannier bands,
  // probably U matrices are rectangles, and not squares, and can't be computed
  // from the hamiltonian of the disentangled subspace

  //----------------------------------------------------
  // Find the lattice vectors for the Fourier transforms

  auto t3 = crystal.buildWignerSeitzVectors(kMesh);
  Eigen::MatrixXd elBravaisVectors = std::get<0>(t3);
  Eigen::VectorXd elDegeneracies = std::get<1>(t3);

  auto t4 = crystal.buildWignerSeitzVectors(qMesh);
  Eigen::MatrixXd phBravaisVectors = std::get<0>(t4);
  Eigen::VectorXd phDegeneracies = std::get<1>(t4);

  //--------------------------------------------------------
  // Start the Bloch to Wannier transformation

  Eigen::Tensor<std::complex<double>, 5> g_wannier =
      blochToWannier(elBravaisVectors, phBravaisVectors, g_full,
                     elBandStructure, phBandStructure);

  // Now I can dump g_wannier to file

  mpi->barrier();
}

void ElPhBlochToWannierApp::checkRequirements(Context &context) {
  (void)context;
  //  throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");
  //  throwErrorIfUnset(context.getKMesh(), "kMesh");
  //  throwErrorIfUnset(context.getEpwFileName(), "EpwFileName");
  //  throwErrorIfUnset(context.getTemperatures(), "temperatures");
  //  throwErrorIfUnset(context.getSmearingMethod(), "smearingMethod");
  //  if ( context.getSmearingMethod() == DeltaFunction::gaussian) {
  //    throwErrorIfUnset(context.getSmearingWidth(), "smearingWidth");
  //  }
  //
  //  if ( context.getDopings().size() == 0 &&
  //      context.getChemicalPotentials().size() == 0) {
  //    Error e("Either chemical potentials or dopings must be set");
  //  }
}

Eigen::Tensor<std::complex<double>, 5> ElPhBlochToWannierApp::blochToWannier(
    const Eigen::MatrixXd &elBravaisVectors,
    const Eigen::MatrixXd &phBravaisVectors,
    const Eigen::Tensor<std::complex<double>, 5> &g_full,
    FullBandStructure &elBandStructure, FullBandStructure &phBandStructure) {

  int numQPoints = phBandStructure.getNumPoints();
  int numKPoints = elBandStructure.getNumPoints();

  auto kPoints = elBandStructure.getPoints();
  auto qPoints = phBandStructure.getPoints();

  int numElBravaisVectors = elBravaisVectors.cols();
  int numPhBravaisVectors = phBravaisVectors.cols();

  int numBands = elBandStructure.getNumBands();
  int numModes = phBandStructure.getNumBands();

  std::array<Eigen::Index, 5> zeros;
  //  zeros.data() << 0,0,0,0,0;
  for (auto &s : zeros)
    s = 0;

  Eigen::Tensor<std::complex<double>, 5> g_full_tmp(
      numBands, numBands, numModes, numKPoints, numQPoints);
  g_full_tmp.setZero();
  for (int iq = 0; iq < numQPoints; iq++) {
    auto iqIndex = WavevectorIndex(iq);
    Eigen::Vector3d q = elBandStructure.getWavevector(iqIndex);
    for (int ik = 0; ik < numKPoints; ik++) {
      auto ikIndex = WavevectorIndex(ik);
      Eigen::Vector3d k = elBandStructure.getWavevector(ikIndex);

      // Coordinates and index of k+q point
      Eigen::Vector3d kq = k + q;
      Eigen::Vector3d kqCrystal = kPoints.cartesianToCrystal(kq);
      int ikq = kPoints.getIndex(kqCrystal);
      auto ikqIndex = WavevectorIndex(ikq);

      // First we transform from the Bloch to Wannier Gauge

      Eigen::MatrixXcd uK = elBandStructure.getEigenvectors(ikIndex);
      Eigen::MatrixXcd uKq = elBandStructure.getEigenvectors(ikqIndex);
      uKq = uKq.adjoint();

      Eigen::Tensor<std::complex<double>, 3> tmp(numBands, numBands, numModes);
      tmp.setZero();
      for (int nu = 0; nu < numModes; nu++) {
        for (int i = 0; i < numBands; i++) {
          for (int j = 0; j < numBands; j++) {
            for (int l = 0; l < numBands; l++) {
              tmp(i, j, nu) += uKq(i, l) * g_full(l, j, nu, ik, iq);
            }
          }
        }
      }
      for (int nu = 0; nu < numModes; nu++) {
        for (int i = 0; i < numBands; i++) {
          for (int j = 0; j < numBands; j++) {
            for (int l = 0; l < numBands; l++) {
              g_full_tmp(i, j, nu, ik, iq) += tmp(i, l, nu) * uK(l, j);
            }
          }
        }
      }
    } // ik
  }   // iq
  g_full.reshape(zeros);

  // Fourier transform on the electronic coordinates
  Eigen::Tensor<std::complex<double>, 5> g_mixed(
      numBands, numBands, numModes, numElBravaisVectors, numQPoints);
  g_mixed.setZero();
  for (int iR = 0; iR < numElBravaisVectors; iR++) {
    for (int ik = 0; ik < numKPoints; ik++) {
      auto ikIndex = WavevectorIndex(ik);
      Eigen::Vector3d k = elBandStructure.getWavevector(ikIndex);
      double arg = k.dot(elBravaisVectors.col(iR));
      std::complex<double> phase = exp(-complexI * arg) / double(numKPoints);
      for (int iq = 0; iq < numQPoints; iq++) {
        for (int i = 0; i < numBands; i++) {
          for (int j = 0; j < numBands; j++) {
            for (int nu = 0; nu < numModes; nu++) {
              g_mixed(i, j, nu, iR, iq) += g_full_tmp(i, j, nu, ik, iq) * phase;
            }
          }
        }
      }
    }
  } // iq
  g_full_tmp.reshape(zeros);

  Eigen::Tensor<std::complex<double>, 5> g_wannier_tmp(
      numBands, numBands, numModes, numElBravaisVectors, numQPoints);
  g_wannier_tmp.setZero();
  for (int iq = 0; iq < numQPoints; iq++) {
    auto iqIndex = WavevectorIndex(iq);
    Eigen::MatrixXcd uQ = phBandStructure.getEigenvectors(iqIndex);
    uQ = uQ.inverse();
    for (int nu = 0; nu < numModes; nu++) {
      for (int nu2 = 0; nu2 < numModes; nu2++) {
        for (int irE = 0; irE < numElBravaisVectors; irE++) {
          for (int i = 0; i < numBands; i++) {
            for (int j = 0; j < numBands; j++) {
              g_wannier_tmp(i, j, nu, irE, iq) +=
                  g_mixed(i, j, nu2, irE, iq) * uQ(nu2, nu);
            }
          }
        }
      }
    }
  }
  g_mixed.reshape(zeros);

  Eigen::Tensor<std::complex<double>, 5> g_wannier(
      numBands, numBands, numModes, numElBravaisVectors, numPhBravaisVectors);
  g_wannier.setZero();
  for (int iq = 0; iq < numQPoints; iq++) {
    auto iqIndex = WavevectorIndex(iq);
    Eigen::Vector3d q = phBandStructure.getWavevector(iqIndex);
    for (int irP = 0; irP < numPhBravaisVectors; irP++) {
      double arg = q.dot(phBravaisVectors.col(irP));
      std::complex<double> phase = exp(-complexI * arg) / double(numQPoints);
      for (int irE = 0; irE < numElBravaisVectors; irE++) {
        for (int i = 0; i < numBands; i++) {
          for (int j = 0; j < numBands; j++) {
            for (int nu = 0; nu < numModes; nu++) {
              g_wannier(i, j, nu, irE, irP) +=
                  phase * g_wannier_tmp(i, j, nu, irE, iq);
            }
          }
        }
      }
    }
  }
  g_wannier_tmp.reshape(zeros);
  return g_wannier;
}

Eigen::Tensor<std::complex<double>, 3>
ElPhBlochToWannierApp::setupRotationMatrices(const std::string &wannierPrefix,
                                             FullPoints &fullPoints) {
  std::string line;

  if (wannierPrefix.empty()) {
    Error e("Must provide an input H0 file name");
  }

  std::string fileName = wannierPrefix + "_u.dat";

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

    double x, y, z;
    infile >> x >> y >> z;
    Eigen::Vector3d thisK;
    thisK << x, y, z; // vector in crystal coords

    int ikk = fullPoints.getIndex(thisK);

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

  std::string fileNameDis = wannierPrefix + "_u_dis.dat";

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

    int ikk = fullPoints.getIndex(thisK);

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
