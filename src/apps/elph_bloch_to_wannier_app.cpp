#include "bandstructure.h"
#include "electron_h0_wannier.h"
#include "elph_bloch_to_wannier_app.h"
#include "phonon_h0.h"
#include "eigen.h"
#include "qe_input_parser.h"

void ElPhBlochToWannierApp::run(Context &context) {
  (void) context;

  auto t2 = QEParser::parsePhHarmonic(context);
  auto crystal = std::get<0>(t2);
  auto phononH0 = std::get<1>(t2);

  auto t1 = QEParser::parseElHarmonicWannier(context, &crystal);
  auto crystalEl = std::get<0>(t1);
  auto electronH0 = std::get<1>(t1);

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
  // probably U matrices are rectangles

  //----------------------------------------------------
  // Find the lattice vectors for the Fourier transforms

  Eigen::MatrixXd elBravaisVectors = getBravaisVectors(crystal, kMesh);
  Eigen::MatrixXd phBravaisVectors = getBravaisVectors(crystal, qMesh);

  //--------------------------------------------------------
  // Start the Bloch to Wannier transformation

  Eigen::Tensor<std::complex<double>, 5> g_wannier =
      blochToWannier(elBravaisVectors, phBravaisVectors, g_full,
                     elBandStructure, phBandStructure);

  // Now I can dump g_wannier to file

  mpi->barrier();
}

void ElPhBlochToWannierApp::checkRequirements(Context &context) {
  (void) context;
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
  for (auto& s: zeros) s = 0;

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
  Eigen::Tensor<std::complex<double>,5> g_mixed(numBands, numBands, numModes,
                                                numElBravaisVectors, numQPoints);
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

  Eigen::Tensor<std::complex<double>,5> g_wannier_tmp(
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

  Eigen::Tensor<std::complex<double>,5> g_wannier(
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

Eigen::MatrixXd ElPhBlochToWannierApp::getBravaisVectors(Crystal &crystal,
                                   const Eigen::Vector3i &kMesh) {
  Eigen::MatrixXd bravaisVectors(3,0);
  bravaisVectors.setZero();
  return bravaisVectors;
}


