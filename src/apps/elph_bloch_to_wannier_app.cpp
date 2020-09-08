#include "elph_bloch_to_wannier_app.h"
#include "bandstructure.h"
#include "eigen.h"
#include "qe_input_parser.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

int computeBandOffset(const Eigen::VectorXd &energiesFull,
                      const Eigen::VectorXd &energiesWannier) {
  // energies Wannier are the entangled bands

  int numFull = energiesFull.size();
  int numWannier = energiesWannier.size();

  // we find the offset by comparing the energy differences
  // the offset which minimizes energy differences is the chosen one
  int possibleValues = numFull - numWannier + 1;
  Eigen::VectorXd difference(possibleValues);
  difference.setZero();
  for (int i = 0; i < possibleValues; i++ ) {
    for (int ib=0; ib<numWannier; ib++) {
      difference(i) = pow(energiesFull(ib) - energiesWannier(ib+i),2);
    }
  }

  // offset = index of min difference
  int offset = -1;
  for ( int i=0; i<possibleValues; i++ ) {
    if ( difference(i) == difference.minCoeff() ) {
      offset = i;
      break;
    }
  }

  return offset;
}

void ElPhBlochToWannierApp::run(Context &context) {
  (void)context;

  std::string phoebePrefixQE = "silicon";
  std::string wannierPrefix = "si";

  // read Hamiltonian of phonons and electrons

  // actually, we only need the crystal
  auto t1 = QEParser::parsePhHarmonic(context);
  auto crystal = std::get<0>(t1);
  auto phononH0 = std::get<1>(t1);

  auto t2 = QEParser::parseElHarmonicWannier(context, &crystal);
  auto crystalEl = std::get<0>(t2);
  auto electronH0 = std::get<1>(t2);

  if (phoebePrefixQE.empty()) {
    Error e("Must provide an input H0 file name");
  }
  // open the 0000 input file
  std::string fileName = phoebePrefixQE + ".phoebe.0000.dat";
  std::ifstream infile(fileName);
  std::string line;
  if (not infile.is_open()) {
    Error e("H0 file not found");
  }
  std::getline(infile, line); // first line is a title
  int numQEBands; // number of Kohn-Sham states
  double numElectrons; // number of electrons (spin degeneracy included)
  int numSpin; // should always be one, without support for spin
  infile >> numQEBands >> numElectrons >> numSpin;
  Eigen::Vector3i kMesh, qMesh;
  infile >> qMesh(0) >> qMesh(1) >> qMesh(2)
      >> kMesh(0) >> kMesh(1) >> kMesh(2);
  int bogusI;
  double bogusD;
  int numAtoms;
  infile >> bogusD >> numAtoms; // alat and nat
  assert(numAtoms == crystal.getNumAtoms());
  // unit cell
  for ( int i=0; i<9; i++ ) {
    infile >> bogusD;
  }
  // reciprocal unit cell
  for ( int i=0; i<9; i++ ) {
    infile >> bogusD;
  }
  // ityp
  for ( int i=0; i<numAtoms; i++ ) {
    infile >> bogusI;
  }
  // positions
  for ( int i=0; i<3*numAtoms; i++ ) {
    infile >> bogusD;
  }

  int numQPoints, numIrrQPoints;
  infile >> numQPoints >> numIrrQPoints;

  Eigen::MatrixXd qgridFull(3,numQPoints);
  for ( int iq=0; iq<numQPoints; iq++ ) {
    infile >> qgridFull(0,iq) >> qgridFull(1,iq) >> qgridFull(2,iq);
  }

  int numKPoints;
  infile >> numKPoints;
  Eigen::MatrixXd kgridFull(3,numKPoints);
  for ( int ik=0; ik<numKPoints; ik++ ) {
    infile >> kgridFull(0,ik) >> kgridFull(1,ik) >> kgridFull(2,ik);
  }

  Eigen::MatrixXd energies(numQEBands,numKPoints);
  for ( int ik=0; ik<numKPoints; ik++ ) {
    for ( int ib=0; ib<numQEBands; ib++ ) {
      infile >> energies(ib, ik);
    }
  }

  infile.close();

  //-------------------------------------
  // read Wannier90 rotation matrices

  FullPoints kPoints(crystal, kMesh);
  FullPoints qPoints(crystal, qMesh);

  Eigen::Tensor<std::complex<double>, 3>
  uMatrices = setupRotationMatrices(wannierPrefix, kPoints);

  // compute the band offset, to align QE bands with the entangled bands
  // to use in input to Wannier90

  int bandsOffset;
  {
    Eigen::VectorXd energiesQEAtZero = energies.col(0); // k = 0

    std::string fileName = wannierPrefix + ".nnkp";
    std::ifstream infile(fileName);
    for ( int i=0; i<18; i++ ) {
      std::getline(infile, line); // skip the first 18 lines
    }
    double kx, ky, kz;
    infile >> kx >> ky >> kz;
    assert(kx*kx+ky*ky+kz*kz < 1.0e-5); // check first point is gamma
    infile.close();

    // read .eig file to get energies

    std::vector<double> ens;
    std::string eigFileName = wannierPrefix + ".eig";
    std::ifstream eigfile(eigFileName);
    int ib, ik;
    double x;
    while ( eigfile >> ib >> ik >> x ) {
      if (ik > 1) {
        break;
      }
      ens.push_back(x);
    }
    eigfile.close();
    int nb = ens.size();
    Eigen::VectorXd energiesWannierAtZero(nb);
    for ( ib =0; ib<nb; ib++) {
      energiesWannierAtZero(ib) = ens[ib];
    }

    bandsOffset = computeBandOffset(energiesQEAtZero, energiesWannierAtZero);
  }

  //---------------------------------------------------

  // read coupling from file

  int numModes = 3 * numAtoms;
  int numBands = uMatrices.dimension(0); // number of entangled bands
  int numWannier = uMatrices.dimension(1); // number of entangled bands
  Eigen::Tensor<std::complex<double>, 5> g_full(numBands, numBands, numModes,
                                                numKPoints, numQPoints);
  Eigen::Tensor<std::complex<double>, 3> phEigenvectors(numModes, numModes, numQPoints);
  g_full.setZero();
  phEigenvectors.setZero();

  // g_full is computed on full grid

  Eigen::VectorXi ikMap(numKPoints);
  for (long ikOld=0; ikOld<numKPoints; ikOld++) {
    Eigen::Vector3d kOld = kgridFull.col(ikOld);
    long ikNew = kPoints.getIndex(kOld);
    ikMap(ikOld) = ikNew;
  }

  for (int iqIrr=0; iqIrr<numIrrQPoints; iqIrr++ ) {
    std::stringstream ss;
    ss << std::setw(4) << std::setfill('0') << iqIrr;
    std::string numString = ss.str();
    std::string fileName = phoebePrefixQE + ".phoebe." + numString + ".dat";
    std::ifstream infile(fileName);

    int nqStar;
    infile >> nqStar;
    std::vector<Eigen::Vector3d> qStar;
    for (int iq = 0; iq < nqStar; iq++) {
      Eigen::Vector3d thisQ; // in crystal coordinates
      infile >> thisQ(0) >> thisQ(1) >> thisQ(2);
      qStar.push_back(thisQ);
    }

    Eigen::VectorXd phononEnergies(numModes);
    for (int nu = 0; nu < numModes; nu++) {
      infile >> phononEnergies(nu);
    }

    Eigen::MatrixXd phononEigenvectors(numModes, numModes);
    for (int j = 0; j < numModes; j++) {
      for (int i = 0; i < numModes; i++) {
        infile >> phononEigenvectors(i, j);
      }
    }

    Eigen::Tensor<std::complex<double>, 5> thisG(numQEBands, numQEBands,
                                                 numModes, numKPoints, nqStar);
    thisG.setZero();
    for (int iq = 0; iq < nqStar; iq++) {
      for (int nu = 0; nu < numModes; nu++) {
        for (int ik = 0; ik < numKPoints; ik++) {
          for (int ib2 = 0; ib2 < numQEBands; ib2++) {
            for (int ib1 = 0; ib1 < numQEBands; ib1++) {
              double re, im;
              infile >> re >> im;
              thisG(ib1,ib2,nu,ik,iq) = {re, im};
            }
          }
        }
      }
    }
    infile.close();

    for ( int iqStar=0; iqStar<nqStar; iqStar++ ) {
      Eigen::Vector3d qVec = qStar[iqStar];
      long iqFull = qPoints.getIndex(qVec);
      for (int nu = 0; nu < numModes; nu++) {
        for (int ik = 0; ik < numKPoints; ik++) {
          for (int ib2 = 0; ib2 < numWannier; ib2++) {
            for (int ib1 = 0; ib1 < numWannier; ib1++) {
              g_full(ib1, ib2, nu, ikMap(ik), iqFull) =
                thisG(bandsOffset+ib1, bandsOffset+ib2, nu, ik, iqStar);
            }
          }
        }
      }
      for (int j = 0; j < numModes; j++) {
        for (int i = 0; i < numModes; i++) {
          phEigenvectors(i, j, iqFull) = phononEigenvectors(i, j);
        }
      }
    }
  }

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
                     uMatrices, phEigenvectors, kPoints, qPoints);

  // Now I can dump g_wannier to file
  // write the total number of electrons to file!

  {
    std::string outFileName = phoebePrefixQE + ".phoebe.elph.dat";
    std::ofstream outfile(outFileName);
    if (not outfile.is_open()) {
      Error e("Output file couldn't be opened");
    }
    outfile << numElectrons << numSpin << "\n";
    outfile << kMesh << "\n";
    outfile << qMesh << "\n";
    outfile << phBravaisVectors.rows() << phBravaisVectors.cols() << "\n";
    outfile << phBravaisVectors << "\n";
    outfile << elBravaisVectors.rows() << elBravaisVectors.cols() << "\n";
    outfile << elBravaisVectors << "\n";
    outfile << g_wannier.dimensions() << "\n";
    outfile << g_wannier << "\n";
    outfile.close();
  }
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
    Eigen::Tensor<std::complex<double>, 5> &g_full,
    const Eigen::Tensor<std::complex<double>,3> &uMatrices,
    const Eigen::Tensor<std::complex<double>,3> &phEigenvectors,
    FullPoints &kPoints, FullPoints &qPoints) {

  int numBands = g_full.dimension(0); // # of entangled bands
  int numModes = g_full.dimension(2);
  int numKPoints = g_full.dimension(3);
  int numQPoints = g_full.dimension(4);
  int numElBravaisVectors = elBravaisVectors.cols();
  int numPhBravaisVectors = phBravaisVectors.cols();
  int numWannier = uMatrices.dimension(1);

  std::array<Eigen::Index, 5> zeros;
  //  zeros.data() << 0,0,0,0,0;
  for (auto &s : zeros) {
    s = 0;
  }

  Eigen::Tensor<std::complex<double>, 5> g_full_tmp(
      numBands, numBands, numModes, numKPoints, numQPoints);
  g_full_tmp.setZero();
  for (long iq = 0; iq < numQPoints; iq++) {
    Eigen::Vector3d q = qPoints.getPointCoords(iq, Points::cartesianCoords);
    for (long ik = 0; ik < numKPoints; ik++) {
      Eigen::Vector3d k = kPoints.getPointCoords(ik, Points::cartesianCoords);

      // Coordinates and index of k+q point
      Eigen::Vector3d kq = k + q;
      Eigen::Vector3d kqCrystal = kPoints.cartesianToCrystal(kq);
      long ikq = kPoints.getIndex(kqCrystal);

      // First we transform from the Bloch to Wannier Gauge

      Eigen::MatrixXcd uK(numWannier,numBands);
      Eigen::MatrixXcd uKq(numWannier,numBands);
      for (int i = 0; i < numBands; i++) {
        for (int j = 0; j < numWannier; j++) {
          uK(i,j) = uMatrices(i,j,ik);
          uKq(i,j) = uMatrices(i,j,ikq);
        }
      }
      uKq = uKq.adjoint();

      Eigen::Tensor<std::complex<double>, 3> tmp(numWannier, numBands, numModes);
      tmp.setZero();
      for (int nu = 0; nu < numModes; nu++) {
        for (int i = 0; i < numWannier; i++) {
          for (int j = 0; j < numBands; j++) {
            for (int l = 0; l < numBands; l++) {
              tmp(i, j, nu) += uKq(i, l) * g_full(l, j, nu, ik, iq);
            }
          }
        }
      }
      for (int nu = 0; nu < numModes; nu++) {
        for (int i = 0; i < numWannier; i++) {
          for (int j = 0; j < numWannier; j++) {
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
    for (long ik = 0; ik < numKPoints; ik++) {
      Eigen::Vector3d k = kPoints.getPointCoords(ik, Points::cartesianCoords);
      double arg = k.dot(elBravaisVectors.col(iR));
      std::complex<double> phase = exp(-complexI * arg) / double(numKPoints);
      for (int iq = 0; iq < numQPoints; iq++) {
        for (int i = 0; i < numWannier; i++) {
          for (int j = 0; j < numWannier; j++) {
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
  for (long iq = 0; iq < numQPoints; iq++) {

    Eigen::MatrixXcd uQ(numModes,numModes);
    for (int nu = 0; nu < numModes; nu++) {
      for (int nu2 = 0; nu2 < numModes; nu2++) {
        uQ(nu,nu2) = phEigenvectors(nu,nu2,iq);
      }
    }
    uQ = uQ.inverse();

    for (int nu = 0; nu < numModes; nu++) {
      for (int nu2 = 0; nu2 < numModes; nu2++) {
        for (int irE = 0; irE < numElBravaisVectors; irE++) {
          for (int i = 0; i < numWannier; i++) {
            for (int j = 0; j < numWannier; j++) {
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
    Eigen::Vector3d q = qPoints.getPointCoords(iq, Points::cartesianCoords);
    for (int irP = 0; irP < numPhBravaisVectors; irP++) {
      double arg = q.dot(phBravaisVectors.col(irP));
      std::complex<double> phase = exp(-complexI * arg) / double(numQPoints);
      for (int irE = 0; irE < numElBravaisVectors; irE++) {
        for (int i = 0; i < numWannier; i++) {
          for (int j = 0; j < numWannier; j++) {
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
