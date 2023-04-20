#include "elel_to_phoebe.h"
#include "bandstructure.h"
#include "eigen.h"
#include "elph_qe_to_phoebe_app.h"
#include "io.h"
#include "parser.h"
#include "qe_input_parser.h"
#include <iomanip>
#include <sstream>
#include <string>

#ifdef HDF5_AVAIL
#include <highfive/H5Easy.hpp>
#endif

void ElElToPhoebeApp::run(Context &context) {

  #ifndef HDF5_AVAIL
  Error("In order to run the el-el to Phoebe app, you must build Phoebe with HDF5.");
  #endif

  // read in the electron Wannier files and set up the el_h0 hamiltonian
  // crystal has to be supplied in input file, because Wannier90 doesn't
  // print it to output
  auto t1 = Parser::parseElHarmonicWannier(context);
  auto crystal = std::get<0>(t1);
  auto electronH0 = std::get<1>(t1);
  // set up k-points
  Eigen::Vector3i kMesh = context.getKMesh();
  Points kPoints(crystal, kMesh);
  int numK = kMesh.prod();

  // yambo's points are low precision, so we need to lower
  // the tolerance on checks if kpoints will be the same
  double tolerance = 1e-5;

  // now we parse the first header, BS_head_Q1
  Eigen::MatrixXd yamboQPoints;
  int numBands, bandOffset, numPoints; // TODO is numPoints duplicate with numK?
  if (mpi->mpiHead()) {

    std::string yamboPrefix = context.getYamboInteractionPrefix();
    // Q1 should always exist
    std::string fileName = yamboPrefix + "BS_head_Q1.hdf5";
    // Open the hdf5 file
    HighFive::File file(fileName, HighFive::File::ReadOnly);
    // Set up hdf5 datasets
    HighFive::DataSet d_head_qpt = file.getDataSet("/HEAD_QPT");
    // read the k/q-points in yambo internal coordinates
    d_head_qpt.read(yamboQPoints);
    if (yamboQPoints.rows() != 3)
        Error("Developer error: qpoints are transposed in yambo?");
    numPoints = yamboQPoints.cols();
    if (numPoints != numK)
        Error("Electronic grid from input is inconsistent with Yambo.");
    yamboQPoints = yamboQPoints.reshaped(numPoints,3).eval();
    yamboQPoints.transposeInPlace();

    // if we have more MPI processes than points, the code will hang
    // on the allReduceSum below. Plus, it's wasteful.
    if (mpi->getSize() > numPoints) {
      Error("Bloch2Wannier transformation cannot be run with more "
          "MPI processes than input Yambo PAR files.");
    }

    // Convert from yambo internal coordiantes to phoebe's coordinates
    //TODO make sure the definition of alat is always the diagonal of the
    //lattice vectors
    {
      Eigen::Matrix3d latt = crystal.getDirectUnitCell();
      Eigen::Matrix3d rlattInv = crystal.getReciprocalUnitCell().inverse();
      Eigen::Vector3d alat = {latt(0,0),latt(1,1),latt(2,2)};

      for(int iq = 0; iq < numPoints; iq++) {
        Eigen::Vector3d temp;
        for( int i : {0,1,2} )
          temp(i) = (abs(yamboQPoints(i, iq)) < 1e-7) ? 0.0 : yamboQPoints(i, iq) * twoPi / alat(i);
        temp = rlattInv * temp;
        for( int i : {0,1,2} )
          yamboQPoints(i, iq) = temp(i);
        std::cout << iq << " " << temp.transpose() << std::endl;
        }
    }


    Eigen::VectorXi bandExtrema;
    HighFive::DataSet d_bands = file.getDataSet("/Bands");
    d_bands.read(bandExtrema);
    numBands = bandExtrema(1) - bandExtrema(0) + 1; // numBands can be larger than numWannier, i.e. disentanglement
    bandOffset = bandExtrema.minCoeff();
  }
  // push these quantities to all processes
  mpi->bcast(&numPoints);
  mpi->bcast(&numBands);
  mpi->bcast(&bandOffset);
  // head process has already allocated this, other ones need to
  if (!mpi->mpiHead()) {
    yamboQPoints.resize(3, numPoints);
  }
  mpi->bcast(&yamboQPoints);

  // read U matrices
  std::string wannierPrefix = context.getWannier90Prefix();
  Eigen::Tensor<std::complex<double>, 3> uMatrices;
  // uMatrices have size (numBands, numWannier, numKPoints)
  uMatrices = ElPhQeToPhoebeApp::setupRotationMatrices(wannierPrefix, kPoints, false);

  // check that the band numbers from Yambo header file are the same as from Wannier90
  if (numBands != uMatrices.dimension(0)) {
    Error("Band number not aligned between Yambo and Wannier90.");
  }
  int numWannier = uMatrices.dimension(1);

  // generate the Bravais vectors and weights
  auto t3 = crystal.buildWignerSeitzVectors(kMesh);
  Eigen::MatrixXd elBravaisVectors = std::get<0>(t3);
  Eigen::VectorXd elDegeneracies = std::get<1>(t3);
  int numR = elDegeneracies.size();
  // TODO use Wigner-Seitz vectors from *_wsvec.dat

  //-------------------------

  // now I can allocate the mega tensor of 4-point interaction
  Eigen::Tensor<std::complex<double>, 7> qpCoupling(numPoints, numPoints,
                                                    numPoints, numBands,
                                                    numBands, numBands, numBands);
  qpCoupling.setZero();

  // now we read the files
  // TODO if we were smarter, we could figure out which files each process had,
  // and which indices they would write to -- ideally, it would be in a certain block
  // of the matrix, so we could parallelize the transform of these blocks
  // and then write them in parallel
  for (int iQ : mpi->divideWorkIter(numPoints)) {

    std::string yamboPrefix = context.getYamboInteractionPrefix();
    // we have to offset these qs by two, and the files start from Q1,
    // and the first file has only header info
    std::string fileName = yamboPrefix + "BS_PAR_Q" + std::to_string(iQ+1) + ".hdf5";
    // Open the hdf5 file
    HighFive::File file(fileName, HighFive::File::ReadOnly);
    // Set up hdf5 datasets
    HighFive::DataSet d_bse_resonant = file.getDataSet("/BSE_RESONANT");
    // read in the data
    Eigen::MatrixXcd yamboKernel;
    d_bse_resonant.read(yamboKernel);

    std::cout << "Reading in " << fileName << std::endl;

    // dimension of yambo W matrix
    int M = int(sqrt(yamboKernel.size()));
    if (M > numPoints * numBands * numBands) {
      Error("BSE kernel size not consistent with the rest of the input");
    }

    // note that the yamboKernel is only the lower triangle of the
    // matrix, and the upper triangle is nonsense that needs to be filled
    // first pass: we make sure the matrix is complete
    //int count = 0;
    #pragma omp parallel for
    for (int i = 0; i < M; ++i) {
      for (int j = i + 1; j < M; ++j) {// in this way j > i
        yamboKernel(j, i) = std::conj(yamboKernel(i, j));
      }
    }
    yamboKernel.adjointInPlace();

    // read this auxiliary mapping for Yambo indices
    Eigen::MatrixXi ikbz_ib1_ib2_isp2_isp1;
    HighFive::DataSet d_ikbz = file.getDataSet("/IKBZ_IB1_IB2_ISP2_ISP1");
    d_ikbz.read(ikbz_ib1_ib2_isp2_isp1);
    // reshape to correct dimension
    ikbz_ib1_ib2_isp2_isp1 = ikbz_ib1_ib2_isp2_isp1.reshaped(M,5).eval();
    ikbz_ib1_ib2_isp2_isp1.transposeInPlace();

    // Get the exciton Q this file is associated with
    Eigen::Vector3d excitonQ = yamboQPoints.col(iQ);

    #pragma omp parallel for
    for (int iYamboBSE = 0; iYamboBSE < M; ++iYamboBSE) {
      int iYamboIk1 = ikbz_ib1_ib2_isp2_isp1(0, iYamboBSE);
      int iYamboIb1 = ikbz_ib1_ib2_isp2_isp1(1, iYamboBSE);
      int iYamboIb2 = ikbz_ib1_ib2_isp2_isp1(2, iYamboBSE);
      //int iYamboIS2 = ikbz_ib1_ib2_isp2_isp1(3, iYamboBSE);
      //int iYamboIS1 = ikbz_ib1_ib2_isp2_isp1(4, iYamboBSE);
      // k1
      Eigen::Vector3d thisiK = yamboQPoints.col(iYamboIk1 - 1);
      for (int jYamboBSE = 0; jYamboBSE < M; ++jYamboBSE) {
        int jYamboIk1 = ikbz_ib1_ib2_isp2_isp1(0, jYamboBSE);
        int jYamboIb1 = ikbz_ib1_ib2_isp2_isp1(1, jYamboBSE);
        int jYamboIb2 = ikbz_ib1_ib2_isp2_isp1(2, jYamboBSE);
        //int jYamboIS2 = ikbz_ib1_ib2_isp2_isp1(3, jYamboBSE);
        //int jYamboIS1 = ikbz_ib1_ib2_isp2_isp1(4, jYamboBSE);
        // k3
        Eigen::Vector3d thisjK = yamboQPoints.col(jYamboIk1 - 1);
        // because k3 - k4 = Q in Yambo's convention, and then we swap k2 and k4 by defining k2 = k3-Q.
        // k2
        Eigen::Vector3d thisK2 = thisjK - excitonQ;

        int ikk1 = kPoints.getIndex(thisiK, tolerance);
        int ikk2 = kPoints.getIndex(thisK2, tolerance);
        int ikk3 = kPoints.getIndex(thisjK, tolerance);

        // band indexing starts from bottom of offset
        // note that here, we swap again s2 and s4 by switching jYamboIb2 (s4 in yambo convention) -> iYamboIb2 (s2 in yambo convention)
        int ib1 = iYamboIb1 -  bandOffset;
        int ib2 = jYamboIb2 -  bandOffset;
        int ib3 = jYamboIb1 -  bandOffset;
        int ib4 = iYamboIb2 -  bandOffset;

        qpCoupling(ikk1, ikk2, ikk3, ib1, ib2, ib3, ib4) = yamboKernel(iYamboBSE, jYamboBSE) * energyHaToEv;
        Eigen::Vector3d temp1 = {0.25,0.0,0.0};
        Eigen::Vector3d temp2 = {0.25,0.25,0.0};
        if(abs((thisiK-temp1).norm()) < 1e-3 && abs((thisK2-temp2).norm()) < 1e-3 && ib1 == 3 && ib2 == 2 && ib3 == 3 && ib4 == 2) {
          std::cout << thisiK.transpose() << " " << thisK2.transpose() << " " << thisjK.transpose() << " " << qpCoupling(ikk1, ikk2, ikk3, ib1, ib2, ib3, ib4) << std::endl;
        }
        /*if(ikk1 == 1 && ikk2 == 2 && ikk3 == 3 && ib1 == 1 && ib2 == 1 && ib3 == 1 && ib4 == 1) {
          std::cout << "k1 " << thisiK.transpose() << std::endl;
          std::cout << "k2 " << thisK2.transpose() << std::endl;
          std::cout << "k3 " << thisjK.transpose() << std::endl;
          std::cout << "Q " << excitonQ.transpose() << std::endl;
          std::cout << qpCoupling(ikk1, ikk2, ikk3, ib1, ib2, ib3, ib4) << std::endl;
        }*/
      }
    }
  }
  // if we don't use bigAllReduce, this will crash, as it often overflows the int arguments of the MPI call
  // NOTE: this is somewhat slow, and also I wish we could reduce after FTing
  mpi->bigAllReduceSum(&qpCoupling);
  mpi->barrier();

  if (mpi->mpiHead()) {
    std::cout << "Starting the Wannier transform.\n";
  }

  // TODO -- for electrons, this procedure is probably good
  // however, for holes, Changpeng wonders if we need to take the complex conjugate of
  // the U matrix (to essentially swap bra and ket)

  // NOTE: this transform+rotation series corresponds to eq 60 in Changpeng's notes
  // as of 3/22/23

  // first we do the rotation on k4, the implicit index
  Eigen::Tensor<std::complex<double>, 7> Gamma1(numK, numK, numK, numBands, numBands, numBands, numWannier);
  {
    for (int ik1 = 0; ik1 < numK; ++ik1) {
      Eigen::Vector3d k1C = kPoints.getPointCoordinates(ik1, Points::crystalCoordinates);
      for (int ik2 = 0; ik2 < numK; ++ik2) {
        Eigen::Vector3d k2C = kPoints.getPointCoordinates(ik2, Points::crystalCoordinates);
        for (int ik3 = 0; ik3 < numK; ++ik3) {
          Eigen::Vector3d k3C = kPoints.getPointCoordinates(ik3, Points::crystalCoordinates);
          // K1-K2 = K3-K4 <-- original momentum conservation relation in Yambo
          // we worked to alter this convention above, so now we just have the typical
          // momentum conservation condition K1 + K2 = K3 + K4
          Eigen::Vector3d k4C = k1C + k2C - k3C;
          int ik4 = kPoints.getIndex(k4C);

          #pragma omp parallel for collapse(4)
          for (int ib1 = 0; ib1 < numBands; ++ib1) {
            for (int ib2 = 0; ib2 < numBands; ++ib2) {
              for (int ib3 = 0; ib3 < numBands; ++ib3) {
                for (int iw4 = 0; iw4 < numWannier; ++iw4) {
                  std::complex<double> tmp = {0., 0.};    // to avoid overwriting in OpenMP parallelism, tmp should be defined with the loop
                  for (int ib4 = 0; ib4 < numBands; ++ib4) {
                    tmp += qpCoupling(ik1, ik2, ik3, ib1, ib2, ib3, ib4)
                        * uMatrices(ib4, iw4, ik4);
                  }
                  Gamma1(ik1, ik2, ik3, ib1, ib2, ib3, iw4) = tmp;
                }
              }
            }
          }
        }
      }
    }
  }
  qpCoupling.resize(0, 0, 0, 0, 0, 0, 0);

  // next, rotation on ik3
  Eigen::Tensor<std::complex<double>, 7> Gamma2(numK, numK, numK, numBands, numBands, numWannier, numWannier);
  {
    #pragma omp parallel for collapse(6)
    for (int ik1 = 0; ik1 < numK; ++ik1) {
      for (int ik2 = 0; ik2 < numK; ++ik2) {
        for (int ik3 = 0; ik3 < numK; ++ik3) {
          for (int ib1 = 0; ib1 < numBands; ++ib1) {
            for (int ib2 = 0; ib2 < numBands; ++ib2) {
              for (int iw4 = 0; iw4 < numWannier; ++iw4) {

                for (int iw3 = 0; iw3 < numWannier; ++iw3) {
                  std::complex<double> tmp = {0., 0.};
                  for (int ib3 = 0; ib3 < numBands; ++ib3) {
                    tmp += Gamma1(ik1, ik2, ik3, ib1, ib2, ib3, iw4)
                        * uMatrices(ib3, iw3, ik3);
                  }
                  Gamma2(ik1, ik2, ik3, ib1, ib2, iw3, iw4) = tmp;
                }
              }
            }
          }
        }
      }
    }
  }
  Gamma1.resize(0, 0, 0, 0, 0, 0, 0);

  // rotation on ik2
  Eigen::Tensor<std::complex<double>, 7> Gamma3(numK, numK, numK, numBands, numWannier, numWannier, numWannier);
  {
    #pragma omp parallel for collapse(6)
    for (int ik1 = 0; ik1 < numK; ++ik1) {
      for (int ik2 = 0; ik2 < numK; ++ik2) {
        for (int ik3 = 0; ik3 < numK; ++ik3) {
          for (int ib1 = 0; ib1 < numBands; ++ib1) {
            for (int iw3 = 0; iw3 < numWannier; ++iw3) {
              for (int iw4 = 0; iw4 < numWannier; ++iw4) {

                for (int iw2 = 0; iw2 < numWannier; ++iw2) {
                  std::complex<double> tmp = {0., 0.};
                  for (int ib2 = 0; ib2 < numBands; ++ib2) {
                    tmp += Gamma2(ik1, ik2, ik3, ib1, ib2, iw3, iw4)
                        * std::conj(uMatrices(ib2, iw2, ik2));
                  }
                  Gamma3(ik1, ik2, ik3, ib1, iw2, iw3, iw4) = tmp;
                }
              }
            }
          }
        }
      }
    }
  }
  Gamma2.resize(0, 0, 0, 0, 0, 0, 0);

  // rotation on ik1
  Eigen::Tensor<std::complex<double>, 7> Gamma4(numK, numK, numK, numWannier, numWannier, numWannier, numWannier);
  {
    #pragma omp parallel for collapse(6)
    for (int ik1 = 0; ik1 < numK; ++ik1) {
      for (int ik2 = 0; ik2 < numK; ++ik2) {
        for (int ik3 = 0; ik3 < numK; ++ik3) {
          for (int iw2 = 0; iw2 < numWannier; ++iw2) {
            for (int iw3 = 0; iw3 < numWannier; ++iw3) {
              for (int iw4 = 0; iw4 < numWannier; ++iw4) {
                for (int iw1 = 0; iw1 < numWannier; ++iw1) {
                  std::complex<double> tmp = {0., 0.};
                  for (int ib1 = 0; ib1 < numBands; ++ib1) {
                     tmp += Gamma3(ik1, ik2, ik3, ib1, iw2, iw3, iw4)
                        * std::conj(uMatrices(ib1, iw1, ik1));
                  }
                  Gamma4(ik1, ik2, ik3, iw1, iw2, iw3, iw4) = tmp;
                  //Gamma4(ik1, ik2, ik3, iw1, iw2, iw3, iw4) = qpCoupling(ik1, ik2, ik3, iw1+2, iw2+2, iw3+2, iw4+2);
                }
              }
            }
          }
        }
      }
    }
  }
  Gamma3.resize(0, 0, 0, 0, 0, 0, 0);

  // let's prepare now for the Fourier transform
  // note that the phase factors are the same for all wavevector indices
  Eigen::MatrixXcd phases(numK, numR);
  {
    #pragma omp parallel for collapse(2)
    for (int ik = 0; ik < numK; ++ik) {
      for (int iR = 0; iR < numR; ++iR) {
        Eigen::Vector3d k = kPoints.getPointCoordinates(ik, Points::cartesianCoordinates);
        Eigen::Vector3d R = elBravaisVectors.col(iR);
        double arg = k.dot(R);
        phases(ik, iR) = exp(complexI * arg) / double(numK);
      }
    }
  }

  //if (mpi->mpiHead()) {
  //    for (int iR = 0; iR < numR; ++iR) {
  //        std::cout << elBravaisVectors.col(iR) << std::endl;}
  //}

  // Fourier transform on ik1
  Eigen::Tensor<std::complex<double>, 7> Gamma5(numR, numK, numK, numWannier, numWannier, numWannier, numWannier);
  {
    #pragma omp parallel for collapse(6)
    for (int ik2 = 0; ik2 < numK; ++ik2) {
      for (int ik3 = 0; ik3 < numK; ++ik3) {
        for (int iw1 = 0; iw1 < numWannier; ++iw1) {
          for (int iw2 = 0; iw2 < numWannier; ++iw2) {
            for (int iw3 = 0; iw3 < numWannier; ++iw3) {
              for (int iw4 = 0; iw4 < numWannier; ++iw4) {

                for (int iR1 = 0; iR1 < numR; ++iR1) {
                  std::complex<double> tmp = {0., 0.};
                  for (int ik1 = 0; ik1 < numK; ++ik1) {
                    tmp += phases(ik1, iR1) * Gamma4(ik1, ik2, ik3, iw1, iw2, iw3, iw4);
                  }
                  Gamma5(iR1, ik2, ik3, iw1, iw2, iw3, iw4) = tmp;
                }
              }
            }
          }
        }
      }
    }
    Gamma4.resize(0, 0, 0, 0, 0, 0, 0);
  }

  // Fourier transform on ik2
  Eigen::Tensor<std::complex<double>, 7> Gamma6(numR, numR, numK, numWannier, numWannier, numWannier, numWannier);
  {
    #pragma omp parallel for collapse(6)
    for (int iR1 = 0; iR1 < numR; ++iR1) {
      for (int ik3 = 0; ik3 < numK; ++ik3) {
        for (int iw1 = 0; iw1 < numWannier; ++iw1) {
          for (int iw2 = 0; iw2 < numWannier; ++iw2) {
            for (int iw3 = 0; iw3 < numWannier; ++iw3) {
              for (int iw4 = 0; iw4 < numWannier; ++iw4) {

                for (int iR2 = 0; iR2 < numR; ++iR2) {
                  std::complex<double> tmp = {0., 0.};
                  for (int ik2 = 0; ik2 < numK; ++ik2) {
                    tmp += phases(ik2, iR2) * Gamma5(iR1, ik2, ik3, iw1, iw2, iw3, iw4);
                  }
                  Gamma6(iR1, iR2, ik3, iw1, iw2, iw3, iw4) = tmp;
                }
              }
            }
          }
        }
      }
    }
    Gamma5.resize(0, 0, 0, 0, 0, 0, 0);
  }

  // Fourier transform on ik3
  Eigen::Tensor<std::complex<double>, 7> GammaW(numR, numR, numR, numWannier, numWannier, numWannier, numWannier);
  {
    #pragma omp parallel for collapse(6)
    for (int iR1 = 0; iR1 < numR; ++iR1) {
      for (int iR2 = 0; iR2 < numR; ++iR2) {
        for (int iw1 = 0; iw1 < numWannier; ++iw1) {
          for (int iw2 = 0; iw2 < numWannier; ++iw2) {
            for (int iw3 = 0; iw3 < numWannier; ++iw3) {
              for (int iw4 = 0; iw4 < numWannier; ++iw4) {

                for (int iR3 = 0; iR3 < numR; ++iR3) {
                  std::complex<double> tmp = {0., 0.};
                  for (int ik3 = 0; ik3 < numK; ++ik3) {
                    tmp += std::conj(phases(ik3, iR3)) * Gamma6(iR1, iR2, ik3, iw1, iw2, iw3, iw4);
                  }
                  GammaW(iR1, iR2, iR3, iw1, iw2, iw3, iw4) = tmp;
                }
              }
            }
          }
        }
      }
    }
    Gamma6.resize(0, 0, 0, 0, 0, 0, 0);
  }
  mpi->barrier();

  // Now, let's write it to file
  writeWannierCoupling(context, GammaW, elDegeneracies, elBravaisVectors, kMesh);
}

void ElElToPhoebeApp::writeWannierCoupling(
    Context &context,
    Eigen::Tensor<std::complex<double>, 7> &gWannier,
    const Eigen::VectorXd &elDegeneracies,
    const Eigen::MatrixXd &elBravaisVectors,
    const Eigen::Vector3i &kMesh) {

#ifdef HDF5_AVAIL
  if (mpi->getSize() < 4) {
    // Note: this HDF5 had already been reported and being worked on.
    // It's beyond the purpose of Phoebe's project.
    Warning("HDF5 with <4 MPI process may crash (due to a "
            "library bug),\nuse more MPI processes if that happens");
  }

  int numR = gWannier.dimension(0);
  int numWannier = gWannier.dimension(6);

  std::string phoebePrefixQE = context.getQuantumEspressoPrefix();
  std::string outFileName = phoebePrefixQE + ".phoebe.elel.hdf5";
  // if the hdf5 file is there already, we want to delete it. Occasionally
  // these files seem to get stuck open when a process dies while writing to
  // them, (even if a python script dies) and then they can't be overwritten
  // properly.
  std::remove(&outFileName[0]);
  mpi->barrier();

  // TODO may be a problem here if run with more than 1 MPI process, Jenny should check this
  try {
  {

    mpi->barrier();// wait for file to be overwritten

    // now open the file in serial mode
    // because we do the parallelization by hand
    HighFive::File file(outFileName, HighFive::File::Overwrite);

    // create buffer to save a slice of the tensor, at fixed ir1
    Eigen::Tensor<std::complex<double>, 6> slice;
    slice.resize(numR, numR, numWannier, numWannier, numWannier, numWannier);

    // All processes have the same gWannier tensor, but each one
    // will write a few R1s to file
    auto iR1Iterator = mpi->divideWorkIter(numR);

    for (int iR1 : iR1Iterator) {

      // select a slice
      for (int iR2 = 0; iR2 < numR; ++iR2) {
        for (int iR3 = 0; iR3 < numR; ++iR3) {
          for (int iw1 = 0; iw1 < numWannier; iw1++) {
            for (int iw2 = 0; iw2 < numWannier; iw2++) {
              for (int iw3 = 0; iw3 < numWannier; iw3++) {
                for (int iw4 = 0; iw4 < numWannier; iw4++) {
                  slice(iR2, iR3, iw1, iw2, iw3, iw4) = gWannier(iR1, iR2, iR3, iw1, iw2, iw3, iw4);
                  //if (iR1 == 0 && iR2 <= 1 && iR3 <= 1) {
                  //    std::cout << iw1 << " " << iw2 << " " << iw3 << " " << iw4 << " " << gWannier(iR1, iR2, iR3, iw1, iw2, iw3, iw4) << std::endl;}
                }
              }
            }
          }
        }
      }

      // flatten the tensor in a vector
      Eigen::VectorXcd flatSlice =
          Eigen::Map<Eigen::VectorXcd, Eigen::Unaligned>(slice.data(), slice.size());
      size_t flatSize = std::pow(numR, 2) * std::pow(numWannier,4);

      auto maxSize = int(pow(1000, 3)) / sizeof(std::complex<double>);
      if(flatSize > maxSize) {
        if(mpi->mpiHead()) {
          Warning("You tried to write a slice of data to HDF5 which is greater than the allowed\n"
                "maximum HDF5 write size of ~2GB. If this fails, please contact the developers to let them know\n"
                "you need bunched write functionality for this app.");
          }
      }

      std::string datasetName = "/gWannier_" + std::to_string(iR1);

      // create dataset
      HighFive::DataSet dslice = file.createDataSet<std::complex<double>>(
          datasetName, HighFive::DataSpace::From(flatSlice));

      // write to hdf5
      dslice.write(flatSlice);
    }
    mpi->barrier();
  } // closing braces
  } catch (std::exception &error) {
    Error("Issue writing el-l Wannier representation matrix elements to hdf5.");
  }
  mpi->barrier();

  // now we write some extra info
  if (mpi->mpiHead()) {
    try {

      HighFive::File file(outFileName, HighFive::File::ReadWrite);

      HighFive::DataSet dnWannier = file.createDataSet<int>(
          "/numWannier", HighFive::DataSpace::From(numWannier));
      dnWannier.write(numWannier);

      HighFive::DataSet dkMesh =
          file.createDataSet<int>("/kMesh", HighFive::DataSpace::From(kMesh));
      dkMesh.write(kMesh);

      HighFive::DataSet dElBravais = file.createDataSet<double>(
          "/elBravaisVectors", HighFive::DataSpace::From(elBravaisVectors));
      dElBravais.write(elBravaisVectors);

      HighFive::DataSet dElDegeneracies = file.createDataSet<double>(
          "/elDegeneracies", HighFive::DataSpace::From(elDegeneracies));
      dElDegeneracies.write(elDegeneracies);

    } catch (std::exception &error) {
      Error("Issue writing elel Wannier header to hdf5.");
    }
  }
  mpi->barrier();
#else
  Error("Phoebe requires HDF5 support to write el-el interaction to file.\n"
        "Please recompile with HDF5.");
#endif
}

void ElElToPhoebeApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");
  throwErrorIfUnset(context.getQuantumEspressoPrefix(), "QuantumEspressoPrefix");
  throwErrorIfUnset(context.getWannier90Prefix(), "Wannier90Prefix");
  throwErrorIfUnset(context.getYamboInteractionPrefix(), "YamboInteractionPrefix");
  throwErrorIfUnset(context.getKMesh(), "kMesh");
  // check that crystal structure was provided
  std::string crystalMsg = "crystal structure";
  throwErrorIfUnset(context.getInputAtomicPositions(), crystalMsg);
  throwErrorIfUnset(context.getInputSpeciesNames(), crystalMsg);
  throwErrorIfUnset(context.getInputAtomicSpecies(), crystalMsg);

}

