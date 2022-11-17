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
//#include "H5Cpp.h"
#endif

void ElElToPhoebeApp::run(Context &context) {

  #ifndef HDF5_AVAIL
  Error("In order to run the el-el to Phoebe app, you must build Phoebe with HDF5.");
  #endif

  auto t1 = Parser::parseElHarmonicWannier(context);
  auto crystal = std::get<0>(t1);
  auto electronH0 = std::get<1>(t1);

  // now we parse the first header
  Eigen::MatrixXd yamboQPoints;
  int numBands, bandOffset, numPoints;
  if (mpi->mpiHead()) {

    std::string yamboPrefix = context.getYamboInteractionPrefix();
    // Q1 should always exist
    std::string fileName = yamboPrefix + "BS_head_Q1.hdf5";
    // Open the hdf5 file
    HighFive::File file(fileName, HighFive::File::ReadOnly);
    // Set up hdf5 datasets
    HighFive::DataSet d_head_qpt = file.getDataSet("/HEAD_QPT");
    // read the k/q-points in crystal coordinates
    Eigen::MatrixXd yamboQPoints;
    d_head_qpt.read(yamboQPoints);
    int numDim = yamboQPoints.rows();
    if (numDim != 3) Error("Developer error: qpoints are transposed in yambo?");
    numPoints = yamboQPoints.cols();

    Eigen::MatrixXi bandExtrema;
    HighFive::DataSet d_bands = file.getDataSet("/Bands");
    d_bands.read(bandExtrema);
    std::cout << "band extrema " << bandExtrema << std::endl;
    numBands = bandExtrema(1) - bandExtrema(0) + 1; // if band range is 4,5, that's two bands
    bandOffset = bandExtrema.minCoeff();
    std::cout << electronH0.getNumBands() << " band extrema " << bandExtrema << " offset " << bandOffset << " numBands " << numBands <<  std::endl;
  }

  mpi->bcast(&numPoints);
  mpi->bcast(&numBands);
  mpi->bcast(&bandOffset);
  if (!mpi->mpiHead()) {
    yamboQPoints.resize(3, numPoints);
  }
  mpi->bcast(&yamboQPoints);

  //if (electronH0.getNumBands() != numBands) {
  //  Error("Yambo and Wannier have run with different band number");
  //}

  //----------------
  // set up k-points
  Eigen::Vector3i kMesh = {8, 8, 1}; // TODO: get this from somewhere
  Points kPoints(crystal, kMesh);
  int numK = kMesh.prod();

  // read U matrices
  std::string wannierPrefix = context.getWannier90Prefix();
  Eigen::Tensor<std::complex<double>, 3> uMatrices;
  // uMatrices have size (numBands, numWannier, numKPoints)
  uMatrices = ElPhQeToPhoebeApp::setupRotationMatrices(wannierPrefix, kPoints, false);

  if (numBands != uMatrices.dimension(0)) {
    Error("Band number not aligned between Yambo and Wannier90.");
  }
  int numWannier = uMatrices.dimension(1);

  // generate the Bravais vectors
  auto t3 = crystal.buildWignerSeitzVectors(kMesh);
  Eigen::MatrixXd elBravaisVectors = std::get<0>(t3);
  Eigen::VectorXd elDegeneracies = std::get<1>(t3);
  int numR = elDegeneracies.size();

  //-------------------------

  // now I can allocate the mega tensor of 4-point interaction
  Eigen::Tensor<std::complex<double>, 7> qpCoupling(numPoints, numPoints,
                                                    numPoints, numBands,
                                                    numBands, numBands, numBands);
  qpCoupling.setZero();

  // now we read the files
  for (int iQ : mpi->divideWorkIter(numPoints)) {

    std::string yamboPrefix = context.getYamboInteractionPrefix();
    // tentative reading of the BSE kernel
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

    std::cout << "read in " << fileName << std::endl;

    //int M = sqrt(yamboKernel_.size() / 2);
    int M = sqrt(yamboKernel.size());

    std::cout << "M " << M << " numPoints " << numPoints << " " << numBands << " " << std::endl;

    if (M != numPoints * numBands * numBands) {
      Error("BSE kernel size not consistent with the rest of the input");
    }

    // note che in yamboKernel, solo la lower triangle is filled
    // i.e. kernel[4,1] ha un numero ragionevole, ma k[1,4] no

    // first pass: we make sure the matrix is complete
    #pragma omp parallel for
    for (int i = 0; i < M; ++i) {
      for (int j = i + 1; j < M; ++j) {// in this way j > i
        // TODO: check if this has to be a complex conjugate
        yamboKernel(i, j) = yamboKernel(j, i);
      }
    }

    // read this auxiliary mapping for Yambo indices
    Eigen::MatrixXi ikbz_ib1_ib2_isp2_isp1;
    HighFive::DataSet d_ikbz = file.getDataSet("/IKBZ_IB1_IB2_ISP2_ISP1");
    d_ikbz.read(ikbz_ib1_ib2_isp2_isp1);

    // now we substitute this back into the big coupling tensor
    Eigen::Vector3d excitonQ = yamboQPoints.col(iQ);

  //  #pragma omp parallel for
    for (int iYamboBSE = 0; iYamboBSE < M; ++iYamboBSE) {
      int iYamboIk1 = ikbz_ib1_ib2_isp2_isp1(0, iYamboBSE);
      int iYamboIb1 = ikbz_ib1_ib2_isp2_isp1(1, iYamboBSE);
      int iYamboIb2 = ikbz_ib1_ib2_isp2_isp1(2, iYamboBSE);
      int iYamboIS2 = ikbz_ib1_ib2_isp2_isp1(3, iYamboBSE);
      int iYamboIS1 = ikbz_ib1_ib2_isp2_isp1(4, iYamboBSE);
      std::cout << "iYamboIk1 " << iYamboIk1  << std::endl;
      Eigen::Vector3d thisiK = yamboQPoints.col(iYamboIk1 - 1);
      for (int jYamboBSE = 0; jYamboBSE < M; ++jYamboBSE) {
        int jYamboIk1 = ikbz_ib1_ib2_isp2_isp1(0, jYamboBSE);
        int jYamboIb1 = ikbz_ib1_ib2_isp2_isp1(1, jYamboBSE);
        int jYamboIb2 = ikbz_ib1_ib2_isp2_isp1(2, jYamboBSE);
        int jYamboIS2 = ikbz_ib1_ib2_isp2_isp1(3, jYamboBSE);
        int jYamboIS1 = ikbz_ib1_ib2_isp2_isp1(4, jYamboBSE);
        std::cout << "jYamboIk1 " << jYamboIk1  << std::endl;
        Eigen::Vector3d thisjK = yamboQPoints.col(jYamboIk1 - 1);
        Eigen::Vector3d thisK2 = thisiK - excitonQ;

        int ikk1 = kPoints.getIndex(thisiK);
        int ikk2 = kPoints.getIndex(thisK2);
        int ikk3 = kPoints.getIndex(thisjK);
        int ib1 = iYamboIb1 - bandOffset;// because iYambo starts from 4
        int ib2 = iYamboIb2 - bandOffset;
        int ib3 = jYamboIb1 - bandOffset;
        int ib4 = jYamboIb2 - bandOffset;

        std::complex<double> z = {0,0}; //yamboKernel(iYamboBSE, jYamboBSE); //{yamboKernel(iYamboBSE, jYamboBSE, 0),
                                 // yamboKernel(iYamboBSE, jYamboBSE, 1)};
        //qpCoupling(ikk1, ikk2, ikk3, ib1, ib2, ib3, ib4) = z;
      }
    }
  }
  mpi->bcast(&qpCoupling);

  if (mpi->mpiHead()) {
    std::cout << "Starting the Wannier transform.\n";
  }

  if(mpi->mpiHead()) std::cout << "points " << numK << " " << numPoints << std::endl;

  // first we do the rotation on k4, the implicit index
  Eigen::Tensor<std::complex<double>, 7> Gamma1(numK, numK, numK, numBands, numBands, numBands, numWannier);
  {
    std::complex<double> tmp;
    for (int ik1 = 0; ik1 < numK; ++ik1) {
      Eigen::Vector3d k1C = kPoints.getPointCoordinates(ik1, Points::crystalCoordinates);
      for (int ik2 = 0; ik2 < numK; ++ik2) {
        Eigen::Vector3d k2C = kPoints.getPointCoordinates(ik2, Points::crystalCoordinates);
        for (int ik3 = 0; ik3 < numK; ++ik3) {
          Eigen::Vector3d k3C = kPoints.getPointCoordinates(ik3, Points::crystalCoordinates);
          Eigen::Vector3d k4C = k1C + k2C - k3C;
          int ik4 = kPoints.getIndex(k4C);

          #pragma omp parallel for collapse(4)
          for (int ib1 = 0; ib1 < numBands; ++ib1) {
            for (int ib2 = 0; ib2 < numBands; ++ib2) {
              for (int ib3 = 0; ib3 < numBands; ++ib3) {
                for (int iw4 = 0; iw4 < numWannier; ++iw4) {
                  tmp = {0., 0.};
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
    std::complex<double> tmp;
#pragma omp parallel for collapse(6)
    for (int ik1 = 0; ik1 < numK; ++ik1) {
      for (int ik2 = 0; ik2 < numK; ++ik2) {
        for (int ik3 = 0; ik3 < numK; ++ik3) {
          for (int ib1 = 0; ib1 < numBands; ++ib1) {
            for (int ib2 = 0; ib2 < numBands; ++ib2) {
              for (int iw4 = 0; iw4 < numWannier; ++iw4) {

                for (int iw3 = 0; iw3 < numWannier; ++iw3) {
                  tmp = {0., 0.};
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
    std::complex<double> tmp;
#pragma omp parallel for collapse(6)
    for (int ik1 = 0; ik1 < numK; ++ik1) {
      for (int ik2 = 0; ik2 < numK; ++ik2) {
        for (int ik3 = 0; ik3 < numK; ++ik3) {
          for (int ib1 = 0; ib1 < numBands; ++ib1) {
            for (int iw3 = 0; iw3 < numWannier; ++iw3) {
              for (int iw4 = 0; iw4 < numWannier; ++iw4) {

                for (int iw2 = 0; iw2 < numWannier; ++iw2) {
                  tmp = {0., 0.};
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
    std::complex<double> tmp;
#pragma omp parallel for collapse(6)
    for (int ik1 = 0; ik1 < numK; ++ik1) {
      for (int ik2 = 0; ik2 < numK; ++ik2) {
        for (int ik3 = 0; ik3 < numK; ++ik3) {
          for (int iw2 = 0; iw2 < numWannier; ++iw2) {
            for (int iw3 = 0; iw3 < numWannier; ++iw3) {
              for (int iw4 = 0; iw4 < numWannier; ++iw4) {
                for (int iw1 = 0; iw1 < numWannier; ++iw1) {
                  tmp = {0., 0.};
                  for (int ib1 = 0; ib1 < numBands; ++ib1) {
                     tmp += Gamma3(ik1, ik2, ik3, ib1, iw2, iw3, iw4)
                        * std::conj(uMatrices(ib1, iw1, ik1));
                  }
                  Gamma4(ik1, ik2, ik3, iw1, iw2, iw3, iw4) = tmp;
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
        Eigen::Vector3d R = elBravaisVectors.row(iR);
        double arg = k.dot(R);
        phases(ik, iR) = exp(complexI * arg);
      }
    }
  }

  // Fourier transform on ik1
  Eigen::Tensor<std::complex<double>, 7> Gamma5(numR, numK, numK, numWannier, numWannier, numWannier, numWannier);
  {
    std::complex<double> tmp;
    #pragma omp parallel for collapse(6)
    for (int ik2 = 0; ik2 < numK; ++ik2) {
      for (int ik3 = 0; ik3 < numK; ++ik3) {
        for (int iw1 = 0; iw1 < numWannier; ++iw1) {
          for (int iw2 = 0; iw2 < numWannier; ++iw2) {
            for (int iw3 = 0; iw3 < numWannier; ++iw3) {
              for (int iw4 = 0; iw4 < numWannier; ++iw4) {

                for (int iR = 0; iR < numR; ++iR) {
                  tmp = {0., 0.};
                  for (int ik1 = 0; ik1 < numK; ++ik1) {
                    tmp += phases(ik1, iR) * Gamma4(ik1, ik2, ik3, iw1, iw2, iw3, iw4);
                  }
                  Gamma5(iR, ik2, ik3, iw1, iw2, iw3, iw4) = tmp;
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
    std::complex<double> tmp;
    #pragma omp parallel for collapse(6)
    for (int iR1 = 0; iR1 < numR; ++iR1) {
      for (int ik3 = 0; ik3 < numK; ++ik3) {
        for (int iw1 = 0; iw1 < numWannier; ++iw1) {
          for (int iw2 = 0; iw2 < numWannier; ++iw2) {
            for (int iw3 = 0; iw3 < numWannier; ++iw3) {
              for (int iw4 = 0; iw4 < numWannier; ++iw4) {

                for (int iR2 = 0; iR2 < numR; ++iR2) {
                  tmp = {0., 0.};
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
    std::complex<double> tmp;
    #pragma omp parallel for collapse(6)
    for (int iR1 = 0; iR1 < numR; ++iR1) {
      for (int iR2 = 0; iR2 < numR; ++iR2) {
        for (int iw1 = 0; iw1 < numWannier; ++iw1) {
          for (int iw2 = 0; iw2 < numWannier; ++iw2) {
            for (int iw3 = 0; iw3 < numWannier; ++iw3) {
              for (int iw4 = 0; iw4 < numWannier; ++iw4) {

                for (int iR3 = 0; iR3 < numR; ++iR3) {
                  tmp = {0., 0.};
                  for (int ik3 = 0; ik3 < numK; ++ik3) {
                    tmp += std::conj(phases(ik3, iR3))
                        * Gamma6(iR1, iR2, ik3, iw1, iw2, iw3, iw4);
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

  // Now, let's write it to file
  writeWannierCoupling(context, GammaW, elDegeneracies,
                       elBravaisVectors, kMesh);
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

  try {
    // open the hdf5 file and remove existing files
    if (mpi->mpiHead()) {
      HighFive::File file(outFileName, HighFive::File::Overwrite);
    }
    mpi->barrier();// wait for file to be overwritten

    // now open the file in serial mode
    // because we do the parallelization by hand
    HighFive::File file(outFileName, HighFive::File::Overwrite);

    // create buffer to save a slice of the tensor, at fixed ir1
    Eigen::Tensor<std::complex<double>, 6> slice;
    slice.resize(numR, numR, numWannier, numWannier, numWannier, numWannier);

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
                }
              }
            }
          }
        }
      }

      // flatten the tensor in a vector
      Eigen::VectorXcd flatSlice =
          Eigen::Map<Eigen::VectorXcd, Eigen::Unaligned>(
              slice.data(), slice.size());

      std::string datasetName = "/gWannier_" + std::to_string(iR1);

      // create dataset
      HighFive::DataSet dslice = file.createDataSet<std::complex<double>>(
          datasetName, HighFive::DataSpace::From(flatSlice));

      // write to hdf5
      dslice.write(flatSlice);
    }

  } catch (std::exception &error) {
    Error("Issue writing elel Wannier representation to hdf5.");
  }

  // now we write some extra info
  if (mpi->mpiHead()) {
    try {
      HighFive::File file(outFileName, HighFive::File::ReadWrite);

      //      HighFive::DataSet dnElectrons = file.createDataSet<int>(
      //          "/numElectrons", HighFive::DataSpace::From(numElectrons));
      //      dnElectrons.write(numElectrons);// # of occupied wannier functions

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

#else
  Error("Didn't implement writing el-el interaction to file without HDF5");
#endif
}

void ElElToPhoebeApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");
  throwErrorIfUnset(context.getQuantumEspressoPrefix(), "QuantumEspressoPrefix");
  throwErrorIfUnset(context.getWannier90Prefix(), "Wannier90Prefix");
  throwErrorIfUnset(context.getYamboInteractionPrefix(), "YamboInteractionPrefix");
  // check that crystal structure was provided
  std::string crystalMsg = "crystal structure";
  throwErrorIfUnset(context.getInputAtomicPositions(), crystalMsg);
  throwErrorIfUnset(context.getInputSpeciesNames(), crystalMsg);
  throwErrorIfUnset(context.getInputAtomicSpecies(), crystalMsg);

}

