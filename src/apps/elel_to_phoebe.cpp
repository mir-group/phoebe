#include "elel_to_phoebe.h"
#include "bandstructure.h"
#include "eigen.h"
#include "io.h"
#include "qe_input_parser.h"
#include <iomanip>
#include <sstream>
#include <string>
#include "parser.h"

#ifdef HDF5_AVAIL
#include <highfive/H5Easy.hpp>
#endif

void ElElToPhoebeApp::run(Context &context) {
  (void) context;

  auto t1 = Parser::parseElHarmonicWannier(context);
  auto crystal = std::get<0>(t1);
  auto electronH0 = std::get<1>(t1);

  // tentative reading of kpoints
  {
    std::string fileName = "/home/cepe/Desktop/Changpeng/QP_BSE/ndb.BS_head_Q1";
    // Open the hdf5 file
    HighFive::File file(fileName, HighFive::File::ReadOnly);
    // Set up hdf5 datasets
    HighFive::DataSet d_head_qpt = file.getDataSet("/HEAD_QPT");
    // read the k/q-points in crystal coordinates
    Eigen::MatrixXd yamboQPoints;
    d_head_qpt.read(yamboQPoints);
    int numDim = yamboQPoints.rows();
    int numQ = yamboQPoints.cols();
    if (numDim != 3) Error("qpoints are transposed in yambo?");
    // std::cout << yamboQPoints.transpose() << "\n";
  }

  {
    // tentative reading of the BSE kernel
    std::string fileName = "/home/cepe/Desktop/Changpeng/QP_BSE/ndb.BS_PAR_Q1";
    // Open the hdf5 file
    HighFive::File file(fileName, HighFive::File::ReadOnly);
    // Set up hdf5 datasets
    HighFive::DataSet d_bse_resonant = file.getDataSet("/BSE_RESONANT");
    // read in the data

    std::vector<std::vector<std::vector<double>>> yamboKernel_;
    d_bse_resonant.read(yamboKernel_);
    int M = yamboKernel_.size();
    int N = yamboKernel_[0].size();
    if (N != M) Error("qpoints are transposed in yambo?");

    // note: yamboKernel is a tensor (M,M,2)
    // with the last index being real and imaginary part
    // and M combining k-points and bands

    // now I need to rearrange them
    // NOTE! the

  }

  TOMORROW:
    put this in a loop over the exciton wavevector and have a large tensor;
    next, parse bands from the Wannier file

    note che in yamboKernel_, solo la lower triangle is filled
    i.e. kernel[4,1] ha un numero ragionevole, ma k[1,4] no


  std::cout << "OK!\n";

  return;

//  int numStates = sqrt(BS_K_DIM);
//  Eigen::Tensor<std::complex<double>, 4> kernel(numStates, numStates, numStates, numBands);
//  for (int ik1=0; ik1<numPoints; ++ik1) {
//    for (int ik2=0; ik2<numPoints; ++ik2) {
//      for (int ik3 = 0; ik3 < numPoints; ++ik3) {
//        int ik4 = ...; // TBD
//        for (int ib1=0; ib1<numBands; ++ib1) {
//          for (int ib2 = 0; ib2 < numBands; ++ib2) {
//            for (int ib3 = 0; ib3 < numBands; ++ib3) {
//              for (int ib4 = 0; ib4 < numBands; ++ib4) {
//
//                int bse1 = ...;
//                int bse2 = ...;
//                int is1 = ...;
//                int is2 = ...;
//                int is3 = ...;
//
//                kernel(is1, is2, is3, ib4) = kernelIn(bse1, bse2);
//              }
//            }
//          }
//        }
//      }
//    }
//  }

  // these variables need to be fixed
  //  std::string phoebePrefixQE = context.getQuantumEspressoPrefix();
  //  auto t0 = readQEPhoebeHeader(crystal, phoebePrefixQE);

  Eigen::Vector3i kMesh; // = std::get<1>(t0);
  Eigen::MatrixXd kGridFull; // = std::get<2>(t0);

  Eigen::MatrixXd energies; // = std::get<4>(t0);

  int numQEBands; // = std::get<6>(t0);
  int numElectrons; // = std::get<7>(t0);

  int numBands = 0;
  int numWannier = 0;

  auto t2 = electronH0.getVectors();
  Eigen::MatrixXd bravaisVectors = std::get<0>(t2);
  Eigen::VectorXd bravaisDegeneracies = std::get<1>(t2);
  int numR = bravaisDegeneracies.size();

  //---------

  Points kPoints(crystal, kMesh);

  int numK = kMesh.prod();
  Eigen::Tensor<std::complex<double>, 3> U(numK, numBands, numWannier);
  Eigen::Tensor<std::complex<double>, 7> Gamma(numK, numK, numK, numBands, numBands, numBands, numBands);

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
                    tmp += Gamma(ik1, ik2, ik3, ib1, ib2, ib3, ib4)
                        * U(ik4, ib4, iw4);
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
  Gamma.resize(0, 0, 0, 0, 0, 0, 0);

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
                        * U(ik3, ib3, iw3);
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
                        * std::conj(U(ik2, ib2, iw2));
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
                    tmp += Gamma2(ik1, ik2, ik3, ib1, iw2, iw3, iw4)
                        * std::conj(U(ik1, ib1, iw1));
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
        Eigen::Vector3d R = bravaisVectors.row(iR);
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
  writeWannierCoupling(context, GammaW, numElectrons, bravaisDegeneracies,
                       bravaisVectors, kMesh);
}

void ElElToPhoebeApp::writeWannierCoupling(
    Context &context,
    Eigen::Tensor<std::complex<double>, 7> &gWannier,
    const int &numElectrons,
    const Eigen::VectorXd &elDegeneracies,
    const Eigen::MatrixXd &elBravaisVectors,
    const Eigen::Vector3i &kMesh) {

#ifdef HDF5_AVAIL
  if (mpi->getSize() < 4) {
    // Note: this HDF5 had already been reported and being worked on.
    // It's beyond the purpose of Phoebe's project.
    Warning("HDF5 with <4 MPI process may crash (due to a "
            "library's bug),\nuse more MPI processes if that happens");
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
    mpi->barrier(); // wait for file to be overwritten

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

      HighFive::DataSet dnElectrons = file.createDataSet<int>(
          "/numElectrons", HighFive::DataSpace::From(numElectrons));
      dnElectrons.write(numElectrons);// # of occupied wannier functions

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
  throwErrorIfUnset(context.getQuantumEspressoPrefix(),
                    "QuantumEspressoPrefix");
  throwErrorIfUnset(context.getWannier90Prefix(), "Wannier90Prefix");
}
