#include "bandstructure.h"
#include "eigen.h"
#include "elph_qe_to_phoebe_app.h"
#include "interaction_elph.h"
#include "io.h"
#include "qe_input_parser.h"
#include <exception>
#include <iomanip>
#include <sstream>
#include <string>

#ifdef HDF5_AVAIL
#include <highfive/H5Easy.hpp>
#endif

Eigen::Tensor<std::complex<double>, 5>
ElPhQeToPhoebeApp::BlochToWannierEfficient(
    Context &context, const Eigen::MatrixXd &energies,
    const Eigen::MatrixXd &kGridFull, const int &numIrrQPoints,
    const int &numQEBands, const Eigen::MatrixXd &elBravaisVectors,
    const Eigen::MatrixXd &phBravaisVectors,
    const Eigen::Tensor<std::complex<double>, 3> &uMatrices, Points &kPoints,
    Points &qPoints, Crystal &crystal, PhononH0 &phononH0) {

  int numModes = crystal.getNumAtoms() * 3;
  int numBands = uMatrices.dimension(0);
  int numWannier = uMatrices.dimension(1);
  int numKPoints = kPoints.getNumPoints();
  int numQPoints = qPoints.getNumPoints();
  int numElBravaisVectors = elBravaisVectors.cols();
  int numPhBravaisVectors = phBravaisVectors.cols();

  std::string wannierPrefix = context.getWannier90Prefix();
  int bandsOffset = computeOffset(energies, wannierPrefix);

  Eigen::VectorXi ikMap(numKPoints);
#pragma omp parallel for default(none)                                         \
    shared(numKPoints, kGridFull, kPoints, ikMap)
  for (int ikOld = 0; ikOld < numKPoints; ikOld++) {
    Eigen::Vector3d kOld = kGridFull.col(ikOld);
    int ikNew = kPoints.getIndex(kOld);
    ikMap(ikOld) = ikNew;
  }
  std::string phoebePrefixQE = context.getQuantumEspressoPrefix();

  std::array<Eigen::Index, 5> zeros;
  for (auto &s : zeros) {
    s = 0;
  }

  bool usePolarCorrection = false;
  Eigen::Matrix3d dielectricMatrix = phononH0.getDielectricMatrix();
  if (dielectricMatrix.squaredNorm() > 1.0e-10) { // i.e. if dielectricMatrix
    // wasn't computed
    if (crystal.getNumSpecies() > 1) { // otherwise polar correction = 0
      usePolarCorrection = true;
    }
  }
  if (usePolarCorrection && mpi->mpiHead()) {
    std::cout << "Polar correction" << std::endl;
  }

  //  Eigen::Tensor<std::complex<double>, 5> gWannier(numWannier, numWannier,
  //                                                  numModes,
  //                                                  numPhBravaisVectors,
  //                                                  numElBravaisVectors);
  //  gWannier.setZero();

  ParallelMatrix<double> aux(numElBravaisVectors, 1, mpi->getSize(), 1);

  Eigen::Tensor<std::complex<double>, 5> gWannierPara(
      numWannier, numWannier, numModes, numPhBravaisVectors, aux.localRows());
  gWannierPara.setZero();

  LoopPrint loopPrint("Wannier transform of coupling", "irreducible q-points",
                      mpi->divideWorkIter(numIrrQPoints).size());
  for (int iqIrr : mpi->divideWorkIter(numIrrQPoints)) {
    loopPrint.update(false);
    // for (int iqIrr = 0; iqIrr < numIrrQPoints; iqIrr++) {
    //    std::cout << iqIrr << "\n";

    auto t =
        readChunkGFromQE(iqIrr, context, kPoints, numModes, numQEBands, ikMap);
    auto gStarTmp = std::get<0>(t);
    auto phononEigenvectorsStar = std::get<1>(t);
    auto qStar = std::get<3>(t);

    // Note: the tensor read from file contains
    // gFull(ib1, ib2, nu, ik, iq)
    // = < k+q,ib1 |  dV_{q,nu}  |  k,ib2  >
    // where k and q run over the full mesh.
    // ikMap takes care of the fact that k-points in QE have a different
    // order then phoebe.

    // reorder the q/k indices
    Eigen::Tensor<std::complex<double>, 5> gStar(numBands, numBands, numModes,
                                                 numKPoints, numQPoints);
    for (int iqStar = 0; iqStar < qStar.cols(); iqStar++) {
      for (int nu = 0; nu < numModes; nu++) {
        for (int ik = 0; ik < numKPoints; ik++) {
          for (int ib2 = 0; ib2 < numWannier; ib2++) {
            for (int ib1 = 0; ib1 < numWannier; ib1++) {
              gStar(ib1, ib2, nu, ik, iqStar) = gStarTmp(
                  bandsOffset + ib1, bandsOffset + ib2, nu, ik, iqStar);
            }
          }
        }
      }
    }
    gStarTmp.resize(zeros);
    int numQStar = qStar.cols();

    if (usePolarCorrection) {
      // we need to subtract the polar correction
      // this contribution will be reinstated during the interpolation
      auto volume = crystal.getVolumeUnitCell();
      auto reciprocalUnitCell = crystal.getReciprocalUnitCell();
      auto bornCharges = phononH0.getBornCharges();
      auto atomicPositions = crystal.getAtomicPositions();
      auto qCoarseMesh = phononH0.getCoarseGrid();

      for (int iq = 0; iq < numQStar; iq++) {
        Eigen::Vector3d qCrystal = qStar.col(iq);
        Eigen::Vector3d q = qPoints.crystalToCartesian(qCrystal);
        if (q.norm() > 1.0e-8) {

          Eigen::MatrixXcd ev3(numModes, numModes);
          for (int j = 0; j < numModes; j++) {
            for (int i = 0; i < numModes; i++) {
              ev3(i, j) = phononEigenvectorsStar(i, j, iq);
            }
          }

          for (int ik = 0; ik < numKPoints; ik++) {
            Eigen::Vector3d k =
                kPoints.getPointCoordinates(ik, Points::cartesianCoordinates);

            // Coordinates and index of k+q point
            Eigen::Vector3d kq = k + q;
            Eigen::Vector3d kqCrystal = kPoints.cartesianToCrystal(kq);
            int ikq = kPoints.getIndex(kqCrystal);

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
                q, ev1, ev2, ev3, volume, reciprocalUnitCell, dielectricMatrix,
                bornCharges, atomicPositions, qCoarseMesh);
            for (int nu = 0; nu < numModes; nu++) {
              for (int j = 0; j < numBands; j++) {
                for (int i = 0; i < numBands; i++) {
                  gStar(i, j, nu, ik, iq) -= v(i, j, nu);
                }
              }
            }
          }
        }
      }
    }

    Eigen::Tensor<std::complex<double>, 5> gFullTmp(
        numWannier, numWannier, numModes, numKPoints, numQStar);
    gFullTmp.setZero();

    for (int iq = 0; iq < numQStar; iq++) {
      Eigen::Vector3d qCrystal = qStar.col(iq);
      Eigen::Vector3d q = qPoints.crystalToCartesian(qCrystal);
      for (int ik = 0; ik < numKPoints; ik++) {
        Eigen::Vector3d k =
            kPoints.getPointCoordinates(ik, Points::cartesianCoordinates);

        // Coordinates and index of k+q point
        Eigen::Vector3d kq = k + q;
        Eigen::Vector3d kqCrystal = kPoints.cartesianToCrystal(kq);
        int ikq = kPoints.getIndex(kqCrystal);

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

#pragma omp parallel default(none) shared(numModes, numWannier, numBands, uKq, \
                                          gStar, gFullTmp, ik, iq, uKDagger)
        {
          Eigen::Tensor<std::complex<double>, 3> tmp(numWannier, numBands,
                                                     numModes);
          tmp.setZero();
#pragma omp for nowait collapse(4)
          for (int nu = 0; nu < numModes; nu++) {
            for (int i = 0; i < numWannier; i++) {
              for (int j = 0; j < numBands; j++) {
                for (int l = 0; l < numBands; l++) {
                  // ukq has size(numWannier, numBands)
                  // gFull has size numBands, numBands, ...
                  tmp(i, j, nu) += uKq(l, i) * gStar(l, j, nu, ik, iq);
                }
              }
            }
          }
          Eigen::Tensor<std::complex<double>, 3> tmp2(numWannier, numWannier,
                                                      numModes);
          tmp2.setZero();
#pragma omp for nowait collapse(4)
          for (int nu = 0; nu < numModes; nu++) {
            for (int i = 0; i < numWannier; i++) {
              for (int j = 0; j < numWannier; j++) {
                for (int l = 0; l < numBands; l++) {
                  tmp2(i, j, nu) += tmp(i, l, nu) * uKDagger(j, l);
                }
              }
            }
          }

#pragma omp critical
          for (int nu = 0; nu < numModes; nu++) {
            for (int i = 0; i < numWannier; i++) {
              for (int j = 0; j < numWannier; j++) {
                gFullTmp(i, j, nu, ik, iq) += tmp2(i, j, nu);
              }
            }
          }
        }
      } // ik
    }   // iq
    gStar.reshape(zeros);

    // Fourier transform on the electronic coordinates
    Eigen::Tensor<std::complex<double>, 5> gMixed(
        numWannier, numWannier, numModes, numElBravaisVectors, numQStar);
    gMixed.setZero();

    {
      Eigen::MatrixXcd phases(numKPoints, numElBravaisVectors);
      phases.setZero();
#pragma omp parallel for default(none)                                         \
    shared(numKPoints, kPoints, numElBravaisVectors, elBravaisVectors, phases, \
           mpi, complexI)
      for (int ik = 0; ik < numKPoints; ik++) {
        Eigen::Vector3d k =
            kPoints.getPointCoordinates(ik, Points::cartesianCoordinates);
        for (int iR = 0; iR < numElBravaisVectors; iR++) {
          double arg = k.dot(elBravaisVectors.col(iR));
          phases(ik, iR) = exp(-complexI * arg) / double(numKPoints);
        }
      }

      for (int iq = 0; iq < numQStar; iq++) {
#pragma omp parallel default(none)                                             \
    shared(iq, gFullTmp, phases, numElBravaisVectors, numKPoints, numModes,    \
           numWannier, gMixed)
        {
          Eigen::Tensor<std::complex<double>, 4> tmp(
              numWannier, numWannier, numModes, numElBravaisVectors);
          tmp.setZero();
#pragma omp for nowait
          for (int iR = 0; iR < numElBravaisVectors; iR++) {
            for (int ik = 0; ik < numKPoints; ik++) {
              for (int nu = 0; nu < numModes; nu++) {
                for (int j = 0; j < numWannier; j++) {
                  for (int i = 0; i < numWannier; i++) {
                    tmp(i, j, nu, iR) +=
                        gFullTmp(i, j, nu, ik, iq) * phases(ik, iR);
                  }
                }
              }
            }
          }
#pragma omp critical
          for (int iR = 0; iR < numElBravaisVectors; iR++) {
            for (int nu = 0; nu < numModes; nu++) {
              for (int j = 0; j < numWannier; j++) {
                for (int i = 0; i < numWannier; i++) {
                  gMixed(i, j, nu, iR, iq) += tmp(i, j, nu, iR);
                }
              }
            }
          }
        }
      } // iq
    }
    gFullTmp.reshape(zeros);

    Eigen::Tensor<std::complex<double>, 5> gWannierTmp(
        numWannier, numWannier, numModes, numElBravaisVectors, numQStar);
    gWannierTmp.setZero();
    {
      Eigen::Tensor<std::complex<double>, 3> uQM1s(numModes, numModes,
                                                   numQStar);
      uQM1s.setZero();
      for (int iq = 0; iq < numQStar; iq++) {
        Eigen::MatrixXcd uQ(numModes, numModes);
        for (int nu2 = 0; nu2 < numModes; nu2++) {
          for (int nu = 0; nu < numModes; nu++) {
            uQ(nu, nu2) = phononEigenvectorsStar(nu, nu2, iq);
          }
        }
        auto uQM1 = uQ.inverse();
        for (int nu2 = 0; nu2 < numModes; nu2++) {
          for (int nu = 0; nu < numModes; nu++) {
            uQM1s(nu, nu2, iq) = uQM1(nu, nu2);
          }
        }
        // this isn't equal to the adjoint, due to mass renormalization
        // should be parallelized with OMP already
      }
      for (int iq = 0; iq < numQStar; iq++) {
        for (int nu = 0; nu < numModes; nu++) {
          for (int nu2 = 0; nu2 < numModes; nu2++) {
#pragma omp parallel for collapse(3) default(none) shared(                     \
    numElBravaisVectors, numWannier, gMixed, uQM1s, iq, nu, nu2, gWannierTmp)
            for (int irE = 0; irE < numElBravaisVectors; irE++) {
              for (int i = 0; i < numWannier; i++) {
                for (int j = 0; j < numWannier; j++) {
                  gWannierTmp(i, j, nu, irE, iq) +=
                      gMixed(i, j, nu2, irE, iq) * uQM1s(nu2, nu, iq);
                }
              }
            }
          }
        }
      }
    }
    gMixed.reshape(zeros);

    {
      Eigen::MatrixXcd phases(numPhBravaisVectors, numQStar);
      phases.setZero();
#pragma omp parallel for default(none)                                         \
    shared(qStar, mpi, numQPoints, numPhBravaisVectors, complexI, phases,      \
           qPoints, phBravaisVectors, numQStar)
      for (int iq = 0; iq < numQStar; iq++) {
        Eigen::Vector3d qCrystal = qStar.col(iq);
        Eigen::Vector3d q = qPoints.crystalToCartesian(qCrystal);
        for (int irP = 0; irP < numPhBravaisVectors; irP++) {
          double arg = q.dot(phBravaisVectors.col(irP));
          phases(irP, iq) = exp(-complexI * arg) / double(numQPoints);
        }
      }

      for (int irE = 0; irE < numElBravaisVectors; irE++) {
        Eigen::Tensor<std::complex<double>, 4> tmp(
            numWannier, numWannier, numModes, numPhBravaisVectors);
        tmp.setZero();
#pragma omp parallel for collapse(4) default(none)                             \
    shared(numQStar, numPhBravaisVectors, phases, numModes, numWannier,        \
           gWannierTmp, tmp, irE)
        for (int irP = 0; irP < numPhBravaisVectors; irP++) {
          for (int nu = 0; nu < numModes; nu++) {
            for (int j = 0; j < numWannier; j++) {
              for (int i = 0; i < numWannier; i++) {
                for (int iq = 0; iq < numQStar; iq++) {
                  tmp(i, j, nu, irP) +=
                      phases(irP, iq) * gWannierTmp(i, j, nu, irE, iq);
                }
              }
            }
          }
        }
        mpi->allReduceSum(&tmp);

        if (aux.indicesAreLocal(irE, 0)) { // is local
          int irELocal = aux.global2Local(irE, 0);
          for (int irP = 0; irP < numPhBravaisVectors; irP++) {
            for (int nu = 0; nu < numModes; nu++) {
              for (int j = 0; j < numWannier; j++) {
                for (int i = 0; i < numWannier; i++) {
                  gWannierPara(i, j, nu, irP, irELocal) += tmp(i, j, nu, irP);
                }
              }
            }
          }
        }
      }
    }
  }
  loopPrint.close();

  //  // put results back in common tensor
  //  Eigen::Tensor<std::complex<double>, 5> gWannier(numWannier, numWannier,
  //                                                  numModes,
  //                                                  numPhBravaisVectors,
  //                                                  numElBravaisVectors);
  //  gWannier.setZero();
  //  for (int irE = 0; irE < numElBravaisVectors; irE++) {
  //    if (aux.indicesAreLocal(irE, 0)) { // is local
  //      int irELocal = aux.global2Local(irE, 0);
  //      for (int irP = 0; irP < numPhBravaisVectors; irP++) {
  //        for (int nu = 0; nu < numModes; nu++) {
  //          for (int i = 0; i < numWannier; i++) {
  //            for (int j = 0; j < numWannier; j++) {
  //              gWannier(i, j, nu, irP, irE) +=
  //                  gWannierPara(i, j, nu, irP, irELocal);
  //            }
  //          }
  //        }
  //      }
  //    }
  //  }
  //  mpi->allReduceSum(&gWannier);

  if (mpi->mpiHead()) {
    std::cout << "Done Wannier-transform of g\n" << std::endl;
  }

  return gWannierPara;
}

Eigen::Tensor<std::complex<double>, 5> ElPhQeToPhoebeApp::blochToWannier(
    const Eigen::MatrixXd &elBravaisVectors,
    const Eigen::MatrixXd &phBravaisVectors,
    Eigen::Tensor<std::complex<double>, 5> &gFull,
    const Eigen::Tensor<std::complex<double>, 3> &uMatrices,
    const Eigen::Tensor<std::complex<double>, 3> &phEigenvectors,
    Points &kPoints, Points &qPoints, Crystal &crystal, PhononH0 &phononH0) {

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
  Eigen::Matrix3d dielectricMatrix = phononH0.getDielectricMatrix();
  if (dielectricMatrix.squaredNorm() > 1.0e-10) { // i.e. if dielectricMatrix
                                                  // wasn't computed
    if (crystal.getNumSpecies() > 1) { // otherwise polar correction = 0
      usePolarCorrection = true;
    }
  }

  if (usePolarCorrection) {
    if (mpi->mpiHead()) {
      std::cout << "Polar correction" << std::endl;
    }

    // we need to subtract the polar correction
    // this contribution will be reinstated during the interpolation
    auto volume = crystal.getVolumeUnitCell();
    auto reciprocalUnitCell = crystal.getReciprocalUnitCell();
    auto bornCharges = phononH0.getBornCharges();
    auto atomicPositions = crystal.getAtomicPositions();
    auto qCoarseMesh = phononH0.getCoarseGrid();

    for (int iq = 0; iq < numQPoints; iq++) {
      Eigen::Vector3d q =
          qPoints.getPointCoordinates(iq, Points::cartesianCoordinates);
      if (q.norm() > 1.0e-8) {

        Eigen::MatrixXcd ev3(numModes, numModes);
        for (int j = 0; j < numModes; j++) {
          for (int i = 0; i < numModes; i++) {
            ev3(i, j) = phEigenvectors(i, j, iq);
          }
        }

        for (int ik = 0; ik < numKPoints; ik++) {
          Eigen::Vector3d k =
              kPoints.getPointCoordinates(ik, Points::cartesianCoordinates);

          // Coordinates and index of k+q point
          Eigen::Vector3d kq = k + q;
          Eigen::Vector3d kqCrystal = kPoints.cartesianToCrystal(kq);
          int ikq = kPoints.getIndex(kqCrystal);

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
              q, ev1, ev2, ev3, volume, reciprocalUnitCell, dielectricMatrix,
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

  for (int iq : mpi->divideWorkIter(numQPoints)) {
    Eigen::Vector3d q =
        qPoints.getPointCoordinates(iq, Points::cartesianCoordinates);
    for (int ik = 0; ik < numKPoints; ik++) {
      Eigen::Vector3d k =
          kPoints.getPointCoordinates(ik, Points::cartesianCoordinates);

      // Coordinates and index of k+q point
      Eigen::Vector3d kq = k + q;
      Eigen::Vector3d kqCrystal = kPoints.cartesianToCrystal(kq);
      int ikq = kPoints.getIndex(kqCrystal);

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

#pragma omp parallel default(none) shared(numModes, numWannier, numBands, uKq, \
                                          gFull, gFullTmp, ik, iq, uKDagger)
      {
        Eigen::Tensor<std::complex<double>, 3> tmp(numWannier, numBands,
                                                   numModes);
        tmp.setZero();
#pragma omp for nowait collapse(4)
        for (int nu = 0; nu < numModes; nu++) {
          for (int i = 0; i < numWannier; i++) {
            for (int j = 0; j < numBands; j++) {
              for (int l = 0; l < numBands; l++) {
                // ukq has size(numWannier, numBands)
                // gFull has size numBands, numBands, ...
                tmp(i, j, nu) += uKq(l, i) * gFull(l, j, nu, ik, iq);
              }
            }
          }
        }
        Eigen::Tensor<std::complex<double>, 3> tmp2(numWannier, numWannier,
                                                    numModes);
        tmp2.setZero();
#pragma omp for nowait collapse(4)
        for (int nu = 0; nu < numModes; nu++) {
          for (int i = 0; i < numWannier; i++) {
            for (int j = 0; j < numWannier; j++) {
              for (int l = 0; l < numBands; l++) {
                tmp2(i, j, nu) += tmp(i, l, nu) * uKDagger(j, l);
              }
            }
          }
        }

#pragma omp critical
        for (int nu = 0; nu < numModes; nu++) {
          for (int i = 0; i < numWannier; i++) {
            for (int j = 0; j < numWannier; j++) {
              gFullTmp(i, j, nu, ik, iq) += tmp2(i, j, nu);
            }
          }
        }
      }
    } // ik
  }   // iq
  mpi->allReduceSum(&gFullTmp);
  gFull.reshape(zeros);

  if (mpi->mpiHead()) {
    std::cout << "Electronic Fourier Transform" << std::endl;
  }

  // Fourier transform on the electronic coordinates
  Eigen::Tensor<std::complex<double>, 5> gMixed(
      numWannier, numWannier, numModes, numElBravaisVectors, numQPoints);
  gMixed.setZero();

  {
    Eigen::MatrixXcd phases(numKPoints, numElBravaisVectors);
    phases.setZero();
#pragma omp parallel for default(none)                                         \
    shared(numKPoints, kPoints, numElBravaisVectors, elBravaisVectors, phases, \
           mpi, complexI)
    for (int ik : mpi->divideWorkIter(numKPoints)) {
      Eigen::Vector3d k =
          kPoints.getPointCoordinates(ik, Points::cartesianCoordinates);
      for (int iR = 0; iR < numElBravaisVectors; iR++) {
        double arg = k.dot(elBravaisVectors.col(iR));
        phases(ik, iR) = exp(-complexI * arg) / double(numKPoints);
      }
    }
    mpi->allReduceSum(&phases);

    for (int iq : mpi->divideWorkIter(numQPoints)) {
#pragma omp parallel default(none)                                             \
    shared(iq, gFullTmp, phases, numElBravaisVectors, numKPoints, numModes,    \
           numWannier, gMixed)
      {
        Eigen::Tensor<std::complex<double>, 4> tmp(
            numWannier, numWannier, numModes, numElBravaisVectors);
        tmp.setZero();
#pragma omp for nowait
        for (int iR = 0; iR < numElBravaisVectors; iR++) {
          for (int ik = 0; ik < numKPoints; ik++) {
            for (int nu = 0; nu < numModes; nu++) {
              for (int j = 0; j < numWannier; j++) {
                for (int i = 0; i < numWannier; i++) {
                  tmp(i, j, nu, iR) +=
                      gFullTmp(i, j, nu, ik, iq) * phases(ik, iR);
                }
              }
            }
          }
        }
#pragma omp critical
        for (int iR = 0; iR < numElBravaisVectors; iR++) {
          for (int nu = 0; nu < numModes; nu++) {
            for (int j = 0; j < numWannier; j++) {
              for (int i = 0; i < numWannier; i++) {
                gMixed(i, j, nu, iR, iq) += tmp(i, j, nu, iR);
              }
            }
          }
        }
      }
    } // iq
    mpi->allReduceSum(&gMixed);
  }
  gFullTmp.reshape(zeros);

  if (mpi->mpiHead()) {
    std::cout << "Phonon rotation" << std::endl;
  }

  Eigen::Tensor<std::complex<double>, 5> gWannierTmp(
      numWannier, numWannier, numModes, numElBravaisVectors, numQPoints);
  gWannierTmp.setZero();
  {
    Eigen::Tensor<std::complex<double>, 3> uQM1s(numModes, numModes,
                                                 numQPoints);
    uQM1s.setZero();
    for (int iq : mpi->divideWorkIter(numQPoints)) {
      Eigen::MatrixXcd uQ(numModes, numModes);
      for (int nu2 = 0; nu2 < numModes; nu2++) {
        for (int nu = 0; nu < numModes; nu++) {
          uQ(nu, nu2) = phEigenvectors(nu, nu2, iq);
        }
      }
      auto uQM1 = uQ.inverse();
      for (int nu2 = 0; nu2 < numModes; nu2++) {
        for (int nu = 0; nu < numModes; nu++) {
          uQM1s(nu, nu2, iq) = uQM1(nu, nu2);
        }
      }
      // this isn't equal to the adjoint, due to mass renormalization
      // should be parallelized with OMP already
    }
    mpi->allReduceSum(&uQM1s);
    for (int iq : mpi->divideWorkIter(numQPoints)) {
      for (int nu = 0; nu < numModes; nu++) {
        for (int nu2 = 0; nu2 < numModes; nu2++) {
#pragma omp parallel for collapse(3) default(none) shared(                     \
    numElBravaisVectors, numWannier, gMixed, uQM1s, iq, nu, nu2, gWannierTmp)
          for (int irE = 0; irE < numElBravaisVectors; irE++) {
            for (int i = 0; i < numWannier; i++) {
              for (int j = 0; j < numWannier; j++) {
                gWannierTmp(i, j, nu, irE, iq) +=
                    gMixed(i, j, nu2, irE, iq) * uQM1s(nu2, nu, iq);
              }
            }
          }
        }
      }
    }
    mpi->allReduceSum(&gWannierTmp);
  }
  gMixed.reshape(zeros);

  if (mpi->mpiHead()) {
    std::cout << "Phonon Fourier Transform" << std::endl;
  }

  Eigen::Tensor<std::complex<double>, 5> gWannier(numWannier, numWannier,
                                                  numModes, numPhBravaisVectors,
                                                  numElBravaisVectors);
  gWannier.setZero();
  {
    Eigen::MatrixXcd phases(numPhBravaisVectors, numQPoints);
    phases.setZero();
#pragma omp parallel for default(none)                                         \
    shared(mpi, numQPoints, numPhBravaisVectors, complexI, phases, qPoints,    \
           phBravaisVectors)
    for (int iq : mpi->divideWorkIter(numQPoints)) {
      Eigen::Vector3d q =
          qPoints.getPointCoordinates(iq, Points::cartesianCoordinates);
      for (int irP = 0; irP < numPhBravaisVectors; irP++) {
        double arg = q.dot(phBravaisVectors.col(irP));
        phases(irP, iq) = exp(-complexI * arg) / double(numQPoints);
      }
    }
    mpi->allReduceSum(&phases);

    for (int irE : mpi->divideWorkIter(numElBravaisVectors)) {
#pragma omp parallel default(none)                                             \
    shared(numQPoints, numPhBravaisVectors, numModes, numWannier, phases,      \
           gWannier, irE, gWannierTmp)
      {
        Eigen::Tensor<std::complex<double>, 4> tmp(
            numWannier, numWannier, numModes, numPhBravaisVectors);
        tmp.setZero();
#pragma omp for nowait
        for (int iq = 0; iq < numQPoints; iq++) {
          for (int irP = 0; irP < numPhBravaisVectors; irP++) {
            for (int nu = 0; nu < numModes; nu++) {
              for (int j = 0; j < numWannier; j++) {
                for (int i = 0; i < numWannier; i++) {
                  tmp(i, j, nu, irP) +=
                      phases(irP, iq) * gWannierTmp(i, j, nu, irE, iq);
                }
              }
            }
          }
        }
#pragma omp critical
        for (int irP = 0; irP < numPhBravaisVectors; irP++) {
          for (int nu = 0; nu < numModes; nu++) {
            for (int i = 0; i < numWannier; i++) {
              for (int j = 0; j < numWannier; j++) {
                gWannier(i, j, nu, irP, irE) += tmp(i, j, nu, irP);
              }
            }
          }
        }
      }
    }
    mpi->allReduceSum(&gWannier);
  }
  gWannierTmp.reshape(zeros);

  if (mpi->mpiHead()) {
    std::cout << "Done Wannier-transform of g\n" << std::endl;
  }

  return gWannier;
}

Eigen::Tensor<std::complex<double>, 3>
ElPhQeToPhoebeApp::setupRotationMatrices(const std::string &wannierPrefix,
                                         Points &fullPoints) {
  std::string line;

  if (wannierPrefix.empty()) {
    Error("Must provide an input H0 file name");
  }

  std::string fileName = wannierPrefix + "_u.mat";

  // open input file
  std::ifstream infile(fileName);
  if (not infile.is_open()) {
    Error("U-matrix file not found");
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

    Eigen::Vector3d thisK; // vector in crystal coordinates
    infile >> thisK(0) >> thisK(1) >> thisK(2);

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
    thisK << x, y, z; // vector in crystal coordinates

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
  for (int ik = 0; ik < numPoints; ik++) {

    Eigen::MatrixXcd a(numBands, numWannier);
    for (int i = 0; i < numBands; i++) {
      for (int j = 0; j < numWannier; j++) {
        a(i, j) = uMatrixDis(i, j, ik);
      }
    }

    Eigen::MatrixXcd b(numWannier, numWannier);
    for (int i = 0; i < numWannier; i++) {
      for (int j = 0; j < numWannier; j++) {
        b(i, j) = uMatrix(i, j, ik);
      }
    }

    for (int i = 0; i < numBands; i++) {
      for (int j = 0; j < numWannier; j++) {
        for (int k = 0; k < numWannier; k++) {
          u(i, j, ik) += b(k, j) * a(i, k);
        }
      }
    }
  } // ik
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
      Error("Expecting first coarse grid k-point in Wannier90 to be gamma");
    }
  }

  // read .eig file to get energies

  std::vector<double> energiesWannierAtZero;
  {
    std::string eigFileName = wannierPrefix + ".eig";
    std::ifstream eigenFile(eigFileName);
    int ib, ik;
    double x;
    while (eigenFile >> ib >> ik >> x) {
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
  int possibleValues = numFull - numBandsWannier;
  Eigen::VectorXd difference(possibleValues);
  difference.setZero();
  for (int i = 0; i < possibleValues; i++) {
    for (int ib = 0; ib < numBandsWannier; ib++) {
      difference(i) =
          pow(energiesQEAtZero(ib + i) - energiesWannierAtZero[ib], 2);
    }
  }

  // offset = index of min difference
  int offset = -1;
  for (int i = 0; i < possibleValues; i++) {
    if (difference(i) == difference.minCoeff()) {
      offset = i + 1;
      break;
    }
  }

  if (possibleValues == 0) {
    offset = 0;
  }

  if (offset == -1) {
    Error("Bands offset not found");
  }

  return offset;
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
   *    on the original grid of k-points
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
    auto kCrystal = kPoints.getPointCoordinates(ik);
    kCrystal(0) *= kMesh(0);
    kCrystal(1) *= kMesh(1);
    kCrystal(2) *= kMesh(2);

    int ikOld = int(kCrystal[0] * kMesh(2) * kMesh(1) + kCrystal[1] * kMesh(2) +
                    kCrystal[2]);
    {
      std::string eigFileName = wannierPrefix + ".eig";
      std::ifstream eigenFile(eigFileName);
      int ib, ikk;
      double x;
      while (eigenFile >> ib >> ikk >> x) {
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
    auto k1C = kPoints.getPointCoordinates(ik1, Points::cartesianCoordinates);

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
    auto k1C = kPoints.getPointCoordinates(ik, Points::cartesianCoordinates);
    auto t3 = electronH0.diagonalizeFromCoordinates(k1C);
    auto en = std::get<0>(t3);
    auto u = std::get<1>(t3);

    Eigen::MatrixXcd h0K(numWannier, numWannier);
    h0K.setZero();
    for (int iR = 0; iR < elBravaisVectors.cols(); iR++) {
      Eigen::Vector3d R = elBravaisVectors.col(iR);
      double arg = k1C.dot(R);
      std::complex<double> phase = exp(complexI * arg) / elDegeneracies(iR);
      for (int m = 0; m < numWannier; m++) {
        for (int n = 0; n < numWannier; n++) {
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

  int numPhBands = int(phononH0.getNumBands());

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
            int i = compress2Indices(k1, iCart, numAtoms, 3);
            norm(ib1, ib2) +=
                phEigenvectors(i, ib1, iq) * sqrt(atomicMasses(k1)) *
                phEigenvectors(i, ib2, iq) * sqrt(atomicMasses(k1));
          }
        }
      }
      norm(ib1) -= 1.; // It should be an identity matrix
    }
    for (int ib1 = 0; ib1 < numPhBands; ib1++) {
      assert(abs(norm(ib1)) < 1.0e-6);
    }
  }

  // FT to Wannier representation

  Eigen::Tensor<std::complex<double>, 5> h0R(
      numAtoms * numAtoms * phBravaisVectors.size(), numAtoms, numAtoms, 3, 3);
  h0R.setZero();

  for (int iq = 0; iq < qPoints.getNumPoints(); iq++) {
    auto qC = qPoints.getPointCoordinates(iq, Points::cartesianCoordinates);
    qC = qPoints.bzToWs(qC, Points::cartesianCoordinates);

    // u has size (numBands, numWannier, numKPoints)
    Eigen::MatrixXcd uK(numPhBands, numPhBands);
    for (int k1 = 0; k1 < numAtoms; k1++) {
      for (int iCart : {0, 1, 2}) {
        int i = compress2Indices(k1, iCart, numAtoms, 3);
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

    for (int iR = 0; iR < phBravaisVectors.cols(); iR++) {
      Eigen::Vector3d R0 = phBravaisVectors.col(iR);
      for (int k1 = 0; k1 < numAtoms; k1++) {
        for (int k2 = 0; k2 < numAtoms; k2++) {
          // Eigen::Vector3d R = R0; // - atomicPositions.col(k1)
          //+ atomicPositions.col(k2);
          double arg = qC.dot(R0);
          std::complex<double> phase =
              exp(-complexI * arg) / double(qPoints.getNumPoints());
          for (int iCart : {0, 1, 2}) {
            for (int jCart : {0, 1, 2}) {
              int m = compress2Indices(k1, iCart, numAtoms, 3);
              int n = compress2Indices(k2, jCart, numAtoms, 3);
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
    double imaginarySum = 0.;
    for (int iR0 = 0; iR0 < phBravaisVectors.cols(); iR0++) {
      for (int k1 = 0; k1 < numAtoms; k1++) {
        for (int k2 = 0; k2 < numAtoms; k2++) {
          for (int i : {0, 1, 2}) {
            for (int j : {0, 1, 2}) {
              double x = std::real(h0R(iR0, k1, k2, i, j));
              realSum += pow(x, 2);
              imaginarySum += pow(std::imag(h0R(iR0, k1, k2, i, j)), 2);
              // set to zero the imaginary part to clean noise
              // this is also what QE does
              h0R(iR0, k1, k2, i, j) = {x, 0.};
            }
          }
        }
      }
    }
    // I want the imaginary part to be much smaller than the real
    assert(imaginarySum * pow(10, 6) < realSum);
  }

  //--------------------------------------------------------------------------
  // FFT back

  for (int iq = 0; iq < qPoints.getNumPoints(); iq++) {
    // get U
    auto qC = qPoints.getPointCoordinates(iq, Points::cartesianCoordinates);
    auto t = phononH0.diagonalizeFromCoordinates(qC, false);
    auto en = std::get<0>(t);
    auto u = std::get<0>(t);

    Eigen::MatrixXcd hWK(numPhBands, numPhBands);
    hWK.setZero();
    for (int iR = 0; iR < phBravaisVectors.cols(); iR++) {
      Eigen::Vector3d R0 = phBravaisVectors.col(iR);
      for (int k1 = 0; k1 < numAtoms; k1++) {
        for (int k2 = 0; k2 < numAtoms; k2++) {
          // Eigen::Vector3d R = R0; // - atomicPositions.col(k1)
          //+ atomicPositions.col(k2);
          double arg = qC.dot(R0);
          std::complex<double> phase = exp(complexI * arg) / phDegeneracies(iR);
          for (int iCart : {0, 1, 2}) {
            for (int jCart : {0, 1, 2}) {
              int m = compress2Indices(k1, iCart, numAtoms, 3);
              int n = compress2Indices(k2, jCart, numAtoms, 3);
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
  int numKPoints = int(kPoints.getNumPoints());
  int numModes = int(phononH0.getNumBands());

// needed by ::parse()
#ifdef HDF5_AVAIL
  context.setElphFileName(context.getQuantumEspressoPrefix() +
                          ".phoebe.elph.hdf5");
#else
  context.setElphFileName(context.getQuantumEspressoPrefix() +
                          ".phoebe.elph.dat");
#endif

  auto couplingElPh = InteractionElPhWan::parse(context, crystal, &phononH0);

  for (int ik1 = 0; ik1 < numKPoints; ik1++) {
    Eigen::Vector3d k1C =
        kPoints.getPointCoordinates(ik1, Points::cartesianCoordinates);
    for (int ik2 = 0; ik2 < numKPoints; ik2++) {
      Eigen::Vector3d k2C =
          kPoints.getPointCoordinates(ik2, Points::cartesianCoordinates);

      std::vector<Eigen::Vector3d> k2Cs;
      k2Cs.push_back(k2C);

      Eigen::Vector3d q3C = k2C - k1C;
      Eigen::Vector3d q3Crystal = qPoints.cartesianToCrystal(q3C);
      int iq3 = int(qPoints.getIndex(q3Crystal));
      std::vector<Eigen::Vector3d> q3Cs;
      q3Cs.push_back(q3C);

      auto ik1Index = WavevectorIndex(ik1);
      auto ik2Index = WavevectorIndex(ik2);

      Eigen::MatrixXcd eigenVector1 = bandStructure.getEigenvectors(ik1Index);
      Eigen::MatrixXcd eigenVector2 = bandStructure.getEigenvectors(ik2Index);
      std::vector<Eigen::MatrixXcd> eigenVectors2;
      eigenVectors2.push_back(eigenVector2);

      auto t = phononH0.diagonalizeFromCoordinates(q3C);
      auto eigenVector3 = std::get<1>(t);
      std::vector<Eigen::MatrixXcd> eigenVectors3;
      eigenVectors3.push_back(eigenVector3);

      couplingElPh.calcCouplingSquared(eigenVector1, eigenVectors2,
                                       eigenVectors3, k1C, k2Cs, q3Cs);
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
      // note that I change the meaning of the indices
      assert(abs((sum1 - sum2) / sum1) < 0.0001);
    }
  }
}

void ElPhQeToPhoebeApp::postProcessingWannier(
    Context &context, Crystal &crystal, PhononH0 &phononH0, Points &kPoints,
    Points &qPoints, int numQEBands, int numModes, int numIrrQPoints,
    int numElectrons, int numSpin, const Eigen::MatrixXd &energies,
    const Eigen::MatrixXd &kGridFull, const Eigen::Vector3i &kMesh,
    const Eigen::Vector3i &qMesh, bool runTests) {
  if (mpi->mpiHead()) {
    std::cout << "Starting Wannier post-processing\n" << std::endl;
  }

  std::string wannierPrefix = context.getWannier90Prefix();

  auto t2 = QEParser::parseElHarmonicWannier(context, &crystal);
  auto electronH0 = std::get<1>(t2);
  int numWannier = int(electronH0.getNumBands());

  //----------------------------------------------------------------------------
  // read Wannier90 rotation matrices

  Eigen::Tensor<std::complex<double>, 3> uMatrices;
  // uMatrices has size (numBands, numWannier, numKPoints)
  uMatrices = setupRotationMatrices(wannierPrefix, kPoints);
  int numBands = uMatrices.dimension(0);        // number of entangled bands
  assert(numWannier == uMatrices.dimension(1)); // number of entangled bands

  //----------------------------------------------------------------------------
  // figure out the number of occupied Wannier states

  int bandsOffset = computeOffset(energies, wannierPrefix);
  // Note: Quantum-ESPRESSO defines numElectrons as the number of occupied
  // states in the valence manifold (it doesn't tell us about core electrons,
  // which we probably can only infer from the atomic numbers)
  if (numSpin == 2) {
    Error("Spin is not currently supported");
  }
  // the factor 2 is the spin degeneracy factor. Change if spin support added
  int numFilledWannier = numElectrons - bandsOffset * 2;
  // note how we only allow numFilledWannier to be an integer
  // it can be an even or odd number, so be careful if dividing it by 2

  //----------------------------------------------------------------------------

  // Find the lattice vectors for the Fourier transforms

  auto t3 = crystal.buildWignerSeitzVectors(kMesh);
  Eigen::MatrixXd elBravaisVectors = std::get<0>(t3);
  Eigen::VectorXd elDegeneracies = std::get<1>(t3);

  auto t4 = crystal.buildWignerSeitzVectors(qMesh);
  Eigen::MatrixXd phBravaisVectors = std::get<0>(t4);
  Eigen::VectorXd phDegeneracies = std::get<1>(t4);

  //----------------------------------------------------------------------------

  Eigen::Tensor<std::complex<double>, 5> gWannier;
  if (context.getDistributedElPhCoupling()) {
    // read coupling from file
    auto t5 = readGFromQEFile(context, numModes, numBands, numWannier, kPoints,
                              qPoints, kGridFull, numIrrQPoints, numQEBands,
                              energies);
    auto gFull = std::get<0>(t5);          // (nBands,nBands,nModes,numK,numQ)
    auto phEigenvectors = std::get<1>(t5); // (numModes,numModes,numQPoints)
    auto phEnergies = std::get<2>(t5);     // (numModes,numQPoints)
    gWannier =
        blochToWannier(elBravaisVectors, phBravaisVectors, gFull, uMatrices,
                       phEigenvectors, kPoints, qPoints, crystal, phononH0);
  } else {
    gWannier =
        BlochToWannierEfficient(context, energies, kGridFull, numIrrQPoints,
                                numQEBands, elBravaisVectors, phBravaisVectors,
                                uMatrices, kPoints, qPoints, crystal, phononH0);
  }
  writeWannierCoupling(context, gWannier, numFilledWannier, numSpin, numModes,
                       numWannier, phDegeneracies, elDegeneracies,
                       phBravaisVectors, elBravaisVectors, qMesh, kMesh);

  //--------------------------------------------------------------------------

  if (runTests) {
    auto t6 = readGFromQEFile(context, numModes, numBands, numWannier, kPoints,
                              qPoints, kGridFull, numIrrQPoints, numQEBands,
                              energies);
    auto gFull = std::get<0>(t6);          // (nBands,nBands,nModes,numK,numQ)
    auto phEigenvectors = std::get<1>(t6); // (numModes,numModes,numQPoints)
    auto phEnergies = std::get<2>(t6);     // (numModes,numQPoints)

    testElectronicTransform(kPoints, wannierPrefix, elBravaisVectors, uMatrices,
                            elDegeneracies, electronH0);

    testPhononTransform(crystal, phononH0, qPoints, phEigenvectors,
                        phBravaisVectors, phDegeneracies, phEnergies);

    testBackTransform(context, phononH0, kPoints, qPoints, electronH0, crystal,
                      gFull);
  }
}

void ElPhQeToPhoebeApp::writeWannierCoupling(
    Context &context, Eigen::Tensor<std::complex<double>, 5> &gWannier,
    const int &numFilledWannier, const int &numSpin, const int &numModes,
    const int &numWannier, const Eigen::VectorXd &phDegeneracies,
    const Eigen::VectorXd &elDegeneracies,
    const Eigen::MatrixXd &phBravaisVectors,
    const Eigen::MatrixXd &elBravaisVectors, const Eigen::Vector3i &qMesh,
    const Eigen::Vector3i &kMesh) {
  if (mpi->mpiHead())
    std::cout << "\nStart writing el-ph coupling to file." << std::endl;
  std::string phoebePrefixQE = context.getQuantumEspressoPrefix();

  int numPhBravaisVectors = phDegeneracies.size();
  int numElBravaisVectors = elDegeneracies.size();

  bool matrixDistributed;
  if (gWannier.dimension(4) != numElBravaisVectors) {
    matrixDistributed = true;
  } else {
    matrixDistributed = false;
  }

  ParallelMatrix<double> aux(elDegeneracies.size(), 1, mpi->getSize(), 1);

#ifdef HDF5_AVAIL
  std::string outFileName = phoebePrefixQE + ".phoebe.elph.hdf5";
  // if the hdf5 file is there already, we want to delete it. Occasionally
  // these files seem to get stuck open when a process dies while writing to
  // them, (even if a python script dies) and then they can't be overwritten
  // properly.
  std::remove(&outFileName[0]);

  if (mpi->getSize() < 4) {
    // Note: this HDF5 had already been reported and being worked on.
    // It's beyond the purpose of Phoebe's project.
    Warning("HDF5 with 1 MPI process may crash (due to a "
            "library's bug),\nuse more MPI processes if that happens");
  }

  try {
    // need to open the files differently if MPI is available or not
    // NOTE: do not remove the braces inside this if -- the file must
    // go out of scope, so that it can be reopened/written by head for the
    // small quantities as in the next block.
#if defined(MPI_AVAIL) && !defined(HDF5_SERIAL)
    {
      // open the hdf5 file
      HighFive::File file(
          outFileName, HighFive::File::Overwrite,
          HighFive::MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));

      // flatten the tensor (tensor is not supported) and create the data set
      Eigen::VectorXcd gwan = Eigen::Map<Eigen::VectorXcd, Eigen::Unaligned>(
          gWannier.data(), gWannier.size());

      // note: gwan is distributed
      unsigned int globalSize = numWannier * numWannier * numModes *
                                numPhBravaisVectors * numElBravaisVectors;

      // Create the data-space to write gWannier to
      std::vector<size_t> dims(2);
      dims[0] = 1;
      dims[1] = size_t(globalSize);
      HighFive::DataSet dgwannier = file.createDataSet<std::complex<double>>(
          "/gWannier", HighFive::DataSpace(dims));

      Eigen::VectorXcd gwanSlice;
      size_t offset, numElements;
      if (matrixDistributed) {
        // here, no need to slice the gwan tensor (it's already distributed)
        // but we have to compute the right offsets to file.

        // start point and the number of elements to be written by this process
        offset = aux.getAllLocalRows()[0] * pow(numWannier, 2) *
                        numModes * numPhBravaisVectors;
        numElements = aux.getAllLocalRows().size() * pow(numWannier, 2) *
                             numModes * numPhBravaisVectors;
        gwanSlice = gwan;
      } else {
        // here we slice the gwan tensor (it's not distributed)

        // get the start and stop points of elements to be written by this process
        std::vector<int> workDivs = mpi->divideWork(gwan.size());
        offset = workDivs[0];
        numElements = workDivs[1] - workDivs[0];
        gwanSlice = gwan(Eigen::seq(workDivs[0], workDivs[1]));
      }

      // Each process writes to hdf5
      // The format is ((startRow,startCol),(numRows,numCols)).write(data)
      // Because it's a vector (1 row) all processes write to row=0,
      // col=startPoint
      // with nRows = 1, nCols = number of items this process will write.
      dgwannier.select({0, offset}, {1, numElements}).write(gwanSlice);
    }
#else
    { // do not remove these braces, see above note.
      // case where mpi exists, but HDF5 was built serially

      if (matrixDistributed) {
        // so, in this case, we must reduce the tensor before writing it
        Eigen::Tensor<std::complex<double>, 5> gWannierRed(
            numWannier, numWannier, numModes, numPhBravaisVectors,
            numElBravaisVectors);
        gWannierRed.setZero();
        for (int irE = 0; irE < numElBravaisVectors; irE++) {
          if (aux.indicesAreLocal(irE, 0)) { // is local
            int irELocal = aux.global2Local(irE, 0);
            for (int irP = 0; irP < numPhBravaisVectors; irP++) {
              for (int nu = 0; nu < numModes; nu++) {
                for (int i = 0; i < numWannier; i++) {
                  for (int j = 0; j < numWannier; j++) {
                    gWannierRed(i, j, nu, irP, irE) +=
                        gWannier(i, j, nu, irP, irELocal);
                  }
                }
              }
            }
          }
        }
        mpi->allReduceSum(&gWannierRed);
        gWannier = gWannierRed;
      }

      if (mpi->mpiHead()) {
        // open the hdf5 file
        HighFive::File file(outFileName, HighFive::File::Overwrite);

        // flatten the tensor (tensor is not supported) and create the data set
        Eigen::VectorXcd gwan = Eigen::Map<Eigen::VectorXcd, Eigen::Unaligned>(
            gWannier.data(), gWannier.size());
        HighFive::DataSet dgwannier = file.createDataSet<std::complex<double>>(
            "/gWannier", HighFive::DataSpace::From(gwan));

        // write to hdf5
        dgwannier.write(gwan);
      }
    }
#endif

    // we write the small quantities only with MPI head
    if (mpi->mpiHead()) {

      HighFive::File file(outFileName, HighFive::File::ReadWrite);

      // write out the number of electrons and the spin
      HighFive::DataSet dnElectrons = file.createDataSet<int>(
          "/numElectrons", HighFive::DataSpace::From(numFilledWannier));
      HighFive::DataSet dnSpin = file.createDataSet<int>(
          "/numSpin", HighFive::DataSpace::From(numSpin));
      dnElectrons.write(numFilledWannier); // # of occupied wannier functions
      dnSpin.write(numSpin);

      HighFive::DataSet dnElBands = file.createDataSet<int>(
          "/numElBands", HighFive::DataSpace::From(numWannier));
      HighFive::DataSet dnModes = file.createDataSet<int>(
          "/numPhModes", HighFive::DataSpace::From(numModes));
      dnElBands.write(numWannier);
      dnModes.write(numModes);

      // write out the kMesh and qMesh
      HighFive::DataSet dkMesh =
          file.createDataSet<int>("/kMesh", HighFive::DataSpace::From(kMesh));
      HighFive::DataSet dqMesh =
          file.createDataSet<int>("/qMesh", HighFive::DataSpace::From(qMesh));
      dkMesh.write(kMesh);
      dqMesh.write(qMesh);

      // write bravais lattice vectors
      HighFive::DataSet dPhBravais = file.createDataSet<double>(
          "/phBravaisVectors", HighFive::DataSpace::From(phBravaisVectors));
      HighFive::DataSet dElBravais = file.createDataSet<double>(
          "/elBravaisVectors", HighFive::DataSpace::From(elBravaisVectors));
      dPhBravais.write(phBravaisVectors);
      dElBravais.write(elBravaisVectors);

      // write electron and phonon degeneracies
      HighFive::DataSet dPhDegeneracies = file.createDataSet<double>(
          "/phDegeneracies", HighFive::DataSpace::From(phDegeneracies));
      HighFive::DataSet dElDegeneracies = file.createDataSet<double>(
          "/elDegeneracies", HighFive::DataSpace::From(elDegeneracies));
      dPhDegeneracies.write(phDegeneracies);
      dElDegeneracies.write(elDegeneracies);
    }
  } catch (std::exception &error) {
    Error("Issue writing elph Wannier representation to hdf5.");
  }

#else // need a non-hdf5 write option

  std::string outFileName = "./" + phoebePrefixQE + ".phoebe.elph.dat";

  std::ofstream outfile(outFileName);
  if (not outfile.is_open()) {
    Error("Output file couldn't be opened");
  }
  outfile << numFilledWannier << " " << numSpin << "\n";
  outfile << kMesh << "\n";
  outfile << qMesh << "\n";
  outfile << phBravaisVectors.rows() << " " << phBravaisVectors.cols() << "\n";
  outfile << phBravaisVectors << "\n";
  outfile << phDegeneracies << "\n";
  outfile << elBravaisVectors.rows() << " " << elBravaisVectors.cols() << "\n";
  outfile << elBravaisVectors << "\n";
  outfile << elDegeneracies << "\n";
  outfile << "\n";
  for (auto x : gWannier.dimensions()) {
    outfile << x << " ";
  }
  outfile << "\n";

  outfile << std::setprecision(16);
  int numPhBands = gWannier.dimension(2);

  for (int iRank = 0; iRank < mpi->getSize(); iRank++) {

    for (int i5 = 0; i5 < elDegeneracies.size(); i5++) {
      int i5Local = aux.global2Local(i5, 0);
      if (i5Local >= 0) {
        for (int i4 = 0; i4 < phDegeneracies.size(); i4++) {
          for (int i3 = 0; i3 < numPhBands; i3++) {
            for (int i2 = 0; i2 < numWannier; i2++) {
              for (int i1 = 0; i1 < numWannier; i1++) {
                outfile << std::setw(22)
                        << gWannier(i1, i2, i3, i4, i5Local).real() << " "
                        << std::setw(22)
                        << gWannier(i1, i2, i3, i4, i5Local).imag() << "\n";
              }
            }
          }
        }
      }
      mpi->barrier();
    }
  }
#endif

  if (mpi->mpiHead()) {
    std::cout << "Done writing el-ph coupling to file.\n" << std::endl;
  }
}