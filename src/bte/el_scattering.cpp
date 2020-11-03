#include "el_scattering.h"

#include "constants.h"
#include "helper_el_scattering.h"
#include "io.h"
#include "mpiHelper.h"
#include "periodic_table.h"

ElScatteringMatrix::ElScatteringMatrix(Context &context_,
                                       StatisticsSweep &statisticsSweep_,
                                       BaseBandStructure &innerBandStructure_,
                                       BaseBandStructure &outerBandStructure_,
                                       PhononH0 &h0_,
                                       InteractionElPhWan *couplingElPhWan_)
    : ScatteringMatrix(context_, statisticsSweep_, innerBandStructure_,
                       outerBandStructure_),
      couplingElPhWan(couplingElPhWan_), h0(h0_) {

  doBoundary = false;
  boundaryLength = context.getBoundaryLength();
  if (!std::isnan(boundaryLength)) {
    if (boundaryLength > 0.) {
      doBoundary = true;
    }
  }

  highMemory = context.getScatteringMatrixInMemory();
}

ElScatteringMatrix::ElScatteringMatrix(const ElScatteringMatrix &that)
    : ScatteringMatrix(that), couplingElPhWan(that.couplingElPhWan),
      h0(that.h0), boundaryLength(that.boundaryLength),
      doBoundary(that.doBoundary) {}

ElScatteringMatrix &
ElScatteringMatrix::operator=(const ElScatteringMatrix &that) {
  ScatteringMatrix::operator=(that);
  if (this != &that) {
    couplingElPhWan = that.couplingElPhWan;
    h0 = that.h0;
    boundaryLength = that.boundaryLength;
    doBoundary = that.doBoundary;
  }
  return *this;
}

// 3 cases:
// theMatrix and linedith is passed: we compute and store in memory the scatt
//       matrix and the diagonal
// inPopulation+outPopulation is passed: we compute the action of the
//       scattering matrix on the in vector, returning outVec = sMatrix*vector
// only linewidth is passed: we compute only the linewidths
void ElScatteringMatrix::builder(VectorBTE *linewidth,
                                 std::vector<VectorBTE> &inPopulations,
                                 std::vector<VectorBTE> &outPopulations) {
  const double energyCutoff = 0.001 / ryToCmm1; // discard states with small
  // phonon energies (smaller than 0.001 cm^-1

  int switchCase = 0;
  if (theMatrix.rows() != 0 && linewidth != nullptr &&
      inPopulations.size() == 0 && outPopulations.size() == 0) {
    switchCase = 0;
  } else if (theMatrix.rows() == 0 && linewidth == nullptr &&
             inPopulations.size() != 0 && outPopulations.size() != 0) {
    switchCase = 1;
  } else if (theMatrix.rows() == 0 && linewidth != nullptr &&
             inPopulations.size() == 0 && outPopulations.size() == 0) {
    switchCase = 2;
  } else {
    Error e("builder3Ph found a non-supported case");
  }

  if ((linewidth != nullptr) && (linewidth->dimensionality != 1)) {
    Error e("The linewidths shouldn't have dimensionality");
  }

  auto particle = outerBandStructure.getParticle();

  int numCalcs = statisticsSweep.getNumCalcs();

  // note: innerNumFullPoints is the number of points in the full grid
  // may be larger than innerNumPoints, when we use ActiveBandStructure
  int innerNumFullPoints = innerBandStructure.getNumPoints(true);

  // precompute Fermi-Dirac populations
  VectorBTE outerFermi(statisticsSweep, outerBandStructure, 1);
#pragma omp parallel for
  for (int is : mpi->divideWorkIter(outerBandStructure.getNumStates())) {
    double energy = outerBandStructure.getEnergy(is);
    for (int iCalc = 0; iCalc < statisticsSweep.getNumCalcs(); iCalc++) {
      double temperature = statisticsSweep.getCalcStatistics(iCalc).temperature;
      double chemicalPotential =
          statisticsSweep.getCalcStatistics(iCalc).chemicalPotential;
      outerFermi.data(iCalc, is) =
          particle.getPopulation(energy, temperature, chemicalPotential);
    }
  }
  mpi->allReduceSum(&outerFermi.data);

  VectorBTE innerFermi(statisticsSweep, innerBandStructure, 1);
  if (&innerBandStructure == &outerBandStructure) {
    innerFermi = outerFermi;
  } else {
#pragma omp parallel for
    for (int is : mpi->divideWorkIter(innerBandStructure.getNumStates())) {
      double energy = innerBandStructure.getEnergy(is);
      for (int iCalc = 0; iCalc < statisticsSweep.getNumCalcs(); iCalc++) {
        double temperature =
            statisticsSweep.getCalcStatistics(iCalc).temperature;
        double chemicalPotential =
            statisticsSweep.getCalcStatistics(iCalc).chemicalPotential;
        innerFermi.data(iCalc, is) =
            particle.getPopulation(energy, temperature, chemicalPotential);
      }
    }
    mpi->allReduceSum(&innerFermi.data);
  }

  if (smearing->getType() == DeltaFunction::tetrahedron) {
    Error e("Tetrahedron method not supported by electron scattering");
    // that's because it doesn't work with the window the way it's implemented,
    // and we will almost always have a window for electrons
  }

  bool rowMajor = true;
  std::vector<std::tuple<std::vector<long>, long>> kPairIterator =
      getIteratorWavevectorPairs(switchCase, rowMajor);

  HelperElScattering pointHelper(outerBandStructure, innerBandStructure,
                                 statisticsSweep, smearing->getType(), h0);

  LoopPrint loopPrint("computing scattering matrix", "k-points",
                      kPairIterator.size());

  for (auto t1 : kPairIterator) {
    auto ik2Indexes = std::get<0>(t1);
    auto ik1Irr = std::get<1>(t1);

    WavevectorIndex ik1IrrIndex(ik1Irr);
    Eigen::Vector3d k1IrrC = outerBandStructure.getWavevector(ik1IrrIndex);
    Eigen::VectorXd state1Energies =
        outerBandStructure.getEnergies(ik1IrrIndex);
    int nb1 = state1Energies.size();
    Eigen::MatrixXd v1sIrr = outerBandStructure.getGroupVelocities(ik1IrrIndex);
    Eigen::MatrixXcd eigvec1 = outerBandStructure.getEigenvectors(ik1IrrIndex);

    //    auto rotations =
    //    innerBandStructure.getPoints().getRotationsStar(ik1Irr); for
    //    (Eigen::Matrix3d rotation : rotations) {
    //      Eigen::Vector3d k1C = rotation * k1IrrC;
    //      Eigen::MatrixXd v1s = v1sIrr;
    //      for (int ib1 = 0; ib1 < nb1; ib1++) {
    //        Eigen::Vector3d thisV1Irr = v1sIrr.row(ib1);
    //        Eigen::Vector3d thisV1 = rotation * thisV1Irr;
    //        v1s.row(ib1) = thisV1;
    //      }
    Eigen::MatrixXd v1s = v1sIrr;
    Eigen::Vector3d k1C = k1IrrC;

    loopPrint.update();
    pointHelper.prepare(k1C, ik2Indexes);
    int nk2 = ik2Indexes.size();
    std::vector<Eigen::Vector3d> allQ3C(nk2);
    std::vector<Eigen::VectorXd> allStates3Energies(nk2);
    std::vector<int> allNb3(nk2);
    std::vector<Eigen::MatrixXcd> allEigvecs3(nk2);
    std::vector<Eigen::MatrixXd> allV3s(nk2);
    std::vector<Eigen::MatrixXd> allBose3Data(nk2);

    std::vector<Eigen::Vector3d> allK2C(nk2);
    std::vector<Eigen::MatrixXcd> allEigvecs2(nk2);
    std::vector<Eigen::VectorXd> allState2Energies(nk2);
    std::vector<Eigen::MatrixXd> allV2s(nk2);

    int ik2Counter = -1;
    for (long ik2 : ik2Indexes) {
      ik2Counter++;
      WavevectorIndex ik2Index(ik2);
      allK2C[ik2Counter] = innerBandStructure.getWavevector(ik2Index);
      allState2Energies[ik2Counter] = innerBandStructure.getEnergies(ik2Index);
      allV2s[ik2Counter] = innerBandStructure.getGroupVelocities(ik2Index);
      allEigvecs2[ik2Counter] = innerBandStructure.getEigenvectors(ik2Index);
      auto t2 = pointHelper.get(k1C, ik2);
      allQ3C[ik2Counter] = std::get<0>(t2);
      allStates3Energies[ik2Counter] = std::get<1>(t2);
      allNb3[ik2Counter] = std::get<2>(t2);
      allEigvecs3[ik2Counter] = std::get<3>(t2);
      allV3s[ik2Counter] = std::get<4>(t2);
      allBose3Data[ik2Counter] = std::get<5>(t2);
    }

    couplingElPhWan->calcCouplingSquared(eigvec1, allEigvecs2, allEigvecs3, k1C,
                                         allK2C, allQ3C);

    ik2Counter = -1;
    for (auto ik2Irr : ik2Indexes) {
      ik2Counter++;
      auto coupling = couplingElPhWan->getCouplingSquared(ik2Counter);

      Eigen::VectorXd state2Energies = allState2Energies[ik2Counter];
      int nb2 = state2Energies.size();
      Eigen::MatrixXd v2s = allV2s[ik2Counter];

      int nb3 = allNb3[ik2Counter];
      Eigen::VectorXd state3Energies = allStates3Energies[ik2Counter];
      Eigen::VectorXd bose3Data = allBose3Data[ik2Counter];
      Eigen::MatrixXd v3s = allV3s[ik2Counter];

      int ib2, ib3;
//#pragma omp parallel for private(ib2, ib3) collapse(3)
      for (int ib1 = 0; ib1 < nb1; ib1++) {
        for (ib2 = 0; ib2 < nb2; ib2++) {
          for (ib3 = 0; ib3 < nb3; ib3++) {
            double en1 = state1Energies(ib1);
            double en2 = state2Energies(ib2);
            double en3 = state3Energies(ib3);
            // remove small divergent phonon energies
            if (en3 < energyCutoff) {
              continue;
            }

            int ind1 = outerBandStructure.getIndex(WavevectorIndex(ik1Irr),
                                                   BandIndex(ib1));
            int ind2 = innerBandStructure.getIndex(WavevectorIndex(ik2Irr),
                                                   BandIndex(ib2));

            if (switchCase == 0) {
              // note: above we are parallelizing over wavevectors.
              // (for convenience of computing coupling3ph)
              // Not the same way as Matrix() is parallelized.
              // here we check that we don't duplicate efforts
              if (!theMatrix.indecesAreLocal(ind1, ind2)) {
                continue;
              }
            }

            double delta1, delta2;
            if (smearing->getType() == DeltaFunction::gaussian) {
              delta1 = smearing->getSmearing(en1 - en2 + en3);
              delta2 = smearing->getSmearing(en1 - en2 - en3);
            } else {
              // Eigen::Vector3d smear = v1s.row(ib1s) - v2s.row(ib2);
              Eigen::Vector3d smear = v3s.row(ib3);
              delta1 = smearing->getSmearing(en1 - en2 + en3, smear);
              delta2 = smearing->getSmearing(en1 - en2 - en3, smear);
            }

            if (delta1 < 0. && delta2 < 0.)
              continue;

            // loop on temperature
            for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
              double fermi1 = outerFermi.data(iCalc, ind1);
              double fermi2 = innerFermi.data(iCalc, ind2);
              double bose3 = bose3Data(iCalc, ib3);

              // Calculate transition probability W+
              double rate =
                  coupling(ib1, ib2, ib3) *
                  ((fermi2 + bose3) * delta1 + (1. - fermi2 + bose3) * delta2) /
                  innerNumFullPoints * pi / en3;

              double rateOffDiag = coupling(ib1, ib2, ib3) *
                                   (fermi1 * (1. - fermi2) * bose3 * delta1 +
                                    fermi2 * (1. - fermi1) * bose3 * delta2) /
                                   innerNumFullPoints * pi / en3;

              if (switchCase == 0) {
                // case of matrix construction
                // we build the scattering matrix S
                //                  for (int i = 0; i < dimensionality_; i++) {
                //                    for (int j = 0; j < dimensionality_; j++)
                //                    {
                //                      addMatrixElement(rotation(i, j) *
                //                      rateOffDiag, ind1, ind2,
                //                                       i, j);
                //                    }
                //                  }
                theMatrix(ind1, ind2) += rateOffDiag;
                linewidth->operator()(iCalc, 0, ind1) += rate;
              } else if (switchCase == 1) {
                // case of matrix-vector multiplication
                // we build the scattering matrix A = S*n(n+1)

                for (unsigned int iVec = 0; iVec < inPopulations.size();
                     iVec++) {
                  //                    for (int i = 0; i < dimensionality_;
                  //                    i++) {
                  //                      for (int j = 0; j < dimensionality_;
                  //                      j++) {
                  //                        outPopulations[iVec](iCalc, i, ind1)
                  //                        +=
                  //                            rateOffDiag * rotation(i, j) *
                  //                                inPopulations[iVec](iCalc,
                  //                                j, ind1);
                  //                        outPopulations[iVec](iCalc, i, ind1)
                  //                        +=
                  //                            rate * rotation(i, j) *
                  //                                inPopulations[iVec](iCalc,
                  //                                j, ind1);
                  //                      }
                  //                    }
                  for (int i : {0,1,2}) {
                    outPopulations[iVec](iCalc, i, ind1) +=
                        rateOffDiag * inPopulations[iVec](iCalc, i, ind2);
                    outPopulations[iVec](iCalc, i, ind1) +=
                        rate * inPopulations[iVec](iCalc, i, ind1);
                  }
                }
              } else {
                // case of linewidth construction
                linewidth->operator()(iCalc, 0, ind1) += rate;
              }
            }
          }
        }
      }
      //      }
    }
  }
  if (switchCase == 1) {
    for (unsigned int iVec = 0; iVec < inPopulations.size(); iVec++) {
      mpi->allReduceSum(&outPopulations[iVec].data);
    }
  } else {
    mpi->allReduceSum(&linewidth->data);
  }
  // I prefer to close loopPrint after the MPI barrier: all MPI are synced here
  loopPrint.close();

  // Add boundary scattering

  if (doBoundary) {
#pragma omp parallel for
    for (int is1 = 0; is1 < numStates; is1++) {
      StateIndex is1Ind(is1);
      double energy = outerBandStructure.getEnergy(is1Ind);
      auto vel = outerBandStructure.getGroupVelocity(is1Ind);
      for (int iCalc = 0; iCalc < statisticsSweep.getNumCalcs(); iCalc++) {
        double temperature =
            statisticsSweep.getCalcStatistics(iCalc).temperature;
        // n(n+1)
        double termPop = particle.getPopPopPm1(energy, temperature);
        double rate = vel.squaredNorm() / boundaryLength * termPop;

        switch (switchCase) {
        case (0):
          // case of matrix construction
          linewidth->operator()(iCalc, 0, is1) += rate;
          break;
        case (1):
          // case of matrix-vector multiplication
          for (unsigned int iVec = 0; iVec < inPopulations.size(); iVec++) {
            for (int i = 0; i < dimensionality_; i++) {
              outPopulations[iVec](iCalc, i, is1) +=
                  rate * inPopulations[iVec](iCalc, i, is1);
            }
          }
          break;
        case (2):
          // case of linewidth construction
          linewidth->operator()(iCalc, 0, is1) += rate;
          break;
        }
      }
    }
  }

  // we place the linewidths back in the diagonal of the scattering matrix
  // this because we may need an MPI_allreduce on the linewidths
  if (switchCase == 0) { // case of matrix construction
    long iCalc = 0;
#pragma omp parallel for
    for (long is = 0; is < outerBandStructure.getNumStates(); is++) {
      // TODO: check if this is the right assignment
//      for (int alfa = 0; alfa < dimensionality_; alfa++) {
//        for (int beta = 0; beta < dimensionality_; beta++) {
//          setMatrixElement(linewidth->operator()(iCalc, 0, is), is, is, alfa,
//                           beta);
//        }
//      }
      theMatrix(is,is) = linewidth->operator()(iCalc, 0, is);
    }
  }
}

double ElScatteringMatrix::getMatrixElement(const int &m, const int &n,
                                            const int &alfa, const int &beta) {
  int i = m * dimensionality_ + alfa;
  int j = n * dimensionality_ + beta;
  return theMatrix(i, j);
}

void ElScatteringMatrix::setMatrixElement(const double &x, const int &m,
                                          const int &n, const int &alfa,
                                          const int &beta) {
  int i = compress2Indeces(m, alfa, numStates, dimensionality_);
  int j = compress2Indeces(n, beta, numStates, dimensionality_);
  theMatrix(i, j) = x;
}

void ElScatteringMatrix::addMatrixElement(const double &x, const int &m,
                                          const int &n, const int &alfa,
                                          const int &beta) {
  int i = compress2Indeces(m, alfa, numStates, dimensionality_);
  int j = compress2Indeces(n, beta, numStates, dimensionality_);
  theMatrix(i, j) += x;
}
//
// VectorBTE ElScatteringMatrix::dot(VectorBTE &inPopulation) {
//  if (highMemory) {
//    VectorBTE outPopulation(statisticsSweep, outerBandStructure,
//                            inPopulation.dimensionality);
//    outPopulation.data.setZero();
//    // note: we are assuming that ScatteringMatrix has numCalcs = 1
//    for (auto t0 : theMatrix.getAllLocalStates()) {
//      auto i = std::get<0>(t0);
//      auto j = std::get<1>(t0);
//
//      auto t1 = decompress2Indeces(i, numStates, dimensionality_);
//      auto m = std::get<0>(t1);
//      auto alfa = std::get<1>(t1);
//      auto t2 = decompress2Indeces(j, numStates, dimensionality_);
//      auto n = std::get<0>(t2);
//      auto beta = std::get<1>(t2);
//      outPopulation(0, alfa, m) +=
//          getMatrixElement(m, n, alfa, beta) * inPopulation(0, beta, n);
//    }
//    mpi->allReduceSum(&outPopulation.data);
//    return outPopulation;
//  } else {
//    VectorBTE outPopulation(statisticsSweep, outerBandStructure,
//                            inPopulation.dimensionality);
//    builder(nullptr, &inPopulation, &outPopulation);
//    return outPopulation;
//  }
//}
