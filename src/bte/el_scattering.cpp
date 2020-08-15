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
                                       InteractionElPhWan *couplingElPhWan_,
                                       PhononH0 *h0_)
    : ScatteringMatrix(context_, statisticsSweep_, innerBandStructure_,
                       outerBandStructure_),
      couplingElPhWan(couplingElPhWan_),
      h0(h0_) {
  if (&innerBandStructure != &outerBandStructure && h0 == nullptr) {
    Error e("PhScatteringMatrix needs h0 for incommensurate grids");
  }

  doBoundary = false;
  boundaryLength = context.getBoundaryLength();
  if (!std::isnan(boundaryLength)) {
    if (boundaryLength > 0.) {
      doBoundary = true;
    }
  }
}

ElScatteringMatrix::ElScatteringMatrix(const ElScatteringMatrix &that)
    : ScatteringMatrix(that),
      couplingElPhWan(that.couplingElPhWan),
      h0(that.h0),
      boundaryLength(that.boundaryLength),
      doBoundary(that.doBoundary) {}

ElScatteringMatrix &ElScatteringMatrix::operator=(
    const ElScatteringMatrix &that) {
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
                                 std::vector<VectorBTE> inPopulations,
                                 std::vector<VectorBTE> outPopulations) {
  const double energyCutoff = 0.001 / ryToCmm1;  // discard states with small
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
      outerFermi.data(iCalc, is) = particle.getPopulation(energy, temperature);
    }
  }
  mpi->allReduceSum(&outerFermi.data);
  VectorBTE innerFermi(statisticsSweep, outerBandStructure, 1);
  if (&innerBandStructure == &outerBandStructure) {
    innerFermi = outerFermi;
  } else {
#pragma omp parallel for
    for (int is : mpi->divideWorkIter(innerBandStructure.getNumStates())) {
      double energy = innerBandStructure.getEnergy(is);
      for (int iCalc = 0; iCalc < statisticsSweep.getNumCalcs(); iCalc++) {
        double temperature =
            statisticsSweep.getCalcStatistics(iCalc).temperature;
        innerFermi.data(iCalc, is) =
            particle.getPopulation(energy, temperature);
      }
    }
    mpi->allReduceSum(&innerFermi.data);
  }

  std::vector<std::tuple<std::vector<long>, long>> kPairIterator =
      getIteratorWavevectorPairs(switchCase);

  HelperElScattering pointHelper(innerBandStructure, outerBandStructure,
                                 statisticsSweep, smearing->getType(), h0);

  LoopPrint loopPrint("computing scattering matrix", "k-points",
                      kPairIterator.size());

  for (auto t1 : kPairIterator) {
    auto ik1Indexes = std::get<0>(t1);
    auto ik2Irr = std::get<1>(t1);

    WavevectorIndex ik2IrrIndex(ik2Irr);
    Eigen::Vector3d k2IrrC = innerBandStructure.getWavevector(ik2IrrIndex);
    Eigen::VectorXd state2Energies =
        innerBandStructure.getEnergies(ik2IrrIndex);
    int nb2 = state2Energies.size();
    Eigen::MatrixXd v2sIrr = innerBandStructure.getGroupVelocities(ik2IrrIndex);
    Eigen::MatrixXcd eigvec2 = innerBandStructure.getEigenvectors(ik2IrrIndex);

    auto rotations = innerBandStructure.getPoints().getRotationsStar(ik2Irr);
    for (Eigen::Matrix3d rotation : rotations) {
      Eigen::Vector3d k2C = rotation * k2IrrC;
      Eigen::MatrixXd v2s = v2sIrr;
      for (int ib2 = 0; ib2 < nb2; ib2++) {
        Eigen::Vector3d thisV2Irr = v2sIrr.row(ib2);
        Eigen::Vector3d thisV2 = rotation * thisV2Irr;
        v2s.row(ib2) = thisV2;
      }

      loopPrint.update();
      pointHelper.prepare(ik1Indexes, k2C);

      std::vector<Eigen::Vector3d> allQ3C(ik1Indexes.size());
      std::vector<Eigen::VectorXd> allStates3Energies(ik1Indexes.size());
      std::vector<int> allNb3(ik1Indexes.size());
      std::vector<Eigen::MatrixXcd> allEigvecs3(ik1Indexes.size());
      std::vector<Eigen::MatrixXd> allV3s(ik1Indexes.size());
      std::vector<Eigen::MatrixXd> allBose3Data(ik1Indexes.size());

      std::vector<Eigen::Vector3d> allK1C(ik1Indexes.size());
      std::vector<Eigen::MatrixXcd> allEigvecs1(ik1Indexes.size());
      std::vector<Eigen::VectorXd> allState1Energies(ik1Indexes.size());
      std::vector<Eigen::MatrixXd> allV1s(ik1Indexes.size());

      for (long ik1 : ik1Indexes) {
        WavevectorIndex ik1Index(ik1);
        allK1C[ik1] = outerBandStructure.getWavevector(ik1Index);
        allState1Energies[ik1] = outerBandStructure.getEnergies(ik1Index);
        allV1s[ik1] = outerBandStructure.getGroupVelocities(ik1Index);
        allEigvecs1[ik1] = outerBandStructure.getEigenvectors(ik1Index);

        //        auto [q3C, state3Energies, nb3, eigvecs3, v3s, bose3Data] =
        //            pointHelper.get(ik1, k2C);
        //        allQ3C[ik1] = q3C;
        //        allStates3Energies[ik1] = state3Energies;
        //        allNb3[ik1] = nb3;
        //        allEigvecs3[ik1] = eigvecs3;
        //        allV3s[ik1] = v3s;
        //        allBose3Data[ik1] = bose3Data;
        auto t2 = pointHelper.get(ik1, k2C);
        allQ3C[ik1] = std::get<0>(t2);
        allStates3Energies[ik1] = std::get<1>(t2);
        allNb3[ik1] = std::get<2>(t2);
        allEigvecs3[ik1] = std::get<3>(t2);
        allV3s[ik1] = std::get<4>(t2);
        allBose3Data[ik1] = std::get<5>(t2);
      }

      couplingElPhWan->calcCouplingSquared(eigvec2, allEigvecs1, allEigvecs3,
                                           k2C, allK1C, allQ3C);

      for (auto ik1Irr : ik1Indexes) {
        auto coupling = couplingElPhWan->getCouplingSquared(ik1Irr);

        Eigen::VectorXd state1Energies = allState1Energies[ik1Irr];
        int nb1 = state1Energies.size();
        Eigen::MatrixXd v1s = allV1s[ik1Irr];

        int nb3 = allNb3[ik1Irr];
        Eigen::VectorXd state3Energies = allStates3Energies[ik1Irr];
        Eigen::VectorXd bose3Data = allBose3Data[ik1Irr];
        Eigen::MatrixXd v3s = allV3s[ik1Irr];

        int ib2, ib3;
#pragma omp parallel for private(ib2, ib3) collapse(3)
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
              //
              //              if (switchCase == 0) {
              //                // note: above we are parallelizing over
              //                wavevectors.
              //                // (for convenience of computing coupling3ph)
              //                // Not the same way as Matrix() is parallelized.
              //                // here we check that we don't duplicate efforts
              //                if (!theMatrix.indecesAreLocal(ind1, ind2)) {
              //                  continue;
              //                }
              //              }

              double delta1, delta2;
              if (smearing->getType() == DeltaFunction::gaussian) {
                delta1 = smearing->getSmearing(en1 - en2 + en3);
                delta2 = smearing->getSmearing(en1 - en2 - en3);
              } else if (smearing->getType() ==
                         DeltaFunction::adaptiveGaussian) {
                Eigen::Vector3d va = v2s.row(ib2) - v3s.row(ib3);
                Eigen::Vector3d vb = v2s.row(ib2) + v3s.row(ib3);
                delta1 = smearing->getSmearing(en1 - en2 + en3, va);
                delta2 = smearing->getSmearing(en1 - en2 - en3, vb);
              } else {
                delta1 = smearing->getSmearing(en3 + en1, ik2Irr, ib2);
                delta2 = smearing->getSmearing(en3 - en1, ik2Irr, ib2);
              }

              if (delta1 < 0. && delta2 < 0.) continue;

              // loop on temperature
              for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
                double fermi1 = outerFermi.data(iCalc, ind1);
                double fermi2 = innerFermi.data(iCalc, ind2);
                double bose3 = bose3Data(iCalc, ib3);

                // Calculate transition probability W+
                double rate = coupling(ib1, ib2, ib3) *
                              ((fermi2 + bose3) * delta1 +
                               (1. - fermi2 + bose3) * delta2) /
                              innerNumFullPoints * pi * 0.5 / en3;

                double rateOffDiag = coupling(ib1, ib2, ib3) *
                                     (fermi1 * (1. - fermi2) * bose3 * delta1 +
                                      fermi2 * (1. - fermi1) * bose3 * delta2) /
                                     innerNumFullPoints * pi * 0.5 / en3;

                if (switchCase == 0) {
                  // case of matrix construction
                  // we build the scattering matrix S
                  for (int i = 0; i < dimensionality_; i++) {
                    for (int j = 0; j < dimensionality_; j++) {
                      addMatrixElement(rotation(i, j) * rateOffDiag, ind1, ind2,
                                       i, j);
                    }
                  }
                  linewidth->operator()(iCalc, 0, ind1) += rate;
                } else if (switchCase == 1) {
                  // case of matrix-vector multiplication
                  // we build the scattering matrix A = S*n(n+1)

                  for (unsigned int iVec=0; iVec<inPopulations.size(); iVec++) {
                    for (int i = 0; i < dimensionality_; i++) {
                      for (int j = 0; j < dimensionality_; j++) {
                        outPopulations[iVec](iCalc, i, ind1) +=
                            rateOffDiag * rotation(i, j) *
                            inPopulations[iVec](iCalc, j, ind1);
                        outPopulations[iVec](iCalc, i, ind1) +=
                            rate * rotation(i, j) *
                            inPopulations[iVec](iCalc, j, ind1);
                      }
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
      }
    }
  }

  if (switchCase == 1) {
    for (unsigned int iVec=0; iVec<inPopulations.size(); iVec++) {
      mpi->allReduceSum(&outPopulations[iVec]->data);
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
            for (unsigned int iVec=0; iVec<inPopulations.size(); iVec++) {
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

  // some phonons like acoustic modes at the gamma, with omega = 0,
  // might have zero frequencies, and infinite populations. We set those
  // matrix elements to zero.
  if (switchCase == 0) {
    // case of matrix construction
    for (auto is1 : excludeIndeces) {
      linewidth->data.col(is1).setZero();
      for (auto is2 : excludeIndeces) {
        for (int alfa = 0; alfa < dimensionality_; alfa++) {
          for (int beta = 0; beta < dimensionality_; beta++) {
            setMatrixElement(0., is1, is2, alfa, beta);
          }
        }
      }
    }

  } else if (switchCase == 1) {
    // case of matrix-vector multiplication
    for (unsigned int iVec=0; iVec<inPopulations.size(); iVec++) {
      for (auto is1 : excludeIndeces) {
        outPopulations[iVec].data.col(is1).setZero();
      }
    }

  } else if (switchCase == 2) {
    // case of linewidth construction
    for (auto is1 : excludeIndeces) {
      linewidth->data.col(is1).setZero();
    }
  }

  // we place the linewidths back in the diagonal of the scattering matrix
  // this because we may need an MPI_allreduce on the linewidths
  if (switchCase == 0) {  // case of matrix construction
    long iCalc = 0;
#pragma omp parallel for
    for (long is = 0; is < outerBandStructure.getNumStates(); is++) {
      // TODO: check if this is the right assignment
      for (int alfa = 0; alfa < dimensionality_; alfa++) {
        for (int beta = 0; beta < dimensionality_; beta++) {
          setMatrixElement(linewidth->operator()(iCalc, 0, is), is, is, alfa,
                           beta);
        }
      }
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

VectorBTE ElScatteringMatrix::dot(VectorBTE &inPopulation) {
  if (highMemory) {
    VectorBTE outPopulation(statisticsSweep, outerBandStructure,
                            inPopulation.dimensionality);
    outPopulation.data.setZero();
    // note: we are assuming that ScatteringMatrix has numCalcs = 1
    for (auto t0 : theMatrix.getAllLocalStates()) {
      auto i = std::get<0>(t0);
      auto j = std::get<1>(t0);

      auto t1 = decompress2Indeces(i, numStates, dimensionality_);
      auto m = std::get<0>(t1);
      auto alfa = std::get<1>(t1);
      auto t2 = decompress2Indeces(j, numStates, dimensionality_);
      auto n = std::get<0>(t2);
      auto beta = std::get<1>(t2);
      outPopulation(0, alfa, m) +=
          getMatrixElement(m, n, alfa, beta) * inPopulation(0, beta, n);
    }
    mpi->allReduceSum(&outPopulation.data);
    return outPopulation;
  } else {
    VectorBTE outPopulation(statisticsSweep, outerBandStructure,
                            inPopulation.dimensionality);
    builder(nullptr, &inPopulation, &outPopulation);
    return outPopulation;
  }
}
