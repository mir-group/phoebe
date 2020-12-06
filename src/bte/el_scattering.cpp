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

  long numCalcs = statisticsSweep.getNumCalcs();

  // note: innerNumFullPoints is the number of points in the full grid
  // may be larger than innerNumPoints, when we use ActiveBandStructure
  double norm = 1. / context.getKMesh().prod();

  // precompute Fermi-Dirac populations
  long numOuterIrrStates = outerBandStructure.irrStateIterator().size();
  Eigen::MatrixXd outerFermi(numCalcs, numOuterIrrStates);
  outerFermi.setZero();
#pragma omp parallel for
  for (long ibte : mpi->divideWorkIter(numOuterIrrStates)) {
    BteIndex ibteIdx = BteIndex(ibte);
    StateIndex isIdx = outerBandStructure.bteToState(ibteIdx);
    double energy = outerBandStructure.getEnergy(isIdx);
    for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      double temp = calcStat.temperature;
      double chemPot = calcStat.chemicalPotential;
      outerFermi(iCalc, ibte) = particle.getPopulation(energy, temp, chemPot);
    }
  }
  mpi->allReduceSum(&outerFermi);

  long numInnerIrrStates = innerBandStructure.irrStateIterator().size();
  Eigen::MatrixXd innerFermi(numCalcs, numInnerIrrStates);
  innerFermi.setZero();
#pragma omp parallel for
  for (long ibte : mpi->divideWorkIter(numInnerIrrStates)) {
    BteIndex ibteIdx = BteIndex(ibte);
    StateIndex isIdx = innerBandStructure.bteToState(ibteIdx);
    double energy = innerBandStructure.getEnergy(isIdx);
    for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      double temp = calcStat.temperature;
      double chemPot = calcStat.chemicalPotential;
      innerFermi(iCalc, ibte) = particle.getPopulation(energy, temp, chemPot);
    }
  }
  mpi->allReduceSum(&innerFermi);

  if (smearing->getType() == DeltaFunction::tetrahedron) {
    Error e("Tetrahedron method not supported by electron scattering");
    // that's because it doesn't work with the window the way it's implemented,
    // and we will almost always have a window for electrons
  }

  bool rowMajor = true;
  std::vector<std::tuple<std::vector<long>, long>> kPairIterator =
      getIteratorWavevectorPairs(switchCase, rowMajor);

  HelperElScattering pointHelper(innerBandStructure, outerBandStructure,
                                 statisticsSweep, smearing->getType(), h0);

  LoopPrint loopPrint("computing scattering matrix", "k-points",
                      kPairIterator.size());

  for (auto t1 : kPairIterator) {
    auto ik2Indexes = std::get<0>(t1);
    long ik1 = std::get<1>(t1);
    WavevectorIndex ik1Idx(ik1);

    Eigen::Vector3d k1C = outerBandStructure.getWavevector(ik1Idx);
    Eigen::VectorXd state1Energies = outerBandStructure.getEnergies(ik1Idx);
    long nb1 = state1Energies.size();
    Eigen::MatrixXd v1s = outerBandStructure.getGroupVelocities(ik1Idx);
    Eigen::MatrixXcd eigvec1 = outerBandStructure.getEigenvectors(ik1Idx);

    loopPrint.update();
    pointHelper.prepare(k1C, ik2Indexes);
    long nk2 = ik2Indexes.size();
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

    long ik2Counter = -1;
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

    bool withSymmetries = context.getUseSymmetries();

    ik2Counter = -1;
    for (auto ik2 : ik2Indexes) {
      ik2Counter++;
      Eigen::Tensor<double,3> coupling = couplingElPhWan->getCouplingSquared(ik2Counter);

      Eigen::VectorXd state2Energies = allState2Energies[ik2Counter];
      long nb2 = state2Energies.size();
      Eigen::MatrixXd v2s = allV2s[ik2Counter];

      long nb3 = allNb3[ik2Counter];
      Eigen::VectorXd state3Energies = allStates3Energies[ik2Counter];
      Eigen::MatrixXd bose3Data = allBose3Data[ik2Counter];
      Eigen::MatrixXd v3s = allV3s[ik2Counter];

      auto t3 = innerBandStructure.getRotationToIrreducible(
          allK2C[ik2Counter], Points::cartesianCoords);
      long ik2Irr = std::get<0>(t3);
      Eigen::Matrix3d rotation = std::get<1>(t3);

      WavevectorIndex ik2Idx(ik2);
      WavevectorIndex ik2IrrIdx(ik2Irr);

//#pragma omp parallel for
      for (long ib1 = 0; ib1 < nb1; ib1++) {
        double en1 = state1Energies(ib1);
        for (long ib2 = 0; ib2 < nb2; ib2++) {
          double en2 = state2Energies(ib2);
          for (long ib3 = 0; ib3 < nb3; ib3++) {
            double en3 = state3Energies(ib3);
            // remove small divergent phonon energies
            if (en3 < energyCutoff) {
              continue;
            }

            long is1 = outerBandStructure.getIndex(ik1Idx, BandIndex(ib1));
            long is2 = innerBandStructure.getIndex(ik2Idx, BandIndex(ib2));
            long is2Irr = innerBandStructure.getIndex(ik2IrrIdx,
                                                      BandIndex(ib2));
            StateIndex is1Idx(is1);
            StateIndex is2Idx(is2);
            StateIndex is2IrrIdx(is2Irr);
            BteIndex ind1Idx = outerBandStructure.stateToBte(is1Idx);
            BteIndex ind2Idx = innerBandStructure.stateToBte(is2IrrIdx);
            long ibte1 = ind1Idx.get();
            long ibte2 = ind2Idx.get();

            double delta1, delta2;
            if (smearing->getType() == DeltaFunction::gaussian) {
              delta1 = smearing->getSmearing(en1 - en2 + en3);
              delta2 = smearing->getSmearing(en1 - en2 - en3);
            } else if (smearing->getType() == DeltaFunction::adaptiveGaussian) {
              // Eigen::Vector3d smear = v1s.row(ib1s) - v2s.row(ib2);
              Eigen::Vector3d smear = v3s.row(ib3);
              delta1 = smearing->getSmearing(en1 - en2 + en3, smear);
              delta2 = smearing->getSmearing(en1 - en2 - en3, smear);
            } else {
              delta1 = smearing->getSmearing(en3 + en1, is2Idx);
              delta2 = smearing->getSmearing(en3 - en1, is2Idx);
            }

            if (delta1 < 0. && delta2 < 0.) {
              continue;
            }

            // loop on temperature
            for (long iCalc = 0; iCalc < numCalcs; iCalc++) {

              double fermi1 = outerFermi(iCalc, ibte1);
              double fermi2 = innerFermi(iCalc, ibte2);
              double bose3 = bose3Data(iCalc, ib3);

              // Calculate transition probability W+
              double rate =
                  coupling(ib1, ib2, ib3) *
                  ((fermi2 + bose3) * delta1 + (1. - fermi2 + bose3) * delta2) *
                  norm * pow(twoPi, 4) / en3;

              double rateOffDiag = coupling(ib1, ib2, ib3) *
                                   (fermi1 * (1. - fermi2) * bose3 * delta1 +
                                    fermi2 * (1. - fermi1) * bose3 * delta2) *
                                   norm * pow(twoPi, 4) / en3;

              if (switchCase == 0) {

                if (withSymmetries) {
                  for (long i : {0, 1, 2}) {
                    for (long j : {0, 1, 2}) {
                      auto iIndex = CartIndex(i);
                      auto jIndex = CartIndex(j);
                      long iMat1 = getSMatrixIndex(ind1Idx, iIndex);
                      long iMat2 = getSMatrixIndex(ind2Idx, jIndex);
                      if (theMatrix.indecesAreLocal(iMat1, iMat2)) {
                        if (iMat1 == 0 && iMat2 == 0) {
                          linewidth->operator()(iCalc, 0, ibte1) += rate;
                        }
                        theMatrix(iMat1, iMat2) += rotation(i, j) * rateOffDiag;
                      }
                    }
                  }
                } else {
                  if (theMatrix.indecesAreLocal(ibte1,ibte2)) {
                    linewidth->operator()(iCalc, 0, ibte1) += rate;
                    theMatrix(ibte1, ibte2) += rateOffDiag;
                  }
                }
              } else if (switchCase == 1) {
                // case of matrix-vector multiplication
                // we build the scattering matrix A = S*n(n+1)

                for (unsigned int iVec = 0; iVec < inPopulations.size();
                     iVec++) {
                  Eigen::Vector3d inPopRot;
                  inPopRot.setZero();
                  for (long i : {0, 1, 2}) {
                    for (long j : {0, 1, 2}) {
                      inPopRot(i) +=
                          rotation(i, j) * inPopulations[iVec](iCalc, j, ibte2);
                    }
                  }
                  for (long i : {0, 1, 2}) {
                    outPopulations[iVec](iCalc, i, ibte1) +=
                        rateOffDiag * inPopRot(i);
                    outPopulations[iVec](iCalc, i, ibte1) +=
                        rate * inPopulations[iVec](iCalc, i, ibte1);
                  }
                }
              } else {
                // case of linewidth construction
                linewidth->operator()(iCalc, 0, ibte1) += rate;
              }
            }
          }
        }
      }
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
    for (long is1 : outerBandStructure.irrStateIterator()) {
      double energy = outerBandStructure.getEnergy(is1);
      auto vel = outerBandStructure.getGroupVelocity(is1);

      auto is1Idx = StateIndex(is1);
      long ind1 = outerBandStructure.stateToBte(is1Idx).get();

      for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
        auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
        double temp = calcStat.temperature;
        double chemPot = calcStat.chemicalPotential;
        // n(n+1)
        double termPop = particle.getPopPopPm1(energy, temp, chemPot);
        double rate = vel.squaredNorm() / boundaryLength * termPop;

        if (switchCase == 0) { // case of matrix construction
          linewidth->operator()(iCalc, 0, ind1) += rate;

        } else if (switchCase == 1) { // case of matrix-vector multiplication
          for (unsigned int iVec = 0; iVec < inPopulations.size(); iVec++) {
            for (long i : {0, 1, 2}) {
              outPopulations[iVec](iCalc, i, ind1) +=
                  rate * inPopulations[iVec](iCalc, i, ind1);
            }
          }

        } else { // case of linewidth construction
          // case of linewidth construction
          linewidth->operator()(iCalc, 0, ind1) += rate;
        }
      }
    }
  }

  // we place the linewidths back in the diagonal of the scattering matrix
  // this because we may need an MPI_allreduce on the linewidths
  if (switchCase == 0) { // case of matrix construction
    long iCalc = 0;
    if (context.getUseSymmetries()) {
      for (long ibte = 0; ibte < numStates; ibte++) {
        auto ibteIdx = BteIndex(ibte);
        for (long i : {0, 1, 2}) {
          auto iCart = CartIndex(i);
          long iMat1 = getSMatrixIndex(ibteIdx, iCart);
          theMatrix(iMat1, iMat1) = linewidth->operator()(iCalc, 0, ibte);
        }
      }
    } else {
      for (long is = 0; is < numStates; is++) {
        theMatrix(is, is) = linewidth->operator()(iCalc, 0, is);
      }
    }
  }
}
