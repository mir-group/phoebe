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

  isMatrixOmega = true;

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
// theMatrix and linewidth is passed: we compute and store in memory the scattering
//       matrix and the diagonal
// inPopulation+outPopulation is passed: we compute the action of the
//       scattering matrix on the in vector, returning outVec = sMatrix*vector
// only linewidth is passed: we compute only the linewidths
void ElScatteringMatrix::builder(VectorBTE *linewidth,
                                 std::vector<VectorBTE> &inPopulations,
                                 std::vector<VectorBTE> &outPopulations) {

  int switchCase = 0;
  if (theMatrix.rows() != 0 && linewidth != nullptr &&
      inPopulations.empty() && outPopulations.empty()) {
    switchCase = 0;
  } else if (theMatrix.rows() == 0 && linewidth == nullptr &&
             !inPopulations.empty() && !outPopulations.empty()) {
    switchCase = 1;
  } else if (theMatrix.rows() == 0 && linewidth != nullptr &&
             inPopulations.empty() && outPopulations.empty()) {
    switchCase = 2;
  } else {
    Error("builder3Ph found a non-supported case");
  }

  if ((linewidth != nullptr) && (linewidth->dimensionality != 1)) {
    Error("The linewidths shouldn't have dimensionality");
  }

  auto particle = outerBandStructure.getParticle();

  int numCalculations = statisticsSweep.getNumCalculations();

  // note: innerNumFullPoints is the number of points in the full grid
  // may be larger than innerNumPoints, when we use ActiveBandStructure
  double norm = 1. / context.getKMesh().prod();

  // precompute Fermi-Dirac populations
  int numOuterIrrStates = outerBandStructure.irrStateIterator().size();
  Eigen::MatrixXd outerFermi(numCalculations, numOuterIrrStates);
  outerFermi.setZero();
#pragma omp parallel for default(none)                                         \
    shared(mpi, outerBandStructure, numCalculations, statisticsSweep, particle,       \
           outerFermi, numOuterIrrStates)
  for (int iBte : mpi->divideWorkIter(numOuterIrrStates)) {
    BteIndex iBteIdx = BteIndex(iBte);
    StateIndex isIdx = outerBandStructure.bteToState(iBteIdx);
    double energy = outerBandStructure.getEnergy(isIdx);
    for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      double temp = calcStat.temperature;
      double chemPot = calcStat.chemicalPotential;
      outerFermi(iCalc, iBte) = particle.getPopulation(energy, temp, chemPot);
    }
  }
  mpi->allReduceSum(&outerFermi);

  int numInnerIrrStates = innerBandStructure.irrStateIterator().size();
  Eigen::MatrixXd innerFermi(numCalculations, numInnerIrrStates);
  innerFermi.setZero();
#pragma omp parallel for default(none)                                         \
    shared(numInnerIrrStates, mpi, innerBandStructure, statisticsSweep,        \
           particle, innerFermi, numCalculations)
  for (int iBte : mpi->divideWorkIter(numInnerIrrStates)) {
    BteIndex iBteIdx = BteIndex(iBte);
    StateIndex isIdx = innerBandStructure.bteToState(iBteIdx);
    double energy = innerBandStructure.getEnergy(isIdx);
    for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      double temp = calcStat.temperature;
      double chemPot = calcStat.chemicalPotential;
      innerFermi(iCalc, iBte) = particle.getPopulation(energy, temp, chemPot);
    }
  }
  mpi->allReduceSum(&innerFermi);

  if (smearing->getType() == DeltaFunction::tetrahedron) {
    Error("Tetrahedron method not supported by electron scattering");
    // that's because it doesn't work with the window the way it's implemented,
    // and we will almost always have a window for electrons
  }

  bool rowMajor = true;
  std::vector<std::tuple<std::vector<int>, int>> kPairIterator =
      getIteratorWavevectorPairs(switchCase, rowMajor);

  HelperElScattering pointHelper(innerBandStructure, outerBandStructure,
                                 statisticsSweep, smearing->getType(), h0);

  bool withSymmetries = context.getUseSymmetries();

  double phononCutoff = 5. / ryToCmm1; // used to discard small phonon energies

  LoopPrint loopPrint("computing scattering matrix", "k-points",
                      kPairIterator.size());

  for (auto t1 : kPairIterator) {
    loopPrint.update();
    auto ik2Indexes = std::get<0>(t1);
    int ik1 = std::get<1>(t1);
    WavevectorIndex ik1Idx(ik1);

    Eigen::Vector3d k1C = outerBandStructure.getWavevector(ik1Idx);
    Eigen::VectorXd state1Energies = outerBandStructure.getEnergies(ik1Idx);
    int nb1 = state1Energies.size();
    Eigen::MatrixXd v1s = outerBandStructure.getGroupVelocities(ik1Idx);
    Eigen::MatrixXcd eigenVector1 = outerBandStructure.getEigenvectors(ik1Idx);

    pointHelper.prepare(k1C, ik2Indexes);
    int nk2 = ik2Indexes.size();
    std::vector<Eigen::Vector3d> allQ3C(nk2);
    std::vector<Eigen::VectorXd> allStates3Energies(nk2);
    std::vector<int> allNb3(nk2);
    std::vector<Eigen::MatrixXcd> allEigenVectors3(nk2);
    std::vector<Eigen::MatrixXd> allV3s(nk2);
    std::vector<Eigen::MatrixXd> allBose3Data(nk2);

    std::vector<Eigen::Vector3d> allK2C(nk2);
    std::vector<Eigen::MatrixXcd> allEigenVectors2(nk2);
    std::vector<Eigen::VectorXd> allState2Energies(nk2);
    std::vector<Eigen::MatrixXd> allV2s(nk2);

    int ik2Counter = -1;
    for (int ik2 : ik2Indexes) {
      ik2Counter++;
      WavevectorIndex ik2Idx(ik2);
      allK2C[ik2Counter] = innerBandStructure.getWavevector(ik2Idx);
      allState2Energies[ik2Counter] = innerBandStructure.getEnergies(ik2Idx);
      allV2s[ik2Counter] = innerBandStructure.getGroupVelocities(ik2Idx);
      allEigenVectors2[ik2Counter] = innerBandStructure.getEigenvectors(ik2Idx);
      auto t2 = pointHelper.get(k1C, ik2);
      allQ3C[ik2Counter] = std::get<0>(t2);
      allStates3Energies[ik2Counter] = std::get<1>(t2);
      allNb3[ik2Counter] = std::get<2>(t2);
      allEigenVectors3[ik2Counter] = std::get<3>(t2);
      allV3s[ik2Counter] = std::get<4>(t2);
      allBose3Data[ik2Counter] = std::get<5>(t2);
    }

    couplingElPhWan->calcCouplingSquared(eigenVector1, allEigenVectors2,
                                         allEigenVectors3, k1C, allK2C, allQ3C);

    ik2Counter = -1;
    for (int ik2 : ik2Indexes) {
      ik2Counter++;

      Eigen::Tensor<double,3> coupling = couplingElPhWan->getCouplingSquared(ik2Counter);

      Eigen::Vector3d k2C = allK2C[ik2Counter];
      auto t3 = innerBandStructure.getRotationToIrreducible(
          k2C, Points::cartesianCoordinates);
      int ik2Irr = std::get<0>(t3);
      Eigen::Matrix3d rotation = std::get<1>(t3);

      WavevectorIndex ik2Idx(ik2);
      WavevectorIndex ik2IrrIdx(ik2Irr);

      Eigen::VectorXd state2Energies = allState2Energies[ik2Counter];
      Eigen::MatrixXd v2s = allV2s[ik2Counter];

      Eigen::MatrixXd bose3Data = allBose3Data[ik2Counter];
      Eigen::MatrixXd v3s = allV3s[ik2Counter];
      Eigen::VectorXd state3Energies = allStates3Energies[ik2Counter];

      int nb2 = state2Energies.size();
      int nb3 = state3Energies.size();

      for (int ib2 = 0; ib2 < nb2; ib2++) {
        double en2 = state2Energies(ib2);
        int is2 = innerBandStructure.getIndex(ik2Idx, BandIndex(ib2));
        int is2Irr = innerBandStructure.getIndex(ik2IrrIdx, BandIndex(ib2));
        StateIndex is2Idx(is2);
        StateIndex is2IrrIdx(is2Irr);
        BteIndex ind2Idx = innerBandStructure.stateToBte(is2IrrIdx);
        int iBte2 = ind2Idx.get();

        for (int ib1 = 0; ib1 < nb1; ib1++) {
          double en1 = state1Energies(ib1);
          int is1 = outerBandStructure.getIndex(ik1Idx, BandIndex(ib1));
          StateIndex is1Idx(is1);
          BteIndex ind1Idx = outerBandStructure.stateToBte(is1Idx);
          int iBte1 = ind1Idx.get();

          for (int ib3 = 0; ib3 < nb3; ib3++) {
            double en3 = state3Energies(ib3);

            // remove small divergent phonon energies
            if (en3 < phononCutoff) {
              continue;
            }

            double delta1, delta2;
            if (smearing->getType() == DeltaFunction::gaussian) {
              delta1 = smearing->getSmearing(en1 - en2 + en3);
              delta2 = smearing->getSmearing(en1 - en2 - en3);
            } else if (smearing->getType() == DeltaFunction::adaptiveGaussian) {
              // Eigen::Vector3d smear = v3s.row(ib3);
              Eigen::Vector3d smear = v1s.row(ib1) - v2s.row(ib1);
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
            for (int iCalc = 0; iCalc < numCalculations; iCalc++) {

              double fermi1 = outerFermi(iCalc, iBte1);
              double fermi2 = innerFermi(iCalc, iBte2);
              double bose3 = bose3Data(iCalc, ib3);

              // Calculate transition probability W+

              double rate = coupling(ib1, ib2, ib3) *
                  ((fermi2 + bose3) * delta1 + (1. - fermi2 + bose3) * delta2) *
                  norm / en3 * pi;

              double rateOffDiagonal = coupling(ib1, ib2, ib3) *
                  (fermi1 * (1. - fermi2) * bose3 * delta1 +
                   fermi2 * (1. - fermi1) * bose3 * delta2) *
                  norm / en3 * pi;

              if (switchCase == 0) {

                if (withSymmetries) {
                  for (int i : {0, 1, 2}) {
                    CartIndex iIndex(i);
                    int iMat1 = getSMatrixIndex(ind1Idx, iIndex);
                    for (int j : {0, 1, 2}) {
                      CartIndex jIndex(j);
                      int iMat2 = getSMatrixIndex(ind2Idx, jIndex);
                      if (theMatrix.indicesAreLocal(iMat1, iMat2)) {
                        if (i == 0 && j == 0) {
                          linewidth->operator()(iCalc, 0, iBte1) += rate;
                        }
                        if (is1 != is2Irr) {
                          theMatrix(iMat1, iMat2) +=
                              rotation.inverse()(i, j) * rateOffDiagonal;
                        }
                      }
                    }
                  }
                } else {
                  if (theMatrix.indicesAreLocal(iBte1, iBte2)) {
                    linewidth->operator()(iCalc, 0, iBte1) += rate;
                  }
                  theMatrix(iBte1, iBte2) += rateOffDiagonal;
                }
              } else if (switchCase == 1) {
                // case of matrix-vector multiplication
                // we build the scattering matrix A = S*n(n+1)

                for (unsigned int iVec = 0; iVec < inPopulations.size();
                     iVec++) {
                  Eigen::Vector3d inPopRot;
                  inPopRot.setZero();
                  for (int i : {0, 1, 2}) {
                    for (int j : {0, 1, 2}) {
                      inPopRot(i) +=
                          rotation.inverse()(i, j)
                          * inPopulations[iVec](iCalc, j, iBte2);
                    }
                  }
                  for (int i : {0, 1, 2}) {
                    if (is1 != is2Irr) {
                      outPopulations[iVec](iCalc, i, iBte1) +=
                          rateOffDiagonal * inPopRot(i);
                    }
                    outPopulations[iVec](iCalc, i, iBte1) +=
                        rate * inPopulations[iVec](iCalc, i, iBte1);
                  }
                }
              } else {
                // case of linewidth construction
                linewidth->operator()(iCalc, 0, iBte1) += rate;
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

  // Average over degenerate eigenstates.
  // we turn it off for now and leave the code if needed in the future
  if (switchCase == 2) {
    degeneracyAveragingLinewidths(linewidth);
  }

  // Add boundary scattering

  if (doBoundary) {
#pragma omp parallel for default(none)                                         \
    shared(outerBandStructure, numCalculations, statisticsSweep,               \
           boundaryLength, particle, outPopulations, inPopulations, linewidth, \
           switchCase)
    for (int is1 : outerBandStructure.irrStateIterator()) {
      StateIndex is1Idx(is1);
      auto vel = outerBandStructure.getGroupVelocity(is1Idx);
      int iBte1 = outerBandStructure.stateToBte(is1Idx).get();
      double rate = vel.squaredNorm() / boundaryLength;

      for (int iCalc = 0; iCalc < numCalculations; iCalc++) {

        if (switchCase == 0) { // case of matrix construction
          linewidth->operator()(iCalc, 0, iBte1) += rate;

        } else if (switchCase == 1) { // case of matrix-vector multiplication
          for (unsigned int iVec = 0; iVec < inPopulations.size(); iVec++) {
            for (int i : {0, 1, 2}) {
              outPopulations[iVec](iCalc, i, iBte1) +=
                  rate * inPopulations[iVec](iCalc, i, iBte1);
            }
          }

        } else { // case of linewidth construction
          // case of linewidth construction
          linewidth->operator()(iCalc, 0, iBte1) += rate;
        }
      }
    }
  }

  // we place the linewidths back in the diagonal of the scattering matrix
  // this because we may need an MPI_allReduce on the linewidths
  if (switchCase == 0) { // case of matrix construction
    int iCalc = 0;
    if (context.getUseSymmetries()) {
      // numStates is defined in scattering.cpp as # of irrStates
      // from the outer band structure
      for (int iBte = 0; iBte < numStates; iBte++) {
        BteIndex iBteIdx(iBte);
        // zero the diagonal of the matrix
        for (int i : {0, 1, 2}) {
          CartIndex iCart(i);
          int iMati = getSMatrixIndex(iBteIdx, iCart);
          for (int j : {0, 1, 2}) {
            CartIndex jCart(j);
            int iMatj = getSMatrixIndex(iBteIdx, jCart);
            theMatrix(iMati, iMatj) = 0.;
          }
          theMatrix(iMati, iMati) += linewidth->operator()(iCalc, 0, iBte);
        }
      }
    } else {
      for (int is = 0; is < numStates; is++) {
        theMatrix(is, is) = linewidth->operator()(iCalc, 0, is);
      }
    }
  }
}
