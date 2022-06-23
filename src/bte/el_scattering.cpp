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
                                       InteractionElPhWan *couplingElPhWan_,
                                       Interaction4El *coupling4El_)
    : ScatteringMatrix(context_, statisticsSweep_, innerBandStructure_,
                       outerBandStructure_),
      couplingElPhWan(couplingElPhWan_), coupling4El(coupling4El_),
      h0(h0_) {

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
      coupling4El(that.coupling4El),
      h0(that.h0), boundaryLength(that.boundaryLength),
      doBoundary(that.doBoundary) {}

ElScatteringMatrix &
ElScatteringMatrix::operator=(const ElScatteringMatrix &that) {
  ScatteringMatrix::operator=(that);
  if (this != &that) {
    couplingElPhWan = that.couplingElPhWan;
    coupling4El = that.coupling4El;
    h0 = that.h0;
    boundaryLength = that.boundaryLength;
    doBoundary = that.doBoundary;
  }
  return *this;
}

// 3 cases:
// theMatrix and linewidth is passed: we compute and store in memory the
// scattering
//       matrix and the diagonal
// inPopulation+outPopulation is passed: we compute the action of the
//       scattering matrix on the in vector, returning outVec = sMatrix*vector
// only linewidth is passed: we compute only the linewidths
void ElScatteringMatrix::builder(VectorBTE *linewidth,
                                 std::vector<VectorBTE> &inPopulations,
                                 std::vector<VectorBTE> &outPopulations) {

  int switchCase = 0;
  if (theMatrix.rows() != 0 && linewidth != nullptr && inPopulations.empty() && outPopulations.empty()) {
    switchCase = 0;
  } else if (theMatrix.rows() == 0 && linewidth == nullptr && !inPopulations.empty() && !outPopulations.empty()) {
    switchCase = 1;
  } else if (theMatrix.rows() == 0 && linewidth != nullptr && inPopulations.empty() && outPopulations.empty()) {
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
  auto numOuterIrrStates = int(outerBandStructure.irrStateIterator().size());
  Eigen::MatrixXd outerFermi(numCalculations, numOuterIrrStates);
  outerFermi.setZero();
  std::vector<size_t> iBtes = mpi->divideWorkIter(numOuterIrrStates);
  int niBtes = iBtes.size();
#pragma omp parallel for default(none)                                \
    shared(mpi, outerBandStructure, numCalculations, statisticsSweep, \
           particle, outerFermi, numOuterIrrStates, niBtes, iBtes)
  for (int iiBte = 0; iiBte < niBtes; iiBte++) {
    int iBte = iBtes[iiBte];
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

  auto numInnerIrrStates = int(innerBandStructure.irrStateIterator().size());
  Eigen::MatrixXd innerFermi(numCalculations, numInnerIrrStates);
  innerFermi.setZero();
  iBtes = mpi->divideWorkIter(numInnerIrrStates);
  niBtes = iBtes.size();
#pragma omp parallel for default(none)                                  \
    shared(numInnerIrrStates, mpi, innerBandStructure, statisticsSweep, \
           particle, innerFermi, numCalculations, niBtes, iBtes)
  for (int iiBte = 0; iiBte < niBtes; iiBte++) {
    int iBte = iBtes[iiBte];
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
                                 statisticsSweep, smearing->getType(), h0, couplingElPhWan);

  bool withSymmetries = context.getUseSymmetries();

  double phononCutoff = 5. / ryToCmm1;// used to discard small phonon energies

  if (couplingElPhWan != nullptr) {
    LoopPrint loopPrint("computing scattering matrix", "k-points",
                        int(kPairIterator.size()));

    for (auto t1 : kPairIterator) {
      loopPrint.update();
      auto ik2Indexes = std::get<0>(t1);
      int ik1 = std::get<1>(t1);
      WavevectorIndex ik1Idx(ik1);

      // dummy call to make pooled coupling calculation work. We need to make sure
      // calcCouplingSquared is called the same # of times. This is also taken
      // care of while generating the indices. Here we call calcCoupling.
      // This block is useful if e.g. we have a pool of size 2, the 1st MPI
      // process has 7 k-points, the 2nd MPI process has 6. This block makes
      // the 2nd process call calcCouplingSquared 7 times as well.
      if (ik1 == -1) {
        Eigen::Vector3d k1C = Eigen::Vector3d::Zero();
        int numWannier = couplingElPhWan->getCouplingDimensions()(0);
        Eigen::MatrixXcd eigenVector1 = Eigen::MatrixXcd::Zero(numWannier, 1);
        couplingElPhWan->cacheElPh(eigenVector1, k1C);
        // since this is just a dummy call used to help other MPI processes
        // compute the coupling, and not to compute matrix elements, we can skip
        // to the next loop iteration
        continue;
      }

      Eigen::Vector3d k1C = outerBandStructure.getWavevector(ik1Idx);
      Eigen::VectorXd state1Energies = outerBandStructure.getEnergies(ik1Idx);
      auto nb1 = int(state1Energies.size());
      Eigen::MatrixXd v1s = outerBandStructure.getGroupVelocities(ik1Idx);
      Eigen::MatrixXcd eigenVector1 = outerBandStructure.getEigenvectors(ik1Idx);

      couplingElPhWan->cacheElPh(eigenVector1, k1C);

      pointHelper.prepare(k1C, ik2Indexes);

      // prepare batches based on memory usage
      auto nk2 = int(ik2Indexes.size());
      int numBatches = couplingElPhWan->estimateNumBatches(nk2, nb1);

      // loop over batches of q1s
      // later we will loop over the q1s inside each batch
      // this is done to optimize the usage and data transfer of a GPU
      for (int iBatch = 0; iBatch < numBatches; iBatch++) {
        // start and end point for current batch
        int start = nk2 * iBatch / numBatches;
        int end = nk2 * (iBatch + 1) / numBatches;
        int batch_size = end - start;

        std::vector<Eigen::Vector3d> allQ3C(batch_size);
        std::vector<Eigen::VectorXd> allStates3Energies(batch_size);
        std::vector<int> allNb3(batch_size);
        std::vector<Eigen::MatrixXcd> allEigenVectors3(batch_size);
        std::vector<Eigen::MatrixXd> allV3s(batch_size);
        std::vector<Eigen::MatrixXd> allBose3Data(batch_size);
        std::vector<Eigen::VectorXcd> allPolarData(batch_size);

        std::vector<Eigen::Vector3d> allK2C(batch_size);
        std::vector<Eigen::MatrixXcd> allEigenVectors2(batch_size);
        std::vector<Eigen::VectorXd> allState2Energies(batch_size);
        std::vector<Eigen::MatrixXd> allV2s(batch_size);

        // do prep work for all values of q1 in current batch,
        // store stuff needed for couplings later
#pragma omp parallel for default(none) shared(allNb3, allEigenVectors3, allV3s, allBose3Data, ik2Indexes, pointHelper, allQ3C, allStates3Energies, batch_size, start, allK2C, allState2Energies, allV2s, allEigenVectors2, k1C, allPolarData)
        for (int ik2Batch = 0; ik2Batch < batch_size; ik2Batch++) {
          int ik2 = ik2Indexes[start + ik2Batch];
          WavevectorIndex ik2Idx(ik2);
          allK2C[ik2Batch] = innerBandStructure.getWavevector(ik2Idx);
          allState2Energies[ik2Batch] = innerBandStructure.getEnergies(ik2Idx);
          allV2s[ik2Batch] = innerBandStructure.getGroupVelocities(ik2Idx);
          allEigenVectors2[ik2Batch] = innerBandStructure.getEigenvectors(ik2Idx);
          auto t2 = pointHelper.get(k1C, ik2);
          allQ3C[ik2Batch] = std::get<0>(t2);
          allStates3Energies[ik2Batch] = std::get<1>(t2);
          allNb3[ik2Batch] = std::get<2>(t2);
          allEigenVectors3[ik2Batch] = std::get<3>(t2);
          allV3s[ik2Batch] = std::get<4>(t2);
          allBose3Data[ik2Batch] = std::get<5>(t2);
          allPolarData[ik2Batch] = std::get<6>(t2);
        }

        couplingElPhWan->calcCouplingSquared(eigenVector1, allEigenVectors2,
                                             allEigenVectors3, allQ3C, allPolarData);

        // do postprocessing loop with batch of couplings
        for (int ik2Batch = 0; ik2Batch < batch_size; ik2Batch++) {
          int ik2 = ik2Indexes[start + ik2Batch];

          Eigen::Tensor<double, 3> coupling =
              couplingElPhWan->getCouplingSquared(ik2Batch);

          Eigen::Vector3d k2C = allK2C[ik2Batch];
          auto t3 = innerBandStructure.getRotationToIrreducible(
              k2C, Points::cartesianCoordinates);
          int ik2Irr = std::get<0>(t3);
          Eigen::Matrix3d rotation = std::get<1>(t3);

          WavevectorIndex ik2Idx(ik2);
          WavevectorIndex ik2IrrIdx(ik2Irr);

          Eigen::VectorXd state2Energies = allState2Energies[ik2Batch];
          Eigen::MatrixXd v2s = allV2s[ik2Batch];

          Eigen::MatrixXd bose3Data = allBose3Data[ik2Batch];
          Eigen::MatrixXd v3s = allV3s[ik2Batch];
          Eigen::VectorXd state3Energies = allStates3Energies[ik2Batch];

          auto nb2 = int(state2Energies.size());
          auto nb3 = int(state3Energies.size());

          Eigen::MatrixXd sinh3Data(nb3, numCalculations);
#pragma omp parallel for collapse(2)
          for (int ib3 = 0; ib3 < nb3; ib3++) {
            for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
              double en3 = state3Energies(ib3);
              double kT = statisticsSweep.getCalcStatistics(iCalc).temperature;
              sinh3Data(ib3, iCalc) = 0.5 / sinh(0.5 * en3 / kT);
            }
          }

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
                  Eigen::Vector3d smear = v1s.row(ib1) - v2s.row(ib2);
                  delta1 = smearing->getSmearing(en1 - en2 + en3, smear);
                  delta2 = smearing->getSmearing(en1 - en2 - en3, smear);
                } else {
                  delta1 = smearing->getSmearing(en3 + en1, is2Idx);
                  delta2 = smearing->getSmearing(en3 - en1, is2Idx);
                }

                if (delta1 <= 0. && delta2 <= 0.) {
                  continue;
                }

                // loop on temperature
                for (int iCalc = 0; iCalc < numCalculations; iCalc++) {

                  double fermi1 = outerFermi(iCalc, iBte1);
                  double fermi2 = innerFermi(iCalc, iBte2);
                  double bose3 = bose3Data(iCalc, ib3);
                  double bose3Symm = sinh3Data(ib3, iCalc); // 1/2/sinh() term

                  // Calculate transition probability W+

                  double rate =
                      coupling(ib1, ib2, ib3)
                      * ((fermi2 + bose3) * delta1
                         + (1. - fermi2 + bose3) * delta2)
                      * norm / en3 * pi;

                  double rateOffDiagonal = -
                                           coupling(ib1, ib2, ib3) * bose3Symm * (delta1 + delta2)
                      * norm / en3 * pi;

                  // double rateOffDiagonal = -
                  // coupling(ib1, ib2, ib3)
                  // * ((1 + bose3 - fermi1) * delta1 + (bose3 + fermi1) * delta2)
                  // * norm / en3 * pi;

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
                          inPopRot(i) += rotation.inverse()(i, j) * inPopulations[iVec](iCalc, j, iBte2);
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
    }
    mpi->barrier();
    loopPrint.close();

  }// end of el-ph lifetimes

  if (coupling4El != nullptr) {

    int numK = innerBandStructure.getNumPoints();

    double norm2 = norm * norm;// numK^2

    std::vector<int> ik3Indexes(numK);
    // populate vector with integers from 0 to numPoints-1
    std::iota(std::begin(ik3Indexes), std::end(ik3Indexes), 0);

    // precompute all fermi population terms
    Eigen::MatrixXd fermiFactor(numStates, numCalculations);
    Eigen::MatrixXd charges(numStates, numCalculations);
#pragma omp parallel for collapse(2)
    for (int is1 = 0; is1 < numStates; ++is1) {
      for (int iCalc = 0; iCalc < numCalculations; ++iCalc) {
        StateIndex is1Idx(is1);
        double en = outerBandStructure.getEnergy(is1Idx);
        auto calcInfo = statisticsSweep.getCalcStatistics(iCalc);
        double temp = calcInfo.temperature;
        double chemPot = calcInfo.chemicalPotential;
        fermiFactor(is1, iCalc) = particle.getPopulation(en, temp, chemPot);

        charges(is1, iCalc) = -1.;
        if (en < chemPot) {// it's a hole!
          charges(is1, iCalc) = 1.;
          fermiFactor(is1, iCalc) = 1. - fermiFactor(is1, iCalc);
        }
      }
    }

    std::vector<double> temperatures(numCalculations);
    for (int iCalc = 0; iCalc < numCalculations; ++iCalc) {
      auto calcInfo = statisticsSweep.getCalcStatistics(iCalc);
      double temp = calcInfo.temperature;
      temperatures[iCalc] = temp;
    }

    if (withSymmetries) {
      Error("Didn't implement 4el interaction with symmetries");
      // because it's a mess to think about!
    }

    LoopPrint loopPrint("computing 4-el scattering matrix", "k-points",
                        int(kPairIterator.size()));

    for (auto t1 : kPairIterator) {// loop over k
      loopPrint.update();
      auto ik2Indexes = std::get<0>(t1);
      int ik1 = std::get<1>(t1);
      WavevectorIndex ik1Idx(ik1);

      Eigen::Vector3d k1C = outerBandStructure.getWavevector(ik1Idx);
      Eigen::VectorXd energies1 = outerBandStructure.getEnergies(ik1Idx);
      auto nb1 = int(energies1.size());
      Eigen::MatrixXcd eigenVectors1 = outerBandStructure.getEigenvectors(ik1Idx);

      // do the Fourier transform of the coupling on the 1st wavevector
      coupling4El->cache1stEl(eigenVectors1, k1C);

      for (int ik2 : ik2Indexes) {// loop over k'
        WavevectorIndex ik2Idx(ik2);
        Eigen::Vector3d k2C = innerBandStructure.getWavevector(ik2Idx);
        Eigen::VectorXd energies2 = innerBandStructure.getEnergies(ik2Idx);
        auto nb2 = int(energies2.size());
        Eigen::MatrixXcd eigenVectors2 = innerBandStructure.getEigenvectors(ik2Idx);

        // do the Fourier transform of the coupling on the 2nd wavevector
        coupling4El->cache2ndEl(eigenVectors2, k2C);

        std::vector<Eigen::Vector3d> k3Cs(numK), k4Cs(numK);
        std::vector<Eigen::MatrixXcd> eigenVectors3(numK), eigenVectors4(numK);
        std::vector<int> ik4Indexes(numK);

#pragma omp parallel for
        for (int ik3 : ik3Indexes) {
          WavevectorIndex ik3Idx(ik3);
          Eigen::Vector3d k3C = innerBandStructure.getWavevector(ik3Idx);

          Eigen::Vector3d k4CTemp = k1C + k2C - k3C;
          Eigen::Vector3d k4CTempCrys = innerBandStructure.getPoints().cartesianToCrystal(k4CTemp);
          // note: I may need to think about what happens if k4C is not on the grid. Ignore?
          int ik4 = innerBandStructure.getPoints().getIndex(k4CTempCrys);
          ik4Indexes[ik3] = ik4;
          WavevectorIndex ik4Idx(ik4);
          Eigen::Vector3d k4C = innerBandStructure.getWavevector(ik4Idx);
          k3Cs[ik3] = k3C;
          k4Cs[ik3] = k4C;
          eigenVectors3[ik3] = innerBandStructure.getEigenvectors(ik3Idx);
          eigenVectors4[ik3] = innerBandStructure.getEigenvectors(ik4Idx);
        }

        coupling4El->calcCouplingSquared(eigenVectors3, eigenVectors4, k3Cs, k4Cs);

        for (int ik3 : ik3Indexes) {
          WavevectorIndex ik3Idx(ik3);
          Eigen::VectorXd energies3 = innerBandStructure.getEnergies(ik3Idx);
          auto nb3 = int(energies3.size());

          int ik4 = ik4Indexes[ik3];
          WavevectorIndex ik4Idx(ik4);
          Eigen::VectorXd energies4 = innerBandStructure.getEnergies(ik4Idx);
          auto nb4 = int(energies4.size());

          auto tC = coupling4El->getCouplingSquared(ik3);
          Eigen::Tensor<double, 4> couplingA = std::get<0>(tC);
          Eigen::Tensor<double, 4> couplingB = std::get<1>(tC);
          Eigen::Tensor<double, 4> couplingC = std::get<2>(tC);

          // add a delta term
          {
#pragma omp parallel for collapse(4)
            for (int ib1 = 0; ib1 < nb1; ++ib1) {
              for (int ib2 = 0; ib2 < nb2; ++ib2) {
                for (int ib3 = 0; ib3 < nb3; ++ib3) {
                  for (int ib4 = 0; ib4 < nb4; ++ib4) {
                    double en1 = energies1(ib1);
                    double en2 = energies2(ib2);
                    double en3 = energies3(ib3);
                    double en4 = energies4(ib4);
                    double deltaA = smearing->getSmearing(en1 + en2 - en3 - en4);
                    double deltaB = smearing->getSmearing(en1 + en3 - en2 - en4);
                    double deltaC = smearing->getSmearing(en3 + en2 - en1 - en4);
                    couplingA(ib1, ib2, ib3, ib4) *= deltaA;
                    couplingB(ib1, ib3, ib2, ib4) *= deltaB;
                    couplingC(ib3, ib2, ib1, ib4) *= deltaC;
                  }
                }
              }
            }
          }

          for (int ib1 = 0; ib1 < nb1; ++ib1) {
            int is1 = innerBandStructure.getIndex(ik1Idx, BandIndex(ib1));
            for (int ib2 = 0; ib2 < nb2; ++ib2) {
              int is2 = innerBandStructure.getIndex(ik2Idx, BandIndex(ib2));
              for (int ib3 = 0; ib3 < nb3; ++ib3) {
                int is3 = innerBandStructure.getIndex(ik3Idx, BandIndex(ib3));
                for (int ib4 = 0; ib4 < nb4; ++ib4) {
                  int is4 = innerBandStructure.getIndex(ik4Idx, BandIndex(ib4));

                  for (int iCalc = 0; iCalc < numCalculations; ++iCalc) {

                    // here we impose the charge conservation rule
                    double chemPot = statisticsSweep.getCalcStatistics(iCalc).chemicalPotential;

                    double fermi1 = fermiFactor(is1, iCalc);
                    double fermi2 = fermiFactor(is2, iCalc);
                    double fermi3 = fermiFactor(is3, iCalc);
                    double fermi4 = fermiFactor(is4, iCalc);

                    double charge1 = charges(is1, iCalc);
                    double charge2 = charges(is2, iCalc);
                    double charge3 = charges(is3, iCalc);
                    double charge4 = charges(is4, iCalc);

                    double rate = 0.;
                    double rateOffDiagonal = 0.;

                    // condition for e-h creation
                    if (charge1 == charge2 + charge3 + charge3) {
                      rate += norm2 * 2
                          * ((1. - fermi2) * (1. - fermi3) * (1. - fermi4) + fermi2 * fermi3 * fermi4)
                          * couplingA(ib1, ib2, ib3, ib4);
                      rateOffDiagonal -=
                          norm2 * 2 * exp(0.5 / temperatures[iCalc] * (energies1(ib1) - energies2(ib2)))
                          * (1. - fermi3) * (1. - fermi4) * (couplingA(ib1, ib2, ib3, ib4) + couplingB(ib1, ib3, ib2, ib4) + couplingC(ib1, ib4, ib3, ib2));
                    } else// in this case is e-e, e-h or h-h scattering
                    {
                      if (charge1 == charge3 && charge2 == charge4) {
                        rate = norm2 * (fermi2 * (1. - fermi3) * (1. - fermi4) + (1. - fermi2) * fermi3 * fermi4)
                            * couplingA(ib1, ib2, ib3, ib4);
                        rateOffDiagonal +=
                            2. * norm2 * exp(0.5 / temperatures[iCalc] * (energies1(ib1) + energies2(ib2) + 2. * chemPot))
                            * (1. - fermi3) * (1. - fermi4) * couplingA(ib1, ib2, ib3, ib4);
                      }
                      if (charge1 == charge2 && charge3 == charge4) {
                        rateOffDiagonal -=
                            2. * norm2 * exp(0.5 / temperatures[iCalc] * (energies1(ib1) - energies2(ib2)))
                            * fermi3 * (1. - fermi4) * couplingB(ib1, ib3, ib2, ib4);
                      }
                      if (charge1 == charge3 && charge2 == charge4) {
                        rateOffDiagonal -=
                            2. * norm2 * exp(0.5 / temperatures[iCalc] * (energies1(ib1) - energies2(ib2)))
                            * fermi4 * (1. - fermi3) * couplingC(ib1, ib4, ib3, ib2);
                      }
                    }

                    if (switchCase == 0) {

                      if (theMatrix.indicesAreLocal(is1, is2)) {
                        linewidth->operator()(iCalc, 0, is1) += rate;
                      }
                      theMatrix(is1, is2) += rateOffDiagonal;

                    } else if (switchCase == 1) {
                      // case of matrix-vector multiplication
                      // we build the scattering matrix A = S*n(n+1)

                      for (unsigned int iVec = 0; iVec < inPopulations.size();
                           iVec++) {
                        for (int i : {0, 1, 2}) {
                          if (is1 != is2) {
                            outPopulations[iVec](iCalc, i, is1) +=
                                rateOffDiagonal * inPopulations[iVec](iCalc, i, is2);
                          }
                          outPopulations[iVec](iCalc, i, is1) +=
                              rate * inPopulations[iVec](iCalc, i, is1);
                        }
                      }
                    } else {
                      // case of linewidth construction
                      linewidth->operator()(iCalc, 0, is1) += rate;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    mpi->barrier();
    loopPrint.close();
  }

  if (switchCase == 1) {
    for (unsigned int iVec = 0; iVec < inPopulations.size(); iVec++) {
      mpi->allReduceSum(&outPopulations[iVec].data);
    }
  } else {
    mpi->allReduceSum(&linewidth->data);
  }

  // Average over degenerate eigenstates.
  // we turn it off for now and leave the code if needed in the future
  if (switchCase == 2) {
    degeneracyAveragingLinewidths(linewidth);
  }

  // Add boundary scattering
  // this must be done after the mpi->allReduce(linewidth)

  if (doBoundary) {
    std::vector<int> is1s = outerBandStructure.irrStateIterator();
    int nis1s = is1s.size();
#pragma omp parallel for default(none) shared(                            \
    outerBandStructure, numCalculations, statisticsSweep, boundaryLength, \
    particle, outPopulations, inPopulations, linewidth, switchCase, nis1s, is1s)
    for (int iis1 = 0; iis1 < nis1s; iis1++) {
      int is1 = is1s[iis1];
      StateIndex is1Idx(is1);
      auto vel = outerBandStructure.getGroupVelocity(is1Idx);
      int iBte1 = outerBandStructure.stateToBte(is1Idx).get();
      double rate = vel.squaredNorm() / boundaryLength;

      for (int iCalc = 0; iCalc < numCalculations; iCalc++) {

        if (switchCase == 0) {// case of matrix construction
          linewidth->operator()(iCalc, 0, iBte1) += rate;

        } else if (switchCase == 1) {// case of matrix-vector multiplication
          for (unsigned int iVec = 0; iVec < inPopulations.size(); iVec++) {
            for (int i = 0; i < 3; i++) {
              outPopulations[iVec](iCalc, i, iBte1) +=
                  rate * inPopulations[iVec](iCalc, i, iBte1);
            }
          }

        } else {// case of linewidth construction
          // case of linewidth construction
          linewidth->operator()(iCalc, 0, iBte1) += rate;
        }
      }
    }
  }

  // we place the linewidths back in the diagonal of the scattering matrix
  // this because we may need an MPI_allReduce on the linewidths
  if (switchCase == 0) {// case of matrix construction
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
