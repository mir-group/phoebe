#include "phel_scattering.h"
#include "constants.h"
#include "io.h"
#include "mpiHelper.h"
#include "periodic_table.h"

std::vector<std::tuple<int, std::vector<int>>> PhElScatteringMatrix::getIteratorWavevectorPairs() {
  std::vector<std::tuple<int, std::vector<int>>> pairIterator;

  // here I parallelize over ik1
  // which is the outer loop on q-points
  int numKPoints = getElBandStructure().getNumPoints();
  std::vector<size_t> k1Iterator = mpi->divideWorkIter(numKPoints);

  // I don't parallelize the inner band structure, the inner loop
  std::vector<int> k2Iterator(numKPoints);
  // populate vector with integers from 0 to numPoints-1
  std::iota(std::begin(k2Iterator), std::end(k2Iterator), 0);

  for (size_t ik1 : k1Iterator) {
    auto t = std::make_tuple(int(ik1), k2Iterator);
    pairIterator.push_back(t);
  }
  return pairIterator;
}

PhElScatteringMatrix::PhElScatteringMatrix(Context &context_,
                                           StatisticsSweep &statisticsSweep_,
                                           BaseBandStructure &elBandStructure_,
                                           BaseBandStructure &phBandStructure_,
                                           InteractionElPhWan &couplingElPhWan_,
                                           ElectronH0Wannier &electronH0_)
    : ScatteringMatrix(context_, statisticsSweep_, elBandStructure_,
                       phBandStructure_),
      couplingElPhWan(couplingElPhWan_), electronH0(electronH0_) {

  isMatrixOmega = true;
  highMemory = false;
}

PhElScatteringMatrix::PhElScatteringMatrix(const PhElScatteringMatrix &that)
    : ScatteringMatrix(that), couplingElPhWan(that.couplingElPhWan),
      electronH0(that.electronH0) {}

PhElScatteringMatrix &
PhElScatteringMatrix::operator=(const PhElScatteringMatrix &that) {
  ScatteringMatrix::operator=(that);
  if (this != &that) {
    couplingElPhWan = that.couplingElPhWan;
    electronH0 = that.electronH0;
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
void PhElScatteringMatrix::builder(VectorBTE *linewidth,
                                   std::vector<VectorBTE> &inPopulations,
                                   std::vector<VectorBTE> &outPopulations) {
  (void) inPopulations;
  (void) outPopulations;

  if (linewidth == nullptr) {
    Error("builder3Ph found a non-supported case");
  }
  if (linewidth->dimensionality != 1) {
    Error("Linewidths shouldn't have dimensionality");
  }

  Particle elParticle(Particle::electron);

  int numCalculations = statisticsSweep.getNumCalculations();

  // note: innerNumFullPoints is the number of points in the full grid
  // may be larger than innerNumPoints, when we use ActiveBandStructure
  double norm = 1. / context.getQMesh().prod();

  // precompute Fermi-Dirac factors
  int numKPoints = getElBandStructure().getNumPoints();
  int nb1Max = 0;
  for (int ik=0; ik<numKPoints; ++ik) {
    WavevectorIndex ikIdx(ik);
    nb1Max = std::max(nb1Max, int(getElBandStructure().getEnergies(ikIdx).size()));
  }
  // precompute Fermi-Dirac populations
  Eigen::Tensor<double,3> fermiTerm(numCalculations, numKPoints, nb1Max);
  fermiTerm.setZero();
#pragma omp parallel for
  for (int ik : mpi->divideWorkIter(numKPoints)) {
    WavevectorIndex ikIdx(ik);
    Eigen::VectorXd energies = getElBandStructure().getEnergies(ikIdx);
    int nb1 = energies.size();
    for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      double temp = calcStat.temperature;
      double chemPot = calcStat.chemicalPotential;
      for (int ib=0; ib<nb1; ++ib) {
        fermiTerm(iCalc, ik, ib) = elParticle.getPopPopPm1(energies(ib), temp, chemPot);
      }
    }
  }
  mpi->allReduceSum(&fermiTerm);

  Eigen::VectorXd temperatures(numCalculations);
  for (int iCalc=0; iCalc<numCalculations; ++iCalc) {
    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStat.temperature;
    temperatures(iCalc) = temp;
  }

  // precompute the q-dependent part of the polar correction ---------

  int numQPoints = getPhBandStructure().getNumPoints();
  auto qPoints = getPhBandStructure().getPoints();
  int nb3 = getPhBandStructure().getNumBands();
  Eigen::MatrixXcd polarData(numQPoints, nb3);
  polarData.setZero();
#pragma omp parallel for
  for (int iq : mpi->divideWorkIter(numQPoints)){
    WavevectorIndex iqIdx(iq);
    auto q3C = getPhBandStructure().getWavevector(iqIdx);
    auto ev3 = getPhBandStructure().getEigenvectors(iqIdx);
    polarData.row(iq) = couplingElPhWan.polarCorrectionPart1(q3C, ev3);
  }
  mpi->allReduceSum(&polarData);

  //---------------

  if (smearing->getType() == DeltaFunction::tetrahedron) {
    Error("Tetrahedron method not supported by electron scattering");
    // that's because it doesn't work with the window the way it's implemented,
    // and we will almost always have a window for electrons
  }

  auto kPairIterator = getIteratorWavevectorPairs();

  double phononCutoff = 5. / ryToCmm1; // used to discard small phonon energies

  LoopPrint loopPrint("computing scattering matrix", "k-points", numKPoints);

  for (auto t1 : kPairIterator) {
    loopPrint.update();

    int ik1 = std::get<0>(t1);
    WavevectorIndex ik1Idx(ik1);
    auto iq3Indexes = std::get<1>(t1);

    // dummy call to make pooled coupling calculation work. We need to make sure
    // calcCouplingSquared is called the same # of times. This is also taken
    // care of while generating the indices. Here we call calcCoupling.
    // This block is useful if e.g. we have a pool of size 2, the 1st MPI
    // process has 7 k-points, the 2nd MPI process has 6. This block makes
    // the 2nd process call calcCouplingSquared 7 times as well.
    if (ik1 == -1) {
      Eigen::Vector3d k1C = Eigen::Vector3d::Zero();
      int numWannier = couplingElPhWan.getCouplingDimensions()(0);
      Eigen::MatrixXcd eigenVector1 = Eigen::MatrixXcd::Zero(numWannier, 1);
      couplingElPhWan.cacheElPh(eigenVector1, k1C);
      // since this is just a dummy call used to help other MPI processes
      // compute the coupling, and not to compute matrix elements, we can skip
      // to the next loop iteration
      continue;
    }

    Eigen::Vector3d k1C = getElBandStructure().getWavevector(ik1Idx);
    Eigen::VectorXd state1Energies = getElBandStructure().getEnergies(ik1Idx);
    auto nb1 = int(state1Energies.size());
    Eigen::MatrixXd v1s = getElBandStructure().getGroupVelocities(ik1Idx);
    Eigen::MatrixXcd eigenVector1 = getElBandStructure().getEigenvectors(ik1Idx);

    couplingElPhWan.cacheElPh(eigenVector1, k1C);

    // prepare batches based on memory usage
    auto nq3 = int(iq3Indexes.size());
    int numBatches = couplingElPhWan.estimateNumBatches(nq3, nb1);

    // loop over batches of q1s
    // later we will loop over the q1s inside each batch
    // this is done to optimize the usage and data transfer of a GPU
    for (int iBatch = 0; iBatch < numBatches; iBatch++) {
      // start and end point for current batch
      int start = nq3 * iBatch / numBatches;
      int end = nq3 * (iBatch + 1) / numBatches;
      int batch_size = end - start;

      std::vector<Eigen::Vector3d> allQ3C(batch_size);
      std::vector<Eigen::MatrixXcd> allEigenVectors3(batch_size);
      std::vector<Eigen::VectorXcd> allPolarData(batch_size);

      // do prep work for all values of q1 in current batch,
      // store stuff needed for couplings later
#pragma omp parallel for
      for (int iq3Batch = 0; iq3Batch < batch_size; iq3Batch++) {
        int iq3 = iq3Indexes[start + iq3Batch];
        WavevectorIndex iq3Idx(iq3);

        allPolarData[iq3Batch] = polarData.row(iq3);
        allEigenVectors3[iq3Batch] = getPhBandStructure().getEigenvectors(iq3Idx);
        allQ3C[iq3Batch] = getPhBandStructure().getWavevector(iq3Idx);
      }

      // precompute the q indices such that k'-k=q at fixed k
      std::vector<Eigen::Vector3d> allK2C(batch_size);
#pragma omp parallel for
      for (int iq3Batch = 0; iq3Batch < batch_size; iq3Batch++) {
        int iq3 = iq3Indexes[start + iq3Batch];
        auto iq3Index = WavevectorIndex(iq3);
        Eigen::Vector3d q3C = getPhBandStructure().getWavevector(iq3Index);
        // k' = k + q : phonon absorption
        Eigen::Vector3d k2C = q3C + k1C;
        allK2C[iq3Batch] = k2C;
      }

      auto tHelp = electronH0.populate(allK2C);
      std::vector<Eigen::VectorXd> allStates2Energies = std::get<0>(tHelp);
      std::vector<Eigen::MatrixXcd> allEigenVectors2 = std::get<1>(tHelp);
      std::vector<Eigen::Tensor<std::complex<double>,3>> allStates2Velocities = std::get<2>(tHelp);

      couplingElPhWan.calcCouplingSquared(eigenVector1, allEigenVectors2,
                                          allEigenVectors3, allQ3C, allPolarData);

      // do postprocessing loop with batch of couplings
      for (int iq3Batch = 0; iq3Batch < batch_size; iq3Batch++) {
        int iq3 = iq3Indexes[start + iq3Batch];
        WavevectorIndex iq3Idx(iq3);

        Eigen::Tensor<double, 3> coupling =
            couplingElPhWan.getCouplingSquared(iq3Batch);

        Eigen::VectorXd state2Energies = allStates2Energies[iq3Batch];
        Eigen::Tensor<std::complex<double>,3> v2s = allStates2Velocities[iq3Batch];

        Eigen::VectorXd state3Energies = getPhBandStructure().getEnergies(iq3Idx);

        auto nb2 = int(state2Energies.size());

        Eigen::Tensor<double,3> smearing_values(nb1, nb2, nb3);

#pragma omp parallel for collapse(3)
        for (int ib2 = 0; ib2 < nb2; ib2++) {
          for (int ib1 = 0; ib1 < nb1; ib1++) {
            for (int ib3 = 0; ib3 < nb3; ib3++) {
              double en2 = state2Energies(ib2);
              double en1 = state1Energies(ib1);
              double en3 = state3Energies(ib3);
              // remove small divergent phonon energies
              if (en3 < phononCutoff) {
                smearing_values(ib1, ib2, ib3) = 0.;
                continue;
              }
              double delta;
              if (smearing->getType() == DeltaFunction::gaussian) {
                delta = smearing->getSmearing(en1 - en2 + en3);
              } else {
                Eigen::Vector3d smear = v1s.row(ib1);
                for (int i : {0,1,2}) {
                  smear(i) -= v2s(ib2, ib2, i).real();
                }
                delta = smearing->getSmearing(en1 - en2 + en3, smear);
              }
              smearing_values(ib1, ib2, ib3) = std::max(delta, 0.);
            }
          }
        }

#pragma omp parallel for collapse(2)
        for (int ib3 = 0; ib3 < nb3; ib3++) {
          for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
            int is3 = getPhBandStructure().getIndex(iq3Idx, BandIndex(ib3));
            for (int ib1 = 0; ib1 < nb1; ib1++) {
              for (int ib2 = 0; ib2 < nb2; ib2++) {
              // loop on temperature
                // https://arxiv.org/pdf/1409.1268.pdf
                // double rate =
                //    coupling(ib1, ib2, ib3)
                //    * (fermi(iCalc, ik1, ib1) - fermi(iCalc, ik2, ib2))
                //    * smearing_values(ib1, ib2, ib3) * norm / en3 * pi;

                // note: although the expression above is formally correct,
                // fk-fk2 could be negative due to numerical noise.
                // so instead, we do:
                // fk-fk2 ~= dfk/dek dwq

                double rate =
                    coupling(ib1, ib2, ib3) * fermiTerm(iCalc, ik1, ib1)
                    * smearing_values(ib1, ib2, ib3)
                    * norm / temperatures(iCalc) * pi;

                // case of linewidth construction
                linewidth->operator()(iCalc, 0, is3) += rate;
              }
            }
          }
        }
      }
    }
  }

  mpi->allReduceSum(&linewidth->data);
  // I prefer to close loopPrint after the MPI barrier: all MPI are synced here
  loopPrint.close();

  // Average over degenerate eigenstates.
  // we turn it off for now and leave the code if needed in the future
  degeneracyAveragingLinewidths(linewidth);
}
