#include "constants.h"
#include "io.h"
#include "mpiHelper.h"
#include "periodic_table.h"
#include "interaction_elph.h"
#include "phel_scattering_matrix.h"

const double phononCutoff = 5. / ryToCmm1; // used to discard small phonon energies

void addPhElScattering(PhElScatteringMatrix &matrix, Context &context, 
                std::vector<std::tuple<int, std::vector<int>>> kqPairIterator, 
                VectorBTE *linewidth) { 

  // TODO wish there was a good way to not always call matrix.getPh or El bandstruct
  // should we use exclude indices in here? 

  Particle elParticle(Particle::electron);
  InteractionElPhWan *couplingElPhWan = matrix.couplingElPhWan;
  ElectronH0Wannier *electronH0 = matrix.electronH0;
  int numAtoms = matrix.getPhBandStructure().getPoints().getCrystal().getNumAtoms();
  int numCalculations = matrix.statisticsSweep.getNumCalculations();
  DeltaFunction *smearing = matrix.smearing;

  Eigen::VectorXd temperatures(numCalculations);
  for (int iCalc=0; iCalc<numCalculations; ++iCalc) {
    auto calcStat = matrix.statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStat.temperature;
    temperatures(iCalc) = temp;
  }

  // note: innerNumFullPoints is the number of points in the full grid
  // may be larger than innerNumPoints, when we use ActiveBandStructure
  // note: in the equations for this rate, because there's an integraton over k,
  // this rate is actually 1/NK (sometimes written N_eFermi).
  double norm = 1. / context.getKMesh().prod();

  // precompute Fermi-Dirac factors
  int numKPoints = matrix.getElBandStructure().getNumPoints();
  int nb1Max = 0;
  for (int ik=0; ik<numKPoints; ++ik) {
    WavevectorIndex ikIdx(ik);
    nb1Max = std::max(nb1Max, int(matrix.getElBandStructure().getEnergies(ikIdx).size()));
  }

  // NOTE statistics sweep is the one for electrons
  // precompute Fermi-Dirac populations
  // TODO can we fit this into the same format as the other ones
  // to call the helper function instead?
  Eigen::Tensor<double,3> fermiTerm(numCalculations, numKPoints, nb1Max);
  fermiTerm.setZero();
  #pragma omp parallel for
  for (int ik : mpi->divideWorkIter(numKPoints)) {
    WavevectorIndex ikIdx(ik);
    Eigen::VectorXd energies = matrix.getElBandStructure().getEnergies(ikIdx);
    int nb1 = energies.size();
    for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
      auto calcStat = matrix.statisticsSweep.getCalcStatistics(iCalc);
      double temp = calcStat.temperature;
      double chemPot = calcStat.chemicalPotential;
      for (int ib=0; ib<nb1; ++ib) {
        fermiTerm(iCalc, ik, ib) = elParticle.getPopPopPm1(energies(ib), temp, chemPot);
      }
    }
  }
  mpi->allReduceSum(&fermiTerm);

  // precompute the q-dependent part of the polar correction ---------
  int numQPoints = matrix.getPhBandStructure().getNumPoints();
  auto qPoints = matrix.getPhBandStructure().getPoints();
  // we just set this to the largest possible number of phonons
  int nb3Max = 3 * matrix.getPhBandStructure().getPoints().getCrystal().getNumAtoms();
  Eigen::MatrixXcd polarData(numQPoints, nb3Max);
  polarData.setZero();
  #pragma omp parallel for
  for (int iq : mpi->divideWorkIter(numQPoints)){
    WavevectorIndex iqIdx(iq);
    auto q3C = matrix.getPhBandStructure().getWavevector(iqIdx);
    auto ev3 = matrix.getPhBandStructure().getEigenvectors(iqIdx);
    Eigen::VectorXcd thisPolar = couplingElPhWan->polarCorrectionPart1(q3C, ev3);
    for (int i=0; i<thisPolar.size(); ++i) {
      polarData(iq, i) = thisPolar(i);
    }
  }
  mpi->allReduceSum(&polarData);

  //---------------

  // k1, k2 are the electronic states, q3 is the phonon
  // this helper returns pairs of the form vector<idxK1, std::vector(allQ3idxs)>
  LoopPrint loopPrint("computing ph-el linewidths","k,q pairs", kqPairIterator.size());


  // iterate over vector<idxK1, std::vector(allK2idxs)>
  // In this loop, k1 is fixed at the top, and we compute
  // it's electronic properties in the outer loop.
  // q3 is the list of iq3Indices, and k2 is determined using k1 and q3
  for (auto t1 : kqPairIterator) {

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
      int numWannier = couplingElPhWan->getCouplingDimensions()(0);
      Eigen::MatrixXcd eigenVector1 = Eigen::MatrixXcd::Zero(numWannier, 1);
      couplingElPhWan->cacheElPh(eigenVector1, k1C);
      // since this is just a dummy call used to help other MPI processes
      // compute the coupling, and not to compute matrix elements, we can skip
      // to the next loop iteration
      continue;
    }

    // store k1 energies, velocities, eigenvectors from elBandstructure
    Eigen::Vector3d k1C = matrix.getElBandStructure().getWavevector(ik1Idx);
    Eigen::VectorXd state1Energies = matrix.getElBandStructure().getEnergies(ik1Idx);
    auto nb1 = int(state1Energies.size());
    Eigen::MatrixXd v1s = matrix.getElBandStructure().getGroupVelocities(ik1Idx);
    Eigen::MatrixXcd eigenVector1 = matrix.getElBandStructure().getEigenvectors(ik1Idx);
    // unlike in el-ph scattering, here we only loop over irr points.
    // This means we need to multiply by the weights of the irr k1s
    // in our integration over the BZ. This returns the list of kpoints
    // that map to this irr kpoint
    double k1Weight = matrix.getElBandStructure().getPoints().
                                getReducibleStarFromIrreducible(ik1).size();

    // precompute first fourier transform + rotation by k1
    couplingElPhWan->cacheElPh(eigenVector1, k1C);

    // prepare batches of q3s based on memory usage (so that this could be done on gpus)?
    auto nq3 = int(iq3Indexes.size());
    int numBatches = couplingElPhWan->estimateNumBatches(nq3, nb1);

    // loop over batches of q3s
    // later we will loop over the q3s inside each batch
    // this is done to optimize the usage and data transfer of a GPU
    for (int iBatch = 0; iBatch < numBatches; iBatch++) {

      // start and end point for current batch of q3s
      int start = nq3 * iBatch / numBatches;
      int end = nq3 * (iBatch + 1) / numBatches;
      int batch_size = end - start;

      std::vector<Eigen::Vector3d> allQ3C(batch_size);
      std::vector<Eigen::MatrixXcd> allEigenVectors3(batch_size);
      std::vector<Eigen::VectorXcd> allPolarData(batch_size);

      // do prep work for all values of q3 in current batch,
      // store stuff needed for couplings later
      //
      // loop over each iq3 in the batch of q3s
      #pragma omp parallel for
      for (int iq3Batch = 0; iq3Batch < batch_size; iq3Batch++) {

        int iq3 = iq3Indexes[start + iq3Batch];
        WavevectorIndex iq3Idx(iq3);

        allPolarData[iq3Batch] = polarData.row(iq3);
        allEigenVectors3[iq3Batch] = matrix.getPhBandStructure().getEigenvectors(iq3Idx);
        allQ3C[iq3Batch] = matrix.getPhBandStructure().getWavevector(iq3Idx);
      }

      // precompute the k2 indices such that k2-k1=q3, where k1 is fixed
      std::vector<Eigen::Vector3d> allK2C(batch_size);
      #pragma omp parallel for
      for (int iq3Batch = 0; iq3Batch < batch_size; iq3Batch++) {

        int iq3 = iq3Indexes[start + iq3Batch];
        auto iq3Index = WavevectorIndex(iq3);
        Eigen::Vector3d q3C = matrix.getPhBandStructure().getWavevector(iq3Index);

        // k' = k + q : phonon absorption
        Eigen::Vector3d k2C = q3C + k1C;
        allK2C[iq3Batch] = k2C;
      }

      // calculate the state energies, vs, eigs of all k2 points
      bool withVelocities = false;
      if (smearing->getType() == DeltaFunction::adaptiveGaussian) {
        withVelocities = true;
      }
      bool withEigenvectors = true; // we need these below to calculate coupling
      auto tHelp = electronH0->populate(allK2C, withVelocities, withEigenvectors);

      std::vector<Eigen::VectorXd> allStates2Energies = std::get<0>(tHelp);
      std::vector<Eigen::MatrixXcd> allEigenVectors2 = std::get<1>(tHelp);
      std::vector<Eigen::Tensor<std::complex<double>,3>> allStates2Velocities = std::get<2>(tHelp);

      // Generate couplings for fixed k1, all k2s and all Q3Cs
      couplingElPhWan->calcCouplingSquared(eigenVector1, allEigenVectors2,
                                          allEigenVectors3, allQ3C, allPolarData);

      // do postprocessing loop with batch of couplings to calculate the scattering rates
      for (int iq3Batch = 0; iq3Batch < batch_size; iq3Batch++) {

        int iq3 = iq3Indexes[start + iq3Batch];
        WavevectorIndex iq3Idx(iq3);

        Eigen::VectorXd state2Energies = allStates2Energies[iq3Batch];
        auto nb2 = int(state2Energies.size());
        // for gpu would replace with compute OTF
        Eigen::VectorXd state3Energies = matrix.getPhBandStructure().getEnergies(iq3Idx); // iq3Idx

        // NOTE: these loops are already set up to be applicable to gpus
        // the precomputaton of the smearing values and the open mp loops could
        // be converted to GPU relevant version
        int nb3 = state3Energies.size();
        Eigen::Tensor<double,3> smearingValues(nb1, nb2, nb3);

        #pragma omp parallel for collapse(3)
        for (int ib2 = 0; ib2 < nb2; ib2++) {
          for (int ib1 = 0; ib1 < nb1; ib1++) {
            for (int ib3 = 0; ib3 < nb3; ib3++) {

              double en2 = state2Energies(ib2);
              double en1 = state1Energies(ib1);
              double en3 = state3Energies(ib3);

              // remove small divergent phonon energies
              if (en3 < phononCutoff) {
                smearingValues(ib1, ib2, ib3) = 0.;
                continue;
              }
              double delta;

              // NOTE: for gpus, if statements need to be moved outside, as
              // this creates load balancing issues for gpu threads and cpu must
              // dispatch this decision info
              if (smearing->getType() == DeltaFunction::gaussian) {
                delta = smearing->getSmearing(en1 - en2 + en3);
              } else {
                Eigen::Vector3d smear = v1s.row(ib1);
                for (int i : {0,1,2}) {
                  smear(i) -= allStates2Velocities[iq3Batch](ib2, ib2, i).real();
                }
                delta = smearing->getSmearing(en1 - en2 + en3, smear);
              }
              smearingValues(ib1, ib2, ib3) = std::max(delta, 0.);
            }
          }
        }

        Eigen::Tensor<double, 3> coupling = couplingElPhWan->getCouplingSquared(iq3Batch);
        // symmetrize the elph coupling tensor
        matrix.symmetrizeCoupling(coupling, state1Energies, state2Energies, state3Energies);

        #pragma omp parallel for collapse(2)
        for (int ib3 = 0; ib3 < nb3; ib3++) {
          for (int iCalc = 0; iCalc < numCalculations; iCalc++) {

            int is3 = matrix.getPhBandStructure().getIndex(iq3Idx, BandIndex(ib3));
            // the BTE index is an irr point which indexes VectorBTE objects
            // like the linewidths + scattering matrix -- as these are
            // only allocated for irr points when sym is on
            StateIndex isIdx3(is3);
            BteIndex ibteIdx3= matrix.getPhBandStructure().stateToBte(isIdx3);
            int ibte3 = ibteIdx3.get();

            for (int ib1 = 0; ib1 < nb1; ib1++) {
              for (int ib2 = 0; ib2 < nb2; ib2++) {
              // loop on temperature
                // https://arxiv.org/pdf/1409.1268.pdf
                // double rate =
                //    coupling(ib1, ib2, ib3)
                //    * (fermi(iCalc, ik1, ib1) - fermi(iCalc, ik2, ib2))
                //    * smearing_values(ib1, ib2, ib3) * norm / en3 * pi;

                // NOTE: although the expression above is formally correct,
                // fk-fk2 could be negative due to numerical noise.
                // so instead, we do:
                // fk-fk2 ~= dfk/dek dwq
                // However, we don't include the dwq here, as this is gSE^2, which
                // includes a factor of (1/wq)
                double rate =
                    coupling(ib1, ib2, ib3) * fermiTerm(iCalc, ik1, ib1)
                    * smearingValues(ib1, ib2, ib3)
                    * norm / temperatures(iCalc) * pi * k1Weight;

                // case of linewidth construction (the only case, for ph-el)
                linewidth->operator()(iCalc, 0, ibte3) += rate;

                //NOTE: for eliashberg function, we could here add another vectorBTE object
                // as done with the linewidths here, slightly modified

              }
            }
          }
        }
      }
    }
  }
  mpi->barrier(); 
  // I prefer to close loopPrint after the MPI barrier: all MPI are synced here
  loopPrint.close();

}
