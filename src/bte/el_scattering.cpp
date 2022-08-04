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


// alternative copy constructor, which replaces band structure
ElScatteringMatrix::ElScatteringMatrix(ElScatteringMatrix& elScatteringMatrix_,
                        ActiveBandStructure &outerBandStructure_,
                        ActiveBandStructure &innerBandStructure_)
                      : ScatteringMatrix(elScatteringMatrix_, outerBandStructure_, innerBandStructure_),
                        couplingElPhWan(elScatteringMatrix_.couplingElPhWan), h0(elScatteringMatrix_.h0) {

    boundaryLength = elScatteringMatrix_.boundaryLength;
    doBoundary = elScatteringMatrix_.doBoundary;
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
    Error("builderElPh found a non-supported case");
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

#pragma omp parallel for
      for (int ik2Batch = 0; ik2Batch < batch_size; ik2Batch++) {
        symmetrizeCoupling(
            couplingElPhWan->getCouplingSquared(ik2Batch),
            state1Energies, allState2Energies[ik2Batch], allStates3Energies[ik2Batch]
        );
      }

      // do postprocessing loop with batch of couplings
      for (int ik2Batch = 0; ik2Batch < batch_size; ik2Batch++) {
        int ik2 = ik2Indexes[start + ik2Batch];

        Eigen::Tensor<double, 3>& coupling =
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

                //double fermi1 = outerFermi(iCalc, iBte1);
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

/*
 *   Eigen::Vector3d crystalToCartesian(const Eigen::Vector3d &point);
 *
 *   see also functions about folding to ws cell
 *
 * */

void ElScatteringMatrix::addMagneticTerm(const Eigen::Vector3d& magneticField) {

  //bool debug = false;

  // determine the delta k spacing in cartesian coordinates
  auto points = outerBandStructure.getPoints();

  // TODO 0) should we transform the kMesh?
  Eigen::Vector3i kMesh = std::get<0>(points.getMesh());
  Eigen::Vector3d deltaK;
  for (int i : {0,1,2}) {
    deltaK(i) = 1. / kMesh(i);
  }

  auto crystal = points.getCrystal();
  Eigen::Matrix3d R = points.getCrystal().getReciprocalUnitCell();
  // deltaK = R * deltaK;
  //for (int i : {0,1,2} ) {
    // Eigen::Vector3d bi = crystal.getReciprocalUnitCell().col(i);
    // double dilation = bi.norm();
    // deltaK[i] *= dilation;
    //  }


  Eigen::Matrix3d spaceFactor;
  {
    Eigen::Matrix3d mat1, mat2;
    mat1.setZero();
    mat2.setZero();
    for (int i : {0,1,2}) {
      Eigen::Vector3d bi = crystal.getReciprocalUnitCell().col(i);
      double biNorm = bi.norm();
      mat2(i,i) = 1./biNorm;

      bi.array() /= bi.norm();
      mat1.col(i) = bi;
    }
    spaceFactor = mat1.inverse() * mat2;
  }

  // convert the kspacing from crystal to cartesian


  // first, compute all values of the vector Omega=vxB
  // (vxB) is in Cartesian coordinates
  // (the input B vector is in cartesian coordinates too)
  Eigen::MatrixXd correction(numStates,3);
  correction.setZero();
  for (int is=0; is<numStates; ++is) {
    StateIndex isIdx(is);
    // TODO 1) the velocity is cartesian, do we need to transform it?
    Eigen::Vector3d vel = outerBandStructure.getGroupVelocity(isIdx);
    // TODO 2) does magnetic field need to be transformed?

    Eigen::Vector3d y  = vel.cross(magneticField);
    Eigen::Vector3d z = y.transpose() * spaceFactor;

    for (int i : {0,1,2}) {
      correction(is, i) = -z(i) * kMesh(i); // /(deltaK[i]);
    }
  }

  // loop over every state because we do not know
  // if k+1 or k-1 will be on another process. Therefore, we cannot restrict the loop
  // to only local state indices
  int numKpoints = outerBandStructure.getNumPoints();
  for (int is=0; is<numStates; ++is) {

    // generate basic state info
    StateIndex isIdx(is);
    auto t = outerBandStructure.getIndex(is);
    WavevectorIndex ikIdx = std::get<0>(t);
    BandIndex ibIdx = std::get<1>(t);

    // kpoint coordinate values in range [0,1] (crystal coords)
    Eigen::Vector3d k = points.getPointCoordinates(ikIdx.get(), Points::crystalCoordinates);
    Eigen::Vector3d k3idx;

    // kpoint coordinate values in range [0, kMeshX] --> this is i,j,k direction index
    for (int i : {0,1,2}) {
      k3idx(i) = k(i) * kMesh(i);
    }

    // for each direction of the finite differences
    for (int kDisplacement : {0, 1, 2}) {

      // newK has coordinates in range [0, kMeshX]
      Eigen::Vector3d kPlus = k3idx;
      Eigen::Vector3d kMins = k3idx;
      Eigen::Vector3d kPlus2 = k3idx;
      Eigen::Vector3d kMins2 = k3idx;

      // add integer displacement, so that we have the forwards/backwards point indices
      kPlus(kDisplacement) += 1;
      kMins(kDisplacement) -= 1;
      kPlus2(kDisplacement) += 2;
      kMins2(kDisplacement) -= 2;

      // TODO does this work with noncubic things?
      // I think yes for grid of crystal coords.
      // TODO should we also use modulo to wrap around? does band structure wrap for us?
      // kPlus/kMinus has coordinates in range [0, 1]
      for (int i : {0,1,2}) {
        kPlus(i) /= float(kMesh(i));
        kMins(i) /= float(kMesh(i));
        kPlus2(i) /= float(kMesh(i));
        kMins2(i) /= float(kMesh(i));
      }

      // find index of kPlus/kMinus in points
      int ikPlusRed = outerBandStructure.getPointIndex(kPlus, true);
      int ikMinsRed = outerBandStructure.getPointIndex(kMins, true);
      int ikPlusRed2 = outerBandStructure.getPointIndex(kPlus2, true);
      int ikMinsRed2 = outerBandStructure.getPointIndex(kMins2, true);

      // if symmetries are turned on, we need to map the global kmesh
      // index to a irreducible point, as this is what will be used
      // by the interalDiagonal element
      int ikPlus = points.asIrreducibleIndex(ikPlusRed);
      int ikMins = points.asIrreducibleIndex(ikMinsRed);
      int ikPlus2 = points.asIrreducibleIndex(ikPlusRed2);
      int ikMins2 = points.asIrreducibleIndex(ikMinsRed2);

      // TODO -- I think this should be a series of if statements.
      // if we have neither +/-, we're in the middle of a region that
      // was removed. We should expect (hope?) dk f ~ 0
      // TODO the other scenario will be if one or the other kpoint is missing,
      // I think in that case -- could we use forward/backward difference

      // TODO why is it possible for ikPlus to be -1? -- maybe this is if +1/-1 took us off the mesh
      if (ikPlus2 > numKpoints || ikMins2 > numKpoints) { continue; }
      if (ikPlus == -1 || ikMins == -1 || ikMins2 == -1 || ikPlus2 == -1) { continue; }

      WavevectorIndex ikPlusIdx(ikPlus);
      WavevectorIndex ikMinsIdx(ikMins);
      WavevectorIndex ikPlusIdx2(ikPlus2);
      WavevectorIndex ikMinsIdx2(ikMins2);

      // find the kpoint in the kpoint list to get its state index
      // TODO different bands at two kpoints
      int isPlus = outerBandStructure.getIndex(ikPlusIdx, ibIdx);
      int isMins = outerBandStructure.getIndex(ikMinsIdx, ibIdx);
      int isPlus2 = outerBandStructure.getIndex(ikPlusIdx2, ibIdx);
      int isMins2 = outerBandStructure.getIndex(ikMinsIdx2, ibIdx);

      StateIndex istatePlus(isPlus);
      StateIndex istateMins(isMins);
      int ibtePlus = outerBandStructure.stateToBte(istatePlus).get();
      int ibteMins = outerBandStructure.stateToBte(istateMins).get();
      StateIndex istatePlus2(isPlus2);
      StateIndex istateMins2(isMins2);
      int ibtePlus2 = outerBandStructure.stateToBte(istatePlus2).get();
      int ibteMins2 = outerBandStructure.stateToBte(istateMins2).get();

      // before proceeding, we have to use the same coordinates,
      // cartesian or crystal, for both the derivative and omega.
      // TODO update these corrections as below
      if(highMemory) {
        if(theMatrix.indicesAreLocal(ibtePlus,ibtePlus)) {
          theMatrix(ibtePlus, ibtePlus) += correction(is,kDisplacement);
          internalDiagonal(0, 0, ibtePlus) += correction(is,kDisplacement);
        }
        if(theMatrix.indicesAreLocal(ibteMins, ibteMins)) {
          theMatrix(ibteMins, ibteMins) -= correction(is,kDisplacement);
          internalDiagonal(0, 0, ibteMins) -= correction(is,kDisplacement);
        }
        if(theMatrix.indicesAreLocal(ibtePlus2,ibtePlus2)) {
          theMatrix(ibtePlus2, ibtePlus2) += correction(is,kDisplacement);
          internalDiagonal(0, 0, ibtePlus2) += correction(is,kDisplacement);
        }
        if(theMatrix.indicesAreLocal(ibteMins2, ibteMins2)) {
          theMatrix(ibteMins2, ibteMins2) -= correction(is,kDisplacement);
          internalDiagonal(0, 0, ibteMins2) -= correction(is,kDisplacement);
        }
      } else {
        for(int iCalc = 0; iCalc < numCalculations; iCalc++) {
          if(ibtePlus < numStates) {
            internalDiagonal(iCalc, 0, ibtePlus) += (8./12.) * correction(is, kDisplacement);
          }
          if(ibteMins < numStates) {
            internalDiagonal(iCalc, 0, ibteMins) -= (8./12.) * correction(is, kDisplacement);
          }
          if(ibtePlus2 < numStates) {
            internalDiagonal(iCalc, 0, ibtePlus2) -= (1./12.) * correction(is, kDisplacement);
          }
          if(ibteMins2 < numStates) {
            internalDiagonal(iCalc, 0, ibteMins2) += (1./12.) * correction(is, kDisplacement);
          }
        }
      }
    }
  }
}

// =============================================================================================

std::vector<VectorBTE> ElScatteringMatrix::addMagneticTerm(const Eigen::Vector3d& magneticField,
                                                                  std::vector<VectorBTE> n) {

  // get information about how the points are spaced
  auto points = outerBandStructure.getPoints();
  Eigen::Vector3i kMesh = std::get<0>(points.getMesh());
  auto crystal = points.getCrystal();
  Eigen::Matrix3d R = points.getCrystal().getReciprocalUnitCell();

  // objects that will be corrected and returned
  VectorBTE nEB = n[0];
  VectorBTE nTB = n[1];

  // loop over every state because we do not know
  // if k+1 or k-1 will be on another process. Therefore, we cannot restrict the loop
  // to only local state indices
  for (int is=0; is<numStates; ++is) {

    // generate basic state info
    StateIndex isIdx(is);
    auto t = outerBandStructure.getIndex(is);
    WavevectorIndex ikIdx = std::get<0>(t);
    BandIndex ibIdx = std::get<1>(t);

    // kpoint coordinates in crystal basis
    Eigen::Vector3d k = points.getPointCoordinates(ikIdx.get(), Points::crystalCoordinates);
    Eigen::Vector3d k3idx;
    // convert crystal coordintes to an index i,j,k on the mesh
    for (int i : {0, 1, 2}) {
      k3idx(i) = k(i) * kMesh(i);
    }

    // establish the k difference in cartesian coords
    if(kMesh[0] < 4 || kMesh[1] < 4 || kMesh[2] < 4) Error("Cannot perform finite differences with kMesh < 2.");
    Eigen::MatrixXd kdiff(3,3);
    for(int i : {0,1,2}) {
      Eigen::Vector3d k1;
      Eigen::Vector3d k2;
      // define x, y, z increment
      k1(i) = 1./kMesh(i);
      k2(i) = 2./kMesh(i);
      // transform to cartesian
      Eigen::Vector3d k1Cart = R*k1;
      Eigen::Vector3d k2Cart = R*k2;
      auto kdiffDim = k2Cart - k1Cart;
      for(int j : {0,1,2}) {
        kdiff(i, j) = (kdiffDim(j) == 0) ? 0.0 : 1./kdiffDim(j);
      }
      if(mpi->mpiHead()) std::cout << k1Cart.transpose() << " " << k2Cart.transpose() << " " << kdiffDim.transpose() << std::endl;
  }

    // for each direction of the finite differences
    for (int dk : {0, 1, 2}) {

      // find the points which are +2, +1, -1, -2 along this direction
      Eigen::Vector3d k3Plus = k3idx;
      Eigen::Vector3d k3Mins = k3idx;
      Eigen::Vector3d k3Plus2 = k3idx;
      Eigen::Vector3d k3Mins2 = k3idx;

      // add integer displacement, so that we have the forwards/backwards point indices
      k3Plus(dk) += 1;
      k3Mins(dk) -= 1;
      k3Plus2(dk) += 2;
      k3Mins2(dk) -= 2;

      // if one of the points goes out of bounds, wrap around the index
      if(k3Plus(dk) == kMesh(dk))             k3Plus(dk) = 0;
      if(k3Mins(dk) < 0)                      k3Mins(dk) = kMesh(dk) - 1;
      if(k3Plus2(dk) == kMesh(dk))            k3Plus2(dk) = 0;
      if(k3Plus2(dk) == kMesh(dk) + 1)        k3Plus2(dk) = 1;
      if(k3Mins2(dk) == -1)                   k3Mins2(dk) = kMesh(dk) - 1;
      if(k3Mins2(dk) == -2)                   k3Mins2(dk) = kMesh(dk) - 2;

      // convert these points back into crystal coordinates
      for (int i : {0, 1, 2}) {
        k3Plus(i) /= float(kMesh(i));
        k3Mins(i) /= float(kMesh(i));
        k3Plus2(i) /= float(kMesh(i));
        k3Mins2(i) /= float(kMesh(i));
      }

      // find kindex of each point on the points mesh
      int ikPlusRed = outerBandStructure.getPointIndex(k3Plus, true);
      int ikMinsRed = outerBandStructure.getPointIndex(k3Mins, true);
      int ikPlusRed2 = outerBandStructure.getPointIndex(k3Plus2, true);
      int ikMinsRed2 = outerBandStructure.getPointIndex(k3Mins2, true);

      // if symmetries are turned on, we need to map the global kmesh
      // index to a irreducible point, as this is what will be used
      // by the interalDiagonal element
      int ikPlus = points.asIrreducibleIndex(ikPlusRed);
      int ikMins = points.asIrreducibleIndex(ikMinsRed);
      int ikPlus2 = points.asIrreducibleIndex(ikPlusRed2);
      int ikMins2 = points.asIrreducibleIndex(ikMinsRed2);

      // convert these to wavevector index objects
      WavevectorIndex ikPlusIdx(ikPlus);
      WavevectorIndex ikMinsIdx(ikMins);
      WavevectorIndex ikPlusIdx2(ikPlus2);
      WavevectorIndex ikMinsIdx2(ikMins2);

      // get state indices
      int isPlus = outerBandStructure.getIndex(ikPlusIdx, ibIdx);
      int isMins = outerBandStructure.getIndex(ikMinsIdx, ibIdx);
      int isPlus2 = outerBandStructure.getIndex(ikPlusIdx2, ibIdx);
      int isMins2 = outerBandStructure.getIndex(ikMinsIdx2, ibIdx);

      // get vector bte indices
      int ibte = outerBandStructure.stateToBte(isIdx).get();
      StateIndex istatePlus(isPlus);
      StateIndex istateMins(isMins);
      int ibtePlus = outerBandStructure.stateToBte(istatePlus).get();
      int ibteMins = outerBandStructure.stateToBte(istateMins).get();
      StateIndex istatePlus2(isPlus2);
      StateIndex istateMins2(isMins2);
      int ibtePlus2 = outerBandStructure.stateToBte(istatePlus2).get();
      int ibteMins2 = outerBandStructure.stateToBte(istateMins2).get();

      // calculate the corrections for this state
      Eigen::Vector3d vel = outerBandStructure.getGroupVelocity(isIdx);
      Eigen::Vector3d y = -1 * vel.cross(magneticField);
      Eigen::Vector3d correction ={ y(0)*((2*kdiff(dk, 0))), y(1)*((2*kdiff(dk,1))), y(2)*((2*kdiff(dk,2)))};

      for(int iCalc = 0; iCalc < numCalculations; iCalc++) {
        for(int i : {0, 1, 2}) {
          nEB(iCalc, i, ibte) = correction(i) * (n[0](iCalc, i, ibtePlus) - n[0](iCalc, i, ibteMins));
          nTB(iCalc, i, ibte) = correction(i) * (n[1](iCalc, i, ibtePlus) - n[1](iCalc, i, ibteMins));
        }
      }
    }
  }
  std::vector<VectorBTE> nReturn = {nEB,nTB};
  return nReturn;
}

