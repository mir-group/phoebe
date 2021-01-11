#include "ph_scattering.h"
#include "constants.h"
#include "helper_3rd_state.h"
#include "io.h"
#include "mpiHelper.h"
#include "periodic_table.h"
#include <cmath>

PhScatteringMatrix::PhScatteringMatrix(Context &context_,
                                       StatisticsSweep &statisticsSweep_,
                                       BaseBandStructure &innerBandStructure_,
                                       BaseBandStructure &outerBandStructure_,
                                       Interaction3Ph *coupling3Ph_,
                                       PhononH0 *h0_)
    : ScatteringMatrix(context_, statisticsSweep_, innerBandStructure_,
                       outerBandStructure_),
      coupling3Ph(coupling3Ph_), h0(h0_) {
  if (&innerBandStructure != &outerBandStructure && h0 == nullptr) {
    Error e("PhScatteringMatrix needs h0 for incommensurate grids");
  }

  // setup here the isotopic scattering
  if (context.getWithIsotopeScattering()) {
    auto crystal = outerBandStructure.getPoints().getCrystal();
    int numAtoms = crystal.getNumAtoms();

    // create vector with the interaction strength
    massVariance = Eigen::VectorXd::Zero(numAtoms);

    // load the mass variance at natural abundances. Hard coded.
    PeriodicTable periodicTable;
    auto atomsNames = crystal.getAtomicNames();
    int i = 0;
    for (const auto &atomName : atomsNames) {
      double thisMass = periodicTable.getMass(atomName);
      // since the phonon eigenvectors are renormalized with sqrt(mass)
      // we add a correction factor in the coupling here
      massVariance(i) =
          thisMass * thisMass * periodicTable.getMassVariance(atomName);
      i += 1;
    }

    // check for user-defined mass variance
    auto userMassVariance = context.getMassVariance();
    if (userMassVariance.size() > 0) {
      massVariance = userMassVariance;
      if (massVariance.size() != numAtoms) {
        Error e("user mass variance should be set for each atom");
        // i.e. for each atom in the unit cell (not each species)
      }
    }

    doIsotopes = true;

  } else {
    doIsotopes = false;
  }

  doBoundary = false;
  boundaryLength = context.getBoundaryLength();
  if (!std::isnan(boundaryLength)) {
    if (boundaryLength > 0.) {
      doBoundary = true;
    }
  }
}

PhScatteringMatrix::PhScatteringMatrix(const PhScatteringMatrix &that)
    : ScatteringMatrix(that), coupling3Ph(that.coupling3Ph), h0(that.h0),
      massVariance(that.massVariance), doIsotopes(that.doIsotopes),
      boundaryLength(that.boundaryLength), doBoundary(that.doBoundary) {}

PhScatteringMatrix &
PhScatteringMatrix::operator=(const PhScatteringMatrix &that) {
  ScatteringMatrix::operator=(that);
  if (this != &that) {
    coupling3Ph = that.coupling3Ph;
    h0 = that.h0;
    massVariance = that.massVariance;
    doIsotopes = that.doIsotopes;
    boundaryLength = that.boundaryLength;
    doBoundary = that.doBoundary;
  }
  return *this;
}

// 3 cases:
// theMatrix and linewidth is passed: we compute and store in memory the
// scattering matrix and the diagonal
// inPopulation+outPopulation is passed: we compute the action of the
//       scattering matrix on the in vector, returning outVec = sMatrix*vector
// only linewidth is passed: we compute only the linewidths
void PhScatteringMatrix::builder(VectorBTE *linewidth,
                                 std::vector<VectorBTE> &inPopulations,
                                 std::vector<VectorBTE> &outPopulations) {
  // notes: + process is (1+2) -> 3
  //        - processes are (1+3)->2 and (3+2)->1

  const double energyCutoff = 0.001 / ryToCmm1; // discard states with small
  // energies (smaller than 0.001 cm^-1

  int switchCase = 0;
  if (theMatrix.rows() != 0 && linewidth != nullptr && inPopulations.empty() &&
      outPopulations.empty()) {
    switchCase = 0;
  } else if (theMatrix.rows() == 0 && linewidth == nullptr &&
             !inPopulations.empty() && !outPopulations.empty()) {
    switchCase = 1;
  } else if (theMatrix.rows() == 0 && linewidth != nullptr &&
             inPopulations.empty() && outPopulations.empty()) {
    switchCase = 2;
  } else {
    Error e("builder3Ph found a non-supported case");
  }

  if ((linewidth != nullptr) && (linewidth->dimensionality != 1)) {
    Error e("The linewidths shouldn't have dimensionality");
  }

  auto particle = outerBandStructure.getParticle();
  int numAtoms = innerBandStructure.getPoints().getCrystal().getNumAtoms();
  int numCalculations = statisticsSweep.getNumCalculations();

  // note: innerNumFullPoints is the number of points in the full grid
  // may be larger than innerNumPoints, when we use ActiveBandStructure
  double norm = 1. / context.getQMesh().prod();

  // precompute Bose populations
  VectorBTE outerBose(statisticsSweep, outerBandStructure, 1);
#pragma omp parallel for default(none) shared(mpi, particle, outerBose)
  for (int iBte : mpi->divideWorkIter(numStates)) {
    BteIndex iBteIdx(iBte);
    StateIndex isIdx = outerBandStructure.bteToState(iBteIdx);
    double energy = outerBandStructure.getEnergy(isIdx);
    for (int iCalc = 0; iCalc < statisticsSweep.getNumCalculations(); iCalc++) {
      double temperature = statisticsSweep.getCalcStatistics(iCalc).temperature;
      outerBose(iCalc, 0, iBte) = particle.getPopulation(energy, temperature);
    }
  }
  mpi->allReduceSum(&outerBose.data);
  VectorBTE innerBose(statisticsSweep, innerBandStructure, 1);
  if (&innerBandStructure == &outerBandStructure) {
    innerBose = outerBose;
  } else {
#pragma omp parallel for default(none) shared(mpi, particle, innerBose)
    for (int iBte : mpi->divideWorkIter(numStates)) {
      BteIndex iBteIdx(iBte);
      StateIndex isIdx = innerBandStructure.bteToState(iBteIdx);
      double energy = innerBandStructure.getEnergy(isIdx);
      for (int iCalc = 0; iCalc < statisticsSweep.getNumCalculations();
           iCalc++) {
        double temperature =
            statisticsSweep.getCalcStatistics(iCalc).temperature;
        innerBose(iCalc, 0, iBte) = particle.getPopulation(energy, temperature);
      }
    }
    mpi->allReduceSum(&innerBose.data);
  }

  std::vector<std::tuple<std::vector<int>, int>> qPairIterator =
      getIteratorWavevectorPairs(switchCase);

  Helper3rdState pointHelper(innerBandStructure, outerBandStructure, outerBose,
                             smearing->getType(), h0);
  LoopPrint loopPrint("computing scattering matrix", "q-points",
                      qPairIterator.size());

  /** Very important: the code must be executed with a loop over q2 outside
   * and a loop over q1 inside. This is because the 3-ph coupling must compute
   * the Fourier transform on q1 and q2. The 3-ph class splits the Fourier
   * transform in two parts, recycling the FT over q2 for more values of q1.
   * Thus, we have some speed when executing in this order.
   * PointHelper too assumes that order of loop execution.
   */
  // outer loop over q2
  for (auto tup : qPairIterator) {
    std::vector<int> iq1Indexes = std::get<0>(tup);
    int iq2 = std::get<1>(tup);
    WavevectorIndex iq2Index(iq2);

    Point q2Point = innerBandStructure.getPoint(iq2);
    Eigen::VectorXd energies2 = innerBandStructure.getEnergies(iq2Index);
    int nb2 = energies2.size();
    Eigen::MatrixXd v2s = innerBandStructure.getGroupVelocities(iq2Index);
    Eigen::Vector3d q2 = innerBandStructure.getWavevector(iq2Index);
    Eigen::MatrixXcd ev2 = innerBandStructure.getEigenvectors(iq2Index);

    int nq1 = iq1Indexes.size();

    auto t = innerBandStructure.getRotationToIrreducible(
        q2Point.getCoordinates(Points::cartesianCoordinates),
        Points::cartesianCoordinates);
    int iq2Irr = std::get<0>(t);
    Eigen::Matrix3d rotation = std::get<1>(t);
    // rotation such that qIrr = R * qRed

    loopPrint.update();
    pointHelper.prepare(iq1Indexes, iq2);

    // prepare batches based on memory usage
    int numBatches = coupling3Ph->estimateNumBatches(nq1, nb2);

    // precalculate D3cached for current value of q2
    coupling3Ph->cacheD3(q2);

    // loop over batches of q1s
    // later we will loop over the q1s inside each batch
    // this is done to optimize the usage and data transfer of a GPU
    for (int iBatch = 0; iBatch < numBatches; iBatch++) {
      // start and end point for current batch
      int start = nq1 * iBatch / numBatches;
      int end = nq1 * (iBatch + 1) / numBatches;
      int batch_size = end - start;

      std::vector<Eigen::Vector3d> q1_v(batch_size);
      std::vector<Eigen::MatrixXcd> ev1_v(batch_size);
      std::vector<Eigen::MatrixXcd> ev3Plus_v(batch_size);
      std::vector<Eigen::MatrixXcd> ev3Minus_v(batch_size);
      std::vector<int> nb1_v(batch_size);
      std::vector<int> nb3Plus_v(batch_size);
      std::vector<int> nb3Minus_v(batch_size);

      // do prep work for all values of q1 in current batch,
      // store stuff needed for couplings later
      for (int iq1Batch = 0; iq1Batch < batch_size; iq1Batch++) {
        int iq1 = iq1Indexes[start + iq1Batch];
        WavevectorIndex iq1Index(iq1);

        // note: for computing linewidths on a path, we must distinguish
        // that q1 and q2 are on different meshes, and that q3+/- may not
        // fall into known meshes and therefore needs to be computed

        Point q1Point = outerBandStructure.getPoint(iq1);
        Eigen::VectorXd energies1 = outerBandStructure.getEnergies(iq1Index);
        int nb1 = energies1.size();
        Eigen::MatrixXd v1s = outerBandStructure.getGroupVelocities(iq1Index);

        auto tup1 = pointHelper.get(q1Point, q2Point, Helper3rdState::casePlus);
        auto tup2 =
            pointHelper.get(q1Point, q2Point, Helper3rdState::caseMinus);

        auto energies3Plus = std::get<1>(tup1);
        auto ev3Plus = std::get<3>(tup1);
        auto v3ps = std::get<4>(tup1);
        auto bose3PlusData = std::get<5>(tup1);

        auto energies3Minus = std::get<1>(tup2);
        auto ev3Minus = std::get<3>(tup2);
        auto v3ms = std::get<4>(tup2);
        auto bose3MinusData = std::get<5>(tup2);

        q1_v[iq1Batch] = outerBandStructure.getWavevector(iq1Index);
        nb1_v[iq1Batch] = nb1;
        nb3Plus_v[iq1Batch] = energies3Plus.size();
        nb3Minus_v[iq1Batch] = energies3Minus.size();
        ev1_v[iq1Batch] = outerBandStructure.getEigenvectors(iq1Index);
        ev3Plus_v[iq1Batch] = ev3Plus;
        ev3Minus_v[iq1Batch] = ev3Minus;
      }

      // calculate batch of couplings
      auto tuple1 = coupling3Ph->getCouplingsSquared(
          q1_v, q2, ev1_v, ev2, ev3Plus_v, ev3Minus_v, nb1_v, nb2,
          nb3Plus_v, nb3Minus_v);
      auto couplingPlus_v = std::get<0>(tuple1);
      auto couplingMinus_v = std::get<1>(tuple1);

      // do postprocessing loop with batch of couplings
      for (int iq1Batch = 0; iq1Batch < batch_size; iq1Batch++) {
        int iq1 = iq1Indexes[start + iq1Batch];
        WavevectorIndex iq1Index(iq1);
        auto couplingPlus = couplingPlus_v[iq1Batch];
        auto couplingMinus = couplingMinus_v[iq1Batch];
        Point q1Point = outerBandStructure.getPoint(iq1);
        Eigen::VectorXd energies1 =
            outerBandStructure.getEnergies(iq1Index);
        int nb1 = energies1.size();
        Eigen::MatrixXd v1s = outerBandStructure.getGroupVelocities(iq1Index);

        auto tup1 = pointHelper.get(q1Point, q2Point, Helper3rdState::casePlus);
        auto tup2 =
            pointHelper.get(q1Point, q2Point, Helper3rdState::caseMinus);

        auto energies3Plus = std::get<1>(tup1);
        auto nb3Plus = energies3Plus.size();
        auto v3sPlus = std::get<4>(tup1);
        auto bose3PlusData = std::get<5>(tup1);

        auto energies3Minus = std::get<1>(tup2);
        auto nb3Minus = energies3Minus.size();
        auto v3sMinus = std::get<4>(tup2);
        auto bose3MinusData = std::get<5>(tup2);

        for (int ib1 = 0; ib1 < nb1; ib1++) {
          for (int ib2 = 0; ib2 < nb2; ib2++) {
            for (int ib3 = 0; ib3 < nb3Plus; ib3++) {

              double en1 = energies1(ib1);
              double en2 = energies2(ib2);
              double en3Plus = energies3Plus(ib3);
              if (en1 < energyCutoff || en2 < energyCutoff ||
                  en3Plus < energyCutoff) {
                continue;
              }
              double enProd = en1 * en2 * en3Plus;

              int is1 = outerBandStructure.getIndex(iq1Index, BandIndex(ib1));
              int is2 = innerBandStructure.getIndex(iq2Index, BandIndex(ib2));
              int is2Irr = innerBandStructure.getIndex(WavevectorIndex(iq2Irr),
                                                       BandIndex(ib2));
              StateIndex is1Idx(is1);
              StateIndex is2Idx(is2);
              StateIndex is2IrrIdx(is2Irr);
              BteIndex ind1Idx = outerBandStructure.stateToBte(is1Idx);
              BteIndex ind2Idx = innerBandStructure.stateToBte(is2IrrIdx);
              int ind1 = ind1Idx.get();
              int ind2 = ind2Idx.get();

              double deltaPlus;
              if (smearing->getType() == DeltaFunction::gaussian) {
                deltaPlus = smearing->getSmearing(en1 + en2 - en3Plus);
              } else if (smearing->getType() ==
                         DeltaFunction::adaptiveGaussian) {
                Eigen::Vector3d v = v2s.row(ib2) - v3sPlus.row(ib3);
                deltaPlus = smearing->getSmearing(en1 + en2 - en3Plus, v);
              } else {
                deltaPlus = smearing->getSmearing(en3Plus - en1, is2Idx);
              }

              if (deltaPlus < 0) {
                continue;
              }

              // loop on temperature
              for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
                double bose1 = outerBose(iCalc, 0, ind1);
                double bose2 = innerBose(iCalc, 0, ind2);
                double bose3Plus = bose3PlusData(iCalc, ib3);

                // Calculate transition probability W+
                double ratePlus = pi * 0.25 * bose1 * bose2 * (bose3Plus + 1.) *
                                  couplingPlus(ib1, ib2, ib3) * deltaPlus *
                                  norm / enProd;

                if (switchCase == 0) { // case of matrix construction
                  if (context.getUseSymmetries()) {
                    for (int i : {0, 1, 2}) {
                      for (int j : {0, 1, 2}) {
                        CartIndex iIndex(i);
                        CartIndex jIndex(j);
                        int iMat1 = getSMatrixIndex(ind1Idx, iIndex);
                        int iMat2 = getSMatrixIndex(ind2Idx, jIndex);
                        if (theMatrix.indicesAreLocal(iMat1, iMat2)) {
                          if (i == 0 && j == 0) {
                            linewidth->operator()(iCalc, 0, ind1) += ratePlus;
                          }
                          theMatrix(iMat1, iMat2) += rotation(i, j) * ratePlus;
                        }
                      }
                    }
                  } else {
                    linewidth->operator()(iCalc, 0, ind1) += ratePlus;
                    theMatrix(ind1, ind2) += ratePlus;
                  }

                } else if (switchCase ==
                           1) { // case of matrix-vector multiplication
                  // we build the scattering matrix A = S*n(n+1)
                  // here we rotate the populations from the irreducible point
                  for (unsigned int iInput = 0; iInput < inPopulations.size();
                       iInput++) {
                    Eigen::Vector3d inPopRot;
                    inPopRot.setZero();
                    for (int i : {0, 1, 2}) {
                      for (int j : {0, 1, 2}) {
                        inPopRot(i) += rotation(i, j) *
                                       inPopulations[iInput](iCalc, j, ind2);
                      }
                    }

                    for (int i : {0, 1, 2}) {
                      outPopulations[iInput](iCalc, i, ind1) +=
                          ratePlus * inPopRot(i);
                      outPopulations[iInput](iCalc, i, ind1) +=
                          ratePlus * inPopulations[iInput](iCalc, i, ind1);
                    }
                  }

                } else { // case of linewidth construction
                  linewidth->operator()(iCalc, 0, ind1) += ratePlus;
                }
              }
            }
          }
        }

        for (int ib1 = 0; ib1 < nb1; ib1++) {
          for (int ib2 = 0; ib2 < nb2; ib2++) {
            for (int ib3 = 0; ib3 < nb3Minus; ib3++) {
              double en1 = energies1(ib1);
              double en2 = energies2(ib2);
              double en3Minus = energies3Minus(ib3);
              if (en1 < energyCutoff || en2 < energyCutoff ||
                  en3Minus < energyCutoff) {
                continue;
              }
              double enProd = en1 * en2 * en3Minus;

              int is1 = outerBandStructure.getIndex(iq1Index, BandIndex(ib1));
              int is2 = innerBandStructure.getIndex(iq2Index, BandIndex(ib2));
              int is2Irr = innerBandStructure.getIndex(WavevectorIndex(iq2Irr),
                                                       BandIndex(ib2));
              StateIndex is1Idx(is1);
              StateIndex is2Idx(is2);
              StateIndex is2IrrIdx(is2Irr);
              BteIndex ind1Idx = outerBandStructure.stateToBte(is1Idx);
              BteIndex ind2Idx = innerBandStructure.stateToBte(is2IrrIdx);
              int ind1 = ind1Idx.get();
              int ind2 = ind2Idx.get();

              double deltaMinus1, deltaMinus2;
              if (smearing->getType() == DeltaFunction::gaussian) {
                deltaMinus1 = smearing->getSmearing(en1 + en3Minus - en2);
                deltaMinus2 = smearing->getSmearing(en2 + en3Minus - en1);
              } else if (smearing->getType() ==
                         DeltaFunction::adaptiveGaussian) {
                Eigen::Vector3d v = v2s.row(ib2) - v3sMinus.row(ib3);
                deltaMinus1 = smearing->getSmearing(en1 + en3Minus - en2, v);
                deltaMinus2 = smearing->getSmearing(en2 + en3Minus - en1, v);
              } else {
                // Note: here I require inner == outer band structure
                deltaMinus1 = smearing->getSmearing(en1 + en3Minus, is2Idx);
                deltaMinus2 = smearing->getSmearing(en1 - en3Minus, is2Idx);
              }

              if (deltaMinus1 < 0. && deltaMinus2 < 0.)
                continue;
              if (deltaMinus1 < 0.)
                deltaMinus1 = 0.;
              if (deltaMinus2 < 0.)
                deltaMinus2 = 0.;

              for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
                double bose1 = outerBose(iCalc, 0, ind1);
                double bose2 = innerBose(iCalc, 0, ind2);
                double bose3Minus = bose3MinusData(iCalc, ib3);

                // Calculate transition probability W-
                double rateMinus1 =
                    pi * 0.25 * bose3Minus * bose1 * (bose2 + 1.) *
                    couplingMinus(ib1, ib2, ib3) * deltaMinus1 * norm / enProd;
                double rateMinus2 =
                    pi * 0.25 * bose2 * bose3Minus * (bose1 + 1.) *
                    couplingMinus(ib1, ib2, ib3) * deltaMinus2 * norm / enProd;

                if (switchCase == 0) { // case of matrix construction
                  if (context.getUseSymmetries()) {
                    for (int i : {0, 1, 2}) {
                      for (int j : {0, 1, 2}) {
                        CartIndex iIndex(i);
                        CartIndex jIndex(j);
                        int iMat1 = getSMatrixIndex(ind1Idx, iIndex);
                        int iMat2 = getSMatrixIndex(ind2Idx, jIndex);
                        if (theMatrix.indicesAreLocal(iMat1, iMat2)) {
                          if (i == 0 && j == 0) {
                            linewidth->operator()(iCalc, 0, ind1) +=
                                0.5 * rateMinus2;
                          }
                          theMatrix(iMat1, iMat2) -=
                              rotation(i, j) * (rateMinus1 + rateMinus2);
                        }
                      }
                    }
                  } else {
                    linewidth->operator()(iCalc, 0, ind1) += 0.5 * rateMinus2;
                    theMatrix(ind1, ind2) -= rateMinus1 + rateMinus2;
                  }

                } else if (switchCase == 1) { // matrix-vector multiplication
                  for (unsigned int iInput = 0; iInput < inPopulations.size();
                       iInput++) {
                    Eigen::Vector3d inPopRot;
                    inPopRot.setZero();
                    for (int i : {0, 1, 2}) {
                      for (int j : {0, 1, 2}) {
                        inPopRot(i) += rotation(i, j) *
                                       inPopulations[iInput](iCalc, j, ind2);
                      }
                    }

                    for (int i : {0, 1, 2}) {
                      outPopulations[iInput](iCalc, i, ind1) -=
                          (rateMinus1 + rateMinus2) * inPopRot(i);
                      outPopulations[iInput](iCalc, i, ind1) +=
                          0.5 * rateMinus2 *
                          inPopulations[iInput](iCalc, i, ind1);
                    }
                  }
                } else { // case of linewidth construction
                  linewidth->operator()(iCalc, 0, ind1) += 0.5 * rateMinus2;
                }
              }
            }
          }
        }
      }
    }
  }

  // Isotope scattering
  if (doIsotopes) {
    // for (auto [iq1, iq2] : qPairIterator) {
    for (auto tup : qPairIterator) {
      auto iq1Indexes = std::get<0>(tup);
      int iq2 = std::get<1>(tup);

      WavevectorIndex iq2Index(iq2);
      Eigen::VectorXd state2Energies = innerBandStructure.getEnergies(iq2Index);
      int nb2 = state2Energies.size();
      Eigen::Tensor<std::complex<double>, 3> ev2 =
          innerBandStructure.getPhEigenvectors(iq2Index);
      Eigen::MatrixXd v2s = innerBandStructure.getGroupVelocities(iq2Index);

      auto q2Coords = innerBandStructure.getPoint(iq2).getCoordinates(
          Points::cartesianCoordinates);

      auto t = innerBandStructure.getRotationToIrreducible(
          q2Coords, Points::cartesianCoordinates);
      // rotation such that
      int iq2Irr = std::get<0>(t);
      Eigen::Matrix3d rotation = std::get<1>(t);

      for (auto iq1 : iq1Indexes) {
        WavevectorIndex iq1Index(iq1);

        // note: for computing linewidths on a path, we must distinguish
        // that q1 and q2 are on different meshes, and that q3+/- may not
        // fall into known meshes and therefore needs to be computed

        Eigen::VectorXd state1Energies =
            outerBandStructure.getEnergies(iq1Index);
        int nb1 = state1Energies.size();
        Eigen::Tensor<std::complex<double>, 3> ev1 =
            outerBandStructure.getPhEigenvectors(iq1Index);
        Eigen::MatrixXd v1s = outerBandStructure.getGroupVelocities(iq1Index);

        for (int ib1 = 0; ib1 < nb1; ib1++) {
          double en1 = state1Energies(ib1);

          int is1 = outerBandStructure.getIndex(iq1Index, BandIndex(ib1));
          StateIndex is1Idx(is1);
          int ind1 = outerBandStructure.stateToBte(is1Idx).get();

          if (en1 < energyCutoff) {
            continue;
          }

          for (int ib2 = 0; ib2 < nb2; ib2++) {
            double en2 = state2Energies(ib2);
            int is2 = innerBandStructure.getIndex(WavevectorIndex(iq2Irr),
                                                  BandIndex(ib2));
            StateIndex is2Idx(is2);
            int ind2 = innerBandStructure.stateToBte(is2Idx).get();

            if (en2 < energyCutoff) {
              continue;
            }

            double deltaIso;
            if (smearing->getType() == DeltaFunction::gaussian) {
              deltaIso = smearing->getSmearing(en1 - en2);
            } else if (smearing->getType() == DeltaFunction::adaptiveGaussian) {
              deltaIso = smearing->getSmearing(en1 - en2, v2s.row(ib2));
              deltaIso = smearing->getSmearing(en1 - en2, v1s.row(ib1));
              deltaIso *= 0.5;
            } else {
              deltaIso = smearing->getSmearing(en1, is2Idx);
            }

            double termIso = 0.;
            for (int iat = 0; iat < numAtoms; iat++) {
              std::complex<double> zzIso = complexZero;
              for (int kDim : {0, 1, 2}) { // cartesian indices
                zzIso += std::conj(ev1(kDim, iat, ib1)) * ev2(kDim, iat, ib2);
              }
              termIso += std::norm(zzIso) * massVariance(iat);
            }
            termIso *= pi * 0.5 * norm * en1 * en2 * deltaIso;

            for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
              double bose1 = outerBose(iCalc, 0, ind1);
              double bose2 = innerBose(iCalc, 0, ind2);

              double rateIso =
                  termIso * (bose1 * bose2 + 0.5 * (bose1 + bose2));

              if (switchCase == 0) { // case of matrix construction
                if (context.getUseSymmetries()) {
                  for (int i : {0, 1, 2}) {
                    for (int j : {0, 1, 2}) {
                      BteIndex iBte1(ind1);
                      BteIndex iBte2(ind2);
                      CartIndex iIndex(i);
                      CartIndex jIndex(j);
                      int iMat1 = getSMatrixIndex(iBte1, iIndex);
                      int iMat2 = getSMatrixIndex(iBte2, jIndex);
                      if (theMatrix.indicesAreLocal(iMat1, iMat2)) {
                        if (i == 0 && j == 0) {
                          linewidth->operator()(iCalc, 0, ind1) += rateIso;
                        }
                        theMatrix(iMat1, iMat2) += rotation(i, j) * rateIso;
                      }
                    }
                  }
                } else {
                  linewidth->operator()(iCalc, 0, ind1) += rateIso;
                  theMatrix(ind1, ind2) += rateIso;
                }

              } else if (switchCase ==
                         1) { // case of matrix-vector multiplication
                for (unsigned int iInput = 0; iInput < inPopulations.size();
                     iInput++) {

                  // here we rotate the populations from the irreducible point
                  Eigen::Vector3d inPopRot;
                  inPopRot.setZero();
                  for (int i = 0; i < dimensionality_; i++) {
                    for (int j = 0; j < dimensionality_; j++) {
                      inPopRot(i) += rotation(i, j) *
                                     inPopulations[iInput](iCalc, j, ind2);
                    }
                  }
                  for (int i = 0; i < dimensionality_; i++) {
                    outPopulations[iInput](iCalc, i, ind1) +=
                        rateIso * inPopRot(i);
                    outPopulations[iInput](iCalc, i, ind1) +=
                        rateIso * inPopulations[iInput](iCalc, i, ind1);
                  }
                }

              } else { // case of linewidth construction
                linewidth->operator()(iCalc, 0, ind1) += rateIso;
              }
            }
          }
        }
      }
    }
  }

  if (switchCase == 1) {
    for (unsigned int i = 0; i < outPopulations.size(); i++) {
      mpi->allReduceSum(&outPopulations[i].data);
    }
  } else {
    mpi->allReduceSum(&linewidth->data);
  }
  // I prefer to close loopPrint after the MPI barrier: all MPI are synced here
  loopPrint.close();

  // Add boundary scattering

  if (doBoundary) {
#pragma omp parallel for default(none) shared(                                 \
    inPopulations, outPopulations, innerBose, particle, linewidth, switchCase)
    for (int is1 : outerBandStructure.irrStateIterator()) { // in serial!
      StateIndex is1Idx(is1);
      double energy = outerBandStructure.getEnergy(is1Idx);
      auto vel = outerBandStructure.getGroupVelocity(is1Idx);
      int ind1 = outerBandStructure.stateToBte(is1Idx).get();

      for (int iCalc = 0; iCalc < statisticsSweep.getNumCalculations();
           iCalc++) {
        double temperature =
            statisticsSweep.getCalcStatistics(iCalc).temperature;
        // n(n+1)
        double termPop = particle.getPopPopPm1(energy, temperature);
        double rate = vel.squaredNorm() / boundaryLength * termPop;

        if (switchCase == 0) { // case of matrix construction
          linewidth->operator()(iCalc, 0, ind1) += rate;

        } else if (switchCase == 1) { // case of matrix-vector multiplication
          for (unsigned int iVec = 0; iVec < inPopulations.size(); iVec++) {
            for (int i : {0, 1, 2}) {
              outPopulations[iVec](iCalc, i, ind1) +=
                  rate * inPopulations[iVec](iCalc, i, ind1);
            }
          }

        } else { // case of linewidth construction
          linewidth->operator()(iCalc, 0, ind1) += rate;
        }
      }
    }
  }

  // some phonons like acoustic modes at the gamma, with omega = 0,
  // might have zero frequencies, and infinite populations. We set those
  // matrix elements to zero.
  if (switchCase == 0) {
    // case of matrix construction
    if (context.getUseSymmetries()) {
      for (auto iBte1 : excludeIndices) {
        linewidth->data.col(iBte1).setZero();
        for (auto iBte2 : excludeIndices) {
          for (int i : {0, 1, 2}) {
            for (int j : {0, 1, 2}) {
              BteIndex iBte1Idx(iBte1);
              BteIndex iBte2Idx(iBte2);
              CartIndex iCart1(i);
              CartIndex iCart2(j);
              int iMat1 = getSMatrixIndex(iBte1Idx, iCart1);
              int iMat2 = getSMatrixIndex(iBte2Idx, iCart2);
              theMatrix(iMat1, iMat2) = 0.;
            }
          }
        }
      }
    } else {
      for (auto iBte1 : excludeIndices) {
        linewidth->data.col(iBte1).setZero();
        for (auto iBte2 : excludeIndices) {
          theMatrix(iBte1, iBte2) = 0.;
        }
      }
    }
  } else if (switchCase == 1) {
    // case of matrix-vector multiplication
    for (auto iBte1 : excludeIndices) {
      for (unsigned int i = 0; i < outPopulations.size(); i++) {
        outPopulations[i].data.col(iBte1).setZero();
      }
    }

  } else if (switchCase == 2) {
    // case of linewidth construction
    for (auto iBte1 : excludeIndices) {
      linewidth->data.col(iBte1).setZero();
    }
  }

  // we place the linewidths back in the diagonal of the scattering matrix
  // this because we may need an MPI_allReduce on the linewidths
  if (switchCase == 0) { // case of matrix construction
    int iCalc = 0;
    if (context.getUseSymmetries()) {
      for (int iBte = 0; iBte < numStates; iBte++) {
        BteIndex iBteIdx(iBte);
        for (int i : {0, 1, 2}) {
          CartIndex iCart(i);
          int iMat1 = getSMatrixIndex(iBteIdx, iCart);
          theMatrix(iMat1, iMat1) = linewidth->operator()(iCalc, 0, iBte);
        }
      }
    } else {
      for (int is = 0; is < numStates; is++) {
        theMatrix(is, is) = linewidth->operator()(iCalc, 0, is);
      }
    }
  }
}
