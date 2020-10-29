#include "ph_scattering.h"

#include "constants.h"
#include "helper_3rd_state.h"
#include "io.h"
#include "mpiHelper.h"
#include "periodic_table.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

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
    long i = 0;
    for (auto atomName : atomsNames) {
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
// theMatrix and linedith is passed: we compute and store in memory the scatt
//       matrix and the diagonal
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

  long numAtoms = innerBandStructure.getPoints().getCrystal().getNumAtoms();
  int numCalcs = statisticsSweep.getNumCalcs();

  // note: innerNumFullPoints is the number of points in the full grid
  // may be larger than innerNumPoints, when we use ActiveBandStructure
  long innerNumFullPoints = innerBandStructure.getNumPoints(true);

  // precompute Bose populations
  VectorBTE outerBose(statisticsSweep, outerBandStructure, 1);
  for (long is : mpi->divideWorkIter(outerBandStructure.getNumStates())) {
    double energy = outerBandStructure.getEnergy(is);
    for (long iCalc = 0; iCalc < statisticsSweep.getNumCalcs(); iCalc++) {
      double temperature = statisticsSweep.getCalcStatistics(iCalc).temperature;
      outerBose(iCalc, 0, is) = particle.getPopulation(energy, temperature);
    }
  }
  mpi->allReduceSum(&outerBose.data);
  VectorBTE innerBose(statisticsSweep, innerBandStructure, 1);
  if (&innerBandStructure == &outerBandStructure) {
    innerBose = outerBose;
  } else {
    for (long is : mpi->divideWorkIter(innerBandStructure.getNumStates())) {
      double energy = innerBandStructure.getEnergy(is);
      for (long iCalc = 0; iCalc < statisticsSweep.getNumCalcs(); iCalc++) {
        double temperature =
            statisticsSweep.getCalcStatistics(iCalc).temperature;
        innerBose(iCalc, 0, is) = particle.getPopulation(energy, temperature);
      }
    }
    mpi->allReduceSum(&innerBose.data);
  }

  std::vector<std::tuple<std::vector<long>, long>> qPairIterator =
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
    auto iq1Indexes = std::get<0>(tup);
    auto iq2 = std::get<1>(tup);
    auto iq2Index = WavevectorIndex(iq2);

    Point q2 = innerBandStructure.getPoint(iq2);
    Eigen::VectorXd state2Energies = innerBandStructure.getEnergies(iq2Index);
    int nb2 = state2Energies.size();
    Eigen::MatrixXd v2s = innerBandStructure.getGroupVelocities(iq2Index);
    Eigen::Vector3d q2_e = innerBandStructure.getWavevector(iq2Index);
    Eigen::MatrixXcd ev2_e = innerBandStructure.getEigenvectors(iq2Index);

    int nq1 = iq1Indexes.size();

    loopPrint.update();
    pointHelper.prepare(iq1Indexes, iq2);

    // prepare batches based on memory usage
    int numBatches = coupling3Ph->estimateNumBatches(nq1, nb2);

    // precalculate D3cached for current value of q2
    coupling3Ph->cacheD3(q2_e);

    // loop over batches of q1s
    // later we will loop over the q1s inside each batch
    // this is done to optimize the usage and data transfer of a GPU
    for (int ib = 0; ib < numBatches; ib++) {
      // start and end point for current batch
      int start = nq1 * ib / numBatches;
      int end = nq1 * (ib + 1) / numBatches;
      int batch_size = end - start;

      std::vector<Eigen::Vector3d> q1s_e(batch_size);
      std::vector<Eigen::MatrixXcd> ev1s_e(batch_size);
      std::vector<Eigen::MatrixXcd> ev3Pluss_e(batch_size);
      std::vector<Eigen::MatrixXcd> ev3Minss_e(batch_size);
      std::vector<int> nb1s_e(batch_size);
      std::vector<int> nb3Pluss_e(batch_size);
      std::vector<int> nb3Minss_e(batch_size);

      // do prep work for all values of q1 in current batch,
      // store stuff needed for couplings later
      for (int i = 0; i < batch_size; i++) {
        auto iq1 = iq1Indexes[start + i];
        auto iq1Index = WavevectorIndex(iq1);

        // note: for computing linewidths on a path, we must distinguish
        // that q1 and q2 are on different meshes, and that q3+/- may not
        // fall into known meshes and therefore needs to be computed

        Point q1 = outerBandStructure.getPoint(iq1);
        Eigen::VectorXd state1Energies = outerBandStructure.getEnergies(iq1Index);
        int nb1 = state1Energies.size();
        Eigen::MatrixXd v1s = outerBandStructure.getGroupVelocities(iq1Index);

        auto tup1 = pointHelper.get(q1, q2, Helper3rdState::casePlus);
        auto tup2 = pointHelper.get(q1, q2, Helper3rdState::caseMins);

        auto state3PlusEnergies = std::get<1>(tup1);
        auto eigvecs3Plus = std::get<3>(tup1);
        auto v3ps = std::get<4>(tup1);
        auto bose3PlusData = std::get<5>(tup1);

        auto state3MinsEnergies = std::get<1>(tup2);
        auto eigvecs3Mins = std::get<3>(tup2);
        auto v3ms = std::get<4>(tup2);
        auto bose3MinsData = std::get<5>(tup2);

        q1s_e[i] = outerBandStructure.getWavevector(iq1Index);
        nb1s_e[i] = nb1;
        nb3Pluss_e[i] = state3PlusEnergies.size();
        nb3Minss_e[i] = state3MinsEnergies.size();
        ev1s_e[i] = outerBandStructure.getEigenvectors(iq1Index);
        ev3Pluss_e[i] = eigvecs3Plus;
        ev3Minss_e[i] = eigvecs3Mins;
      }

      // calculate batch of couplings
      auto tup = coupling3Ph->getCouplingsSquared(
          q1s_e, q2_e, ev1s_e, ev2_e, ev3Pluss_e, ev3Minss_e, nb1s_e, nb2,
          nb3Pluss_e, nb3Minss_e);
      auto couplingPluss = std::get<0>(tup);
      auto couplingMinss = std::get<1>(tup);

      // do postprocessing loop with batch of couplings
      for (int i = 0; i < batch_size; i++) {
        auto iq1 = iq1Indexes[start + i];
        auto iq1Index = WavevectorIndex(iq1);
        auto couplingPlus = couplingPluss[i];
        auto couplingMins = couplingMinss[i];
        Point q1 = outerBandStructure.getPoint(iq1);
        Eigen::VectorXd state1Energies = outerBandStructure.getEnergies(iq1Index);
        int nb1 = state1Energies.size();
        Eigen::MatrixXd v1s = outerBandStructure.getGroupVelocities(iq1Index);

        auto tup1 = pointHelper.get(q1, q2, Helper3rdState::casePlus);
        auto tup2 = pointHelper.get(q1, q2, Helper3rdState::caseMins);

        auto state3PlusEnergies = std::get<1>(tup1);
        auto nb3Plus = std::get<2>(tup1);
        auto eigvecs3Plus = std::get<3>(tup1);
        auto v3ps = std::get<4>(tup1);
        auto bose3PlusData = std::get<5>(tup1);

        auto state3MinsEnergies = std::get<1>(tup2);
        auto nb3Mins = std::get<2>(tup2);
        auto eigvecs3Mins = std::get<3>(tup2);
        auto v3ms = std::get<4>(tup2);
        auto bose3MinsData = std::get<5>(tup2);

#pragma omp parallel for
        for (long ibbb = 0; ibbb < nb1 * nb2 * nb3Plus; ibbb++) {
          auto tup = decompress3Indeces(ibbb, nb1, nb2, nb3Plus);
          auto ib1 = std::get<0>(tup);
          auto ib2 = std::get<1>(tup);
          auto ib3 = std::get<2>(tup);

          double en1 = state1Energies(ib1);
          double en2 = state2Energies(ib2);
          double en3Plus = state3PlusEnergies(ib3);
          if (en1 < energyCutoff || en2 < energyCutoff ||
              en3Plus < energyCutoff) {
            continue;
          }
          double enProd = en1 * en2 * en3Plus;

          int ind1 =
              outerBandStructure.getIndex(WavevectorIndex(iq1), BandIndex(ib1));
          int ind2 =
              innerBandStructure.getIndex(WavevectorIndex(iq2), BandIndex(ib2));

          if (switchCase == 0) {
            // note: above we are parallelizing over wavevectors.
            // (for convenience of computing coupling3ph)
            // Not the same way as Matrix() is parallelized.
            // here we check that we don't duplicate efforts
            if (!theMatrix.indecesAreLocal(ind1, ind2)) {
              continue;
            }
          }

          double deltaPlus;
          switch (smearing->getType()) {
          case (DeltaFunction::gaussian):
            deltaPlus = smearing->getSmearing(en1 + en2 - en3Plus);
            break;
          case (DeltaFunction::adaptiveGaussian): {
            Eigen::Vector3d v = v2s.row(ib2) - v3ps.row(ib3);
            deltaPlus = smearing->getSmearing(en1 + en2 - en3Plus, v);
          } break;
          default:
            deltaPlus = smearing->getSmearing(en3Plus - en1, iq2, ib2);
            break;
          }

          if (deltaPlus < 0)
            continue;

          // loop on temperature
          for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
            double bose1 = outerBose(iCalc, 0, ind1);
            double bose2 = innerBose(iCalc, 0, ind2);
            double bose3Plus = bose3PlusData(iCalc, ib3);

            // Calculate transition probability W+
            double ratePlus = pi * 0.25 * bose1 * bose2 * (bose3Plus + 1.) *
                              couplingPlus(ib1, ib2, ib3) * deltaPlus /
                              innerNumFullPoints / enProd;

            switch (switchCase) {
            case (0):
              // case of matrix construction
              // we build the scattering matrix S
              theMatrix(ind1, ind2) += ratePlus;
              linewidth->operator()(iCalc, 0, ind1) += ratePlus;
              break;
            case (1):
              // case of matrix-vector multiplication
              // we build the scattering matrix A = S*n(n+1)
              for (unsigned int iInput = 0; iInput < inPopulations.size(); iInput++) {
                for (int i = 0; i < dimensionality_; i++) {
                  outPopulations[iInput](iCalc, i, ind1) +=
                      ratePlus *
                      inPopulations[iInput](iCalc, i, ind2);
                  outPopulations[iInput](iCalc, i, ind1) +=
                      ratePlus *
                      inPopulations[iInput](iCalc, i, ind1);
                }
              }
              break;
            case (2):
              // case of linewidth construction
              linewidth->operator()(iCalc, 0, ind1) += ratePlus;
              break;
            }
          }
        }
#pragma omp parallel for
        for (long ibbb = 0; ibbb < nb1 * nb2 * nb3Mins; ibbb++) {
          auto tup = decompress3Indeces(ibbb, nb1, nb2, nb3Mins);
          auto ib1 = std::get<0>(tup);
          auto ib2 = std::get<1>(tup);
          auto ib3 = std::get<2>(tup);

          double en1 = state1Energies(ib1);
          double en2 = state2Energies(ib2);
          double en3Mins = state3MinsEnergies(ib3);
          if (en1 < energyCutoff || en2 < energyCutoff ||
              en3Mins < energyCutoff) {
            continue;
          }
          double enProd = en1 * en2 * en3Mins;

          int ind1 =
              outerBandStructure.getIndex(WavevectorIndex(iq1), BandIndex(ib1));
          int ind2 =
              innerBandStructure.getIndex(WavevectorIndex(iq2), BandIndex(ib2));

          if (switchCase == 0) {
            // note: above we are parallelizing over wavevectors.
            // (for convenience of computing coupling3ph)
            // Not the same way as Matrix() is parallelized.
            // here we check that we don't duplicate efforts
            if (!theMatrix.indecesAreLocal(ind1, ind2)) {
              continue;
            }
          }

          double deltaMins1, deltaMins2;
          switch (smearing->getType()) {
          case (DeltaFunction::gaussian):
            deltaMins1 = smearing->getSmearing(en1 + en3Mins - en2);
            deltaMins2 = smearing->getSmearing(en2 + en3Mins - en1);
            break;
          case (DeltaFunction::adaptiveGaussian): {
            Eigen::Vector3d v = v2s.row(ib2) - v3ms.row(ib3);
            deltaMins1 = smearing->getSmearing(en1 + en3Mins - en2, v);
            deltaMins2 = smearing->getSmearing(en2 + en3Mins - en1, v);
          } break;
          default:
            deltaMins1 = smearing->getSmearing(en2 - en3Mins, iq1, ib1);
            deltaMins2 = smearing->getSmearing(en1 - en3Mins, iq2, ib2);
            break;
          }

          if (deltaMins1 < 0. && deltaMins2 < 0.)
            continue;
          if (deltaMins1 < 0.)
            deltaMins1 = 0.;
          if (deltaMins2 < 0.)
            deltaMins2 = 0.;

          for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
            double bose1 = outerBose(iCalc, 0, ind1);
            double bose2 = innerBose(iCalc, 0, ind2);
            double bose3Mins = bose3MinsData(iCalc, ib3);

            // Calculatate transition probability W-
            double rateMins1 = pi * 0.25 * bose3Mins * bose1 * (bose2 + 1.) *
                               couplingMins(ib1, ib2, ib3) * deltaMins1 /
                               innerNumFullPoints / enProd;
            double rateMins2 = pi * 0.25 * bose2 * bose3Mins * (bose1 + 1.) *
                               couplingMins(ib1, ib2, ib3) * deltaMins2 /
                               innerNumFullPoints / enProd;

            switch (switchCase) {
            case (0):
              // case of matrix construction
              theMatrix(ind1, ind2) -= rateMins1 + rateMins2;
              linewidth->operator()(iCalc, 0, ind1) += 0.5 * rateMins2;
              break;
            case (1):
              // case of matrix-vector multiplication
              for (unsigned int iInput = 0; iInput < inPopulations.size(); iInput++) {
                for (int i = 0; i < dimensionality_; i++) {
                  outPopulations[iInput](iCalc, i, ind1) -=
                      (rateMins1 + rateMins2) *
                      inPopulations[iInput](iCalc, i, ind2);
                  outPopulations[iInput](iCalc, i, ind1) +=
                      0.5 * rateMins2 *
                      inPopulations[iInput](iCalc, i, ind1);
                }
              }
              break;
            case (2):
              // case of linewidth construction
              linewidth->operator()(iCalc, 0, ind1) += 0.5 * rateMins2;
              break;
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
      auto iq2 = std::get<1>(tup);
#pragma omp parallel for
      for (auto iq1 : iq1Indexes) {
        auto iq1Index = WavevectorIndex(iq1);
        auto iq2Index = WavevectorIndex(iq2);
        Eigen::VectorXd state2Energies =
            innerBandStructure.getEnergies(iq2Index);
        int nb2 = state2Energies.size();

        Eigen::Tensor<std::complex<double>, 3> ev2 =
            innerBandStructure.getPhEigenvectors(iq2Index);

        // note: for computing linewidths on a path, we must distinguish
        // that q1 and q2 are on different meshes, and that q3+/- may not
        // fall into known meshes and therefore needs to be computed

        Eigen::VectorXd state1Energies = outerBandStructure.getEnergies(iq1Index);
        int nb1 = state1Energies.size();

        Eigen::Tensor<std::complex<double>, 3> ev1 =
            outerBandStructure.getPhEigenvectors(iq1Index);

        Eigen::MatrixXd v1s = outerBandStructure.getGroupVelocities(iq1Index);
        Eigen::MatrixXd v2s = innerBandStructure.getGroupVelocities(iq2Index);

        for (int ib1 = 0; ib1 < nb1; ib1++) {
          double en1 = state1Energies(ib1);
          int ind1 =
              outerBandStructure.getIndex(WavevectorIndex(iq1), BandIndex(ib1));
          if (en1 < energyCutoff) {
            continue;
          }

          for (int ib2 = 0; ib2 < nb2; ib2++) {
            double en2 = state2Energies(ib2);
            int ind2 = innerBandStructure.getIndex(WavevectorIndex(iq2),
                                                   BandIndex(ib2));
            if (en2 < energyCutoff) {
              continue;
            }

            if (switchCase == 0) {
              // note: above we are parallelizing over wavevectors.
              // (for convenience of computing coupling3ph)
              // Not the same way as Matrix() is parallelized.
              // here we check that we don't duplicate efforts
              if (!theMatrix.indecesAreLocal(ind1, ind2)) {
                continue;
              }
            }

            double deltaIso;
            switch (smearing->getType()) {
            case (DeltaFunction::gaussian):
              deltaIso = smearing->getSmearing(en1 - en2);
              break;
            case (DeltaFunction::adaptiveGaussian):
              deltaIso = smearing->getSmearing(en1 - en2, v2s.row(ib2));
              deltaIso = smearing->getSmearing(en1 - en2, v1s.row(ib1));
              deltaIso *= 0.5;
              break;
            default:
              deltaIso = smearing->getSmearing(en2, iq2, ib2);
              deltaIso += smearing->getSmearing(en1, iq1, ib1);
              deltaIso *= 0.5;
              break;
            }

            double termIso = 0.;
            for (int iat = 0; iat < numAtoms; iat++) {
              std::complex<double> zzIso = complexZero;
              for (int kdim : {0, 1, 2}) { // cartesian indices
                zzIso += std::conj(ev1(kdim, iat, ib1)) * ev2(kdim, iat, ib2);
              }
              termIso += std::norm(zzIso) * massVariance(iat);
            }
            termIso *= pi * 0.5 / innerNumFullPoints * en1 * en2 * deltaIso;

            for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
              double bose1 = outerBose(iCalc, 0, ind1);
              double bose2 = innerBose(iCalc, 0, ind2);

              double rateIso =
                  termIso * (bose1 * bose2 + 0.5 * (bose1 + bose2));

              switch (switchCase) {
              case (0):
                // case of matrix construction
                // we build the scattering matrix S
                theMatrix(ind1, ind2) += rateIso;
                linewidth->operator()(iCalc, 0, ind1) += rateIso;
                break;
              case (1):
                // case of matrix-vector multiplication
                // we build the scattering matrix A = S*n(n+1)
                for (unsigned int iInput = 0; iInput < inPopulations.size(); iInput++) {
                  for (int i = 0; i < dimensionality_; i++) {
                    outPopulations[iInput](iCalc, i, ind1) +=
                        rateIso *
                        inPopulations[iInput](iCalc, i, ind2);
                    outPopulations[iInput](iCalc, i, ind1) +=
                        rateIso *
                        inPopulations[iInput](iCalc, i, ind1);
                  }
                }
                break;
              case (2):
                // case of linewidth construction
                // there's still a missing norm done later
                linewidth->operator()(iCalc, 0, ind1) += rateIso;
                break;
              }
            }
          }
        }
      }
    }
  }

  if (switchCase == 1) {
    for (unsigned int i=0; i<outPopulations.size(); i++) {
      mpi->allReduceSum(&outPopulations[i].data);
    }
  } else {
    mpi->allReduceSum(&linewidth->data);
  }
  // I prefer to close loopPrint after the MPI barrier: all MPI are synced here
  loopPrint.close();

  // Add boundary scattering

  if (doBoundary) {
#pragma omp parallel for
    for (long is1 = 0; is1 < numStates; is1++) {
      double energy = outerBandStructure.getEnergy(is1);
      auto vel = outerBandStructure.getGroupVelocity(is1);
      for (long iCalc = 0; iCalc < statisticsSweep.getNumCalcs(); iCalc++) {
        double temperature =
            statisticsSweep.getCalcStatistics(iCalc).temperature;
        // n(n+1)
        double termPop = particle.getPopPopPm1(energy, temperature);
        double rate = vel.squaredNorm() / boundaryLength * termPop;

        switch (switchCase) {
        case (0):
          // case of matrix construction
          theMatrix(is1, is1) += rate; // this is redundant, see below
          linewidth->operator()(iCalc, 0, is1) += rate;
          break;
        case (1):
          // case of matrix-vector multiplication
          for (unsigned int iInput = 0; iInput < inPopulations.size(); iInput++) {
            for (int i = 0; i < dimensionality_; i++) {
              outPopulations[iInput](iCalc, i, is1) +=
                  rate * inPopulations[iInput](iCalc, i, is1);
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
        theMatrix(is1, is2) = 0.;
      }
    }

  } else if (switchCase == 1) {
    // case of matrix-vector multiplication
    for (auto is1 : excludeIndeces) {
      for (unsigned int i=0; i<outPopulations.size(); i++) {
        outPopulations[i].data.col(is1).setZero();
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
  if (switchCase == 0) { // case of matrix construction
    long iCalc = 0;
    for (long i = 0; i < outerBandStructure.getNumStates(); i++) {
      theMatrix(i, i) = linewidth->operator()(iCalc, 0, i);
    }
  }
}
