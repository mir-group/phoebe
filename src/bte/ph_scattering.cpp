#include "ph_scattering.h"

#include "constants.h"
#include "io.h"
#include "mpiHelper.h"
#include "periodic_table.h"

PhScatteringMatrix::PhScatteringMatrix(Context &context_,
                                       StatisticsSweep &statisticsSweep_,
                                       BaseBandStructure &innerBandStructure_,
                                       BaseBandStructure &outerBandStructure_,
                                       Interaction3Ph *coupling3Ph_,
                                       PhononH0 *h0_)
    : ScatteringMatrix(context_, statisticsSweep_, innerBandStructure_,
                       outerBandStructure_),
      coupling3Ph(coupling3Ph_), h0(h0_) {
  //	couplingIsotope = couplingIsotope_;
  //	couplingBoundary = couplingBoundary_;
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
void PhScatteringMatrix::builder(ParallelMatrix<double> &matrix,
                                 VectorBTE *linewidth, VectorBTE *inPopulation,
                                 VectorBTE *outPopulation) {
  // notes: + process is (1+2) -> 3
  //        - processes are (1+3)->2 and (3+2)->1

  const double energyCutoff = 0.001 / ryToCmm1; // discard states with small
  // energies (smaller than 0.001 cm^-1

  int switchCase = 0;
  if (matrix.rows() != 0 && linewidth != nullptr && inPopulation == nullptr &&
      outPopulation == nullptr) {
    switchCase = 0;
  } else if (matrix.rows() == 0 && linewidth == nullptr &&
             inPopulation != nullptr && outPopulation != nullptr) {
    switchCase = 1;
  } else if (matrix.rows() == 0 && linewidth != nullptr &&
             inPopulation == nullptr && outPopulation == nullptr) {
    switchCase = 2;
  } else {
    Error e("builder3Ph found a non-supported case");
  }

  if ((linewidth != nullptr) && (linewidth->dimensionality != 1)) {
    Error e("The linewidths shouldn't have dimensionality");
  }

  // three conditions must be met to avoid recomputing q3
  // 1 - q1 and q2 mesh must be the same
  // 2 - the mesh is gamma-centered
  // 3 - the mesh is complete (if q1 and q2 are only around 0, q3 might be
  //     at the border)
  auto tup = outerBandStructure.getPoints().getMesh();
  auto mesh = std::get<0>(tup);
  auto offset = std::get<1>(tup);
  bool dontComputeQ3;
  if ((&innerBandStructure == &outerBandStructure) && (offset.norm() == 0.) &&
      innerBandStructure.hasWindow() == 0) {
    dontComputeQ3 = true;
  } else {
    dontComputeQ3 = false;
  }

  auto particle = outerBandStructure.getParticle();

  long numAtoms = innerBandStructure.getPoints().getCrystal().getNumAtoms();
  long numCalcs = statisticsSweep.getNumCalcs();

  // note: innerNumFullPoints is the number of points in the full grid
  // may be larger than innerNumPoints, when we use ActiveBandStructure
  long innerNumFullPoints = innerBandStructure.getNumPoints(true);

  // precompute Bose populations
  VectorBTE outerBose(statisticsSweep, outerBandStructure, 1);
  for (long is = 0; is < outerBandStructure.getNumStates(); is++) {
    double energy = outerBandStructure.getEnergy(is);
    for (long iCalc = 0; iCalc < statisticsSweep.getNumCalcs(); iCalc++) {
      double temperature = statisticsSweep.getCalcStatistics(iCalc).temperature;
      outerBose.data(iCalc, is) = particle.getPopulation(energy, temperature);
    }
  }
  VectorBTE innerBose(statisticsSweep, outerBandStructure, 1);
  if (&innerBandStructure == &outerBandStructure) {
    innerBose = outerBose;
  } else {
    for (long is = 0; is < innerBandStructure.getNumStates(); is++) {
      double energy = innerBandStructure.getEnergy(is);
      for (long iCalc = 0; iCalc < statisticsSweep.getNumCalcs(); iCalc++) {
        double temperature =
            statisticsSweep.getCalcStatistics(iCalc).temperature;
        outerBose.data(iCalc, is) = particle.getPopulation(energy, temperature);
      }
    }
  }

  // note: these variables are only needed in the loop
  // but since it's an expensive loop, we define them here once and for all
  long nb3Plus, nb3Mins, ind1, ind2;
  double ratePlus, rateMins1, rateMins2, deltaPlus, deltaMins1, deltaMins2;
  double en1, en2, en3Plus, en3Mins, bose1, bose2, bose3Plus, bose3Mins;
  Eigen::Tensor<double, 3> couplingPlus, couplingMins;
  Eigen::VectorXd state3PlusEnergies, state3MinsEnergies;
  Eigen::Vector3d v1, v2, v3, v;
  Eigen::MatrixXd v1s, v2s, v3ps, v3ms;

  Eigen::MatrixXd bose3PlusData, bose3MinsData;
  Eigen::VectorXd eigvals3Plus, eigvals3Mins;
  Eigen::MatrixXcd eigvecs3Plus, eigvecs3Mins;

  // isotopic scattering:
  std::complex<double> zzIso;
  double termIso, rateIso, deltaIso;

  std::vector<std::tuple<long, long>> qPairIterator =
      getIteratorWavevectorPairs(switchCase);

  LoopPrint loopPrint("computing scattering matrix", "q-points",
                      qPairIterator.size());
  // invert qPairIterator into something workable
  int nq2 = 0;
  int previq2 = -1;
  std::vector<long> iq2s;
  std::vector<std::vector<long>> iq1ss;

  for (auto qtup : qPairIterator) {
    auto iq1 = std::get<0>(qtup);
    auto iq2 = std::get<1>(qtup);
    if (iq2 != previq2) {
      nq2 += 1;
      iq2s.push_back(iq2);
      iq1ss.push_back(std::vector<long>());
      previq2 = iq2;
    }
    iq1ss[nq2 - 1].push_back(iq1);
  }
  /*
  for (int i = 0; i < nq2; i++) {
    std::cout << iq2s[i] << " x";
    for (long iq1 : iq1s[i]) {
      std::cout << " " << iq1;
    }
    std::cout << "\n";
  }
  */

  for (int iiq2 = 0; iiq2 < nq2; iiq2++) {
    long iq2 = iq2s[iiq2];

    std::vector<long> iq1s = iq1ss[iiq2];
    int nq1 = iq1s.size();

    std::vector<Point> q1s, q3Pluss, q3Minss;
    std::vector<Eigen::VectorXd> state1Energiess, state3PlusEnergiess,
        state3MinsEnergiess;
    std::vector<Eigen::Vector3d> q1coords, q3PlusCs, q3MinsCs;
    std::vector<Eigen::Tensor<std::complex<double>, 3>> ev1s3d;
    std::vector<Eigen::MatrixXcd> ev1s, ev3ps, ev3ms;
    std::vector<Eigen::MatrixXd> v1ss, v3pss, v3mss;
    std::vector<long> nb1s, nb3Pluss, nb3Minss;
    std::vector<Eigen::Tensor<double, 3>> couplingPluss, couplingMinss;

    Point q2 = innerBandStructure.getPoint(iq2);
    State states2 = innerBandStructure.getState(q2);
    Eigen::VectorXd state2Energies = states2.getEnergies();
    int nb2 = state2Energies.size();

    Eigen::Tensor<std::complex<double>, 3> ev2;
    states2.getEigenvectors(ev2);

    // note: for computing linewidths on a path, we must distinguish
    // that q1 and q2 are on different meshes, and that q3+/- may not
    // fall into known meshes and therefore needs to be computed

    for (int iq1 : iq1s) {
      State states1 = outerBandStructure.getState(iq1);
      Point q1 = states1.getPoint();
      Eigen::VectorXd state1Energies = states1.getEnergies();
      int nb1 = state1Energies.size();

      Eigen::Tensor<std::complex<double>, 3> ev13d;
      Eigen::MatrixXcd ev1;
      states1.getEigenvectors(ev1);
      states1.getEigenvectors(ev13d);
      ev1s3d.push_back(ev13d);
      ev1s.push_back(ev1);

      v1s = states1.getGroupVelocities();
      v2s = states2.getGroupVelocities();

      q1s.push_back(q1);
      v1ss.push_back(v1s);
      state1Energiess.push_back(state1Energies);
      nb1s.push_back(nb1);
      q1coords.push_back(states1.getCoords(Points::cartesianCoords));

      // if the meshes are the same (and gamma centered)
      // q3 will fall into the same grid, and it's easy to get
      if (dontComputeQ3) {
        Point q3Plus = q1 + q2;
        State states3Plus = innerBandStructure.getState(q3Plus);
        state3PlusEnergies = states3Plus.getEnergies();
        nb3Plus = state3PlusEnergies.size();

        Point q3Mins = q1 - q2;
        State states3Mins = innerBandStructure.getState(q3Mins);
        state3MinsEnergies = states3Mins.getEnergies();
        nb3Mins = state3MinsEnergies.size();

        q3Pluss.push_back(q3Plus);
        q3Minss.push_back(q3Mins);

        q3PlusCs.push_back(Eigen::Vector3d(0, 0, 0));
        q3MinsCs.push_back(Eigen::Vector3d(0, 0, 0));

        nb3Pluss.push_back(nb3Plus);
        nb3Minss.push_back(nb3Mins);
        state3PlusEnergiess.push_back(state3PlusEnergies);
        state3MinsEnergiess.push_back(state3MinsEnergies);

        Eigen::MatrixXcd ev3p, ev3m;
        states3Plus.getEigenvectors(ev3p);
        states3Mins.getEigenvectors(ev3m);
        ev3ps.push_back(ev3p);
        ev3ms.push_back(ev3m);

        auto tup3 = coupling3Ph->getCouplingSquared(states1, states2,
                                                    states3Plus, states3Mins);
        auto cp = std::get<0>(tup3);
        auto cm = std::get<1>(tup3);
        couplingPlus = cp;
        couplingMins = cm;
        couplingPluss.push_back(couplingPlus);
        couplingMinss.push_back(couplingMins);

        v3ps = states3Plus.getGroupVelocities();
        v3ms = states3Mins.getGroupVelocities();

        v3pss.push_back(v3ps);
        v3mss.push_back(v3ms);

      } else {
        // otherwise, q3 doesn't fall into the same grid
        // and we must therefore compute it from the hamiltonian

        Eigen::Vector3d q3PlusC = q1.getCoords(Points::cartesianCoords) +
                                  q2.getCoords(Points::cartesianCoords);
        Eigen::Vector3d q3MinsC = q1.getCoords(Points::cartesianCoords) -
                                  q2.getCoords(Points::cartesianCoords);

        auto tup = h0->diagonalizeFromCoords(q3PlusC);
        auto eigvals3Plus = std::get<0>(tup);
        auto eigvecs3Plus = std::get<1>(tup);
        auto tup1 = h0->diagonalizeFromCoords(q3MinsC);
        auto eigvals3Mins = std::get<0>(tup1);
        auto eigvecs3Mins = std::get<1>(tup1);

        nb3Plus = eigvals3Plus.size();
        nb3Mins = eigvals3Mins.size();

        q3PlusCs.push_back(q3PlusC);
        q3MinsCs.push_back(q3MinsC);

        // pray that these are never used
        q3Pluss.push_back(q1);
        q3Minss.push_back(q1);

        DetachedState states3Plus(q3PlusC, eigvals3Plus, numAtoms, nb3Plus,
                                  eigvecs3Plus);
        DetachedState states3Mins(q3MinsC, eigvals3Mins, numAtoms, nb3Mins,
                                  eigvecs3Mins);

        state3PlusEnergies = states3Plus.getEnergies();
        state3MinsEnergies = states3Mins.getEnergies();
        nb3Pluss.push_back(nb3Plus);
        nb3Minss.push_back(nb3Mins);
        state3PlusEnergiess.push_back(state3PlusEnergies);
        state3MinsEnergiess.push_back(state3MinsEnergies);

        Eigen::MatrixXcd ev3p, ev3m;
        states3Plus.getEigenvectors(ev3p);
        states3Mins.getEigenvectors(ev3m);
        ev3ps.push_back(ev3p);
        ev3ms.push_back(ev3m);

        auto tup2 = coupling3Ph->getCouplingSquared(states1, states2,
                                                    states3Plus, states3Mins);
        auto cp = std::get<0>(tup2);
        auto cm = std::get<1>(tup2);
        couplingPlus = cp;
        couplingMins = cm;

        couplingPluss.push_back(couplingPlus);
        couplingMinss.push_back(couplingMins);

        v3pss.push_back(Eigen::MatrixXd(0, 0));
        v3mss.push_back(Eigen::MatrixXd(0, 0));
      }
    }

    for (int iiq1 = 0; iiq1 < nq1; iiq1++) {
      int iq1 = iq1s[iiq1];
      auto state1Energies = state1Energiess[iiq1];
      auto state3PlusEnergies = state3PlusEnergiess[iiq1];
      auto state3MinsEnergies = state3MinsEnergiess[iiq1];
      v1s = v1ss[iiq1];
      int nb3Plus = nb3Pluss[iiq1];
      int nb3Mins = nb3Minss[iiq1];
      auto couplingPlus = couplingPluss[iiq1];
      auto couplingMins = couplingMinss[iiq1];
      int nb1 = nb1s[iiq1];
      auto ev1 = ev1s3d[iiq1];

      if (dontComputeQ3) {
        auto q3Plus = q3Pluss[iiq1];

        auto q3Mins = q3Minss[iiq1];
        bose3PlusData = Eigen::MatrixXd::Zero(numCalcs, nb3Plus);
        bose3MinsData = Eigen::MatrixXd::Zero(numCalcs, nb3Plus);

        for (long ib3 = 0; ib3 < nb3Plus; ib3++) {
          long ind3 = outerBandStructure.getIndex(
              WavevectorIndex(q3Plus.getIndex()), BandIndex(ib3));
          bose3PlusData.col(ib3) = outerBose.data.col(ind3);
        }
        for (long ib3 = 0; ib3 < nb3Mins; ib3++) {
          long ind3 = outerBandStructure.getIndex(
              WavevectorIndex(q3Mins.getIndex()), BandIndex(ib3));
          bose3MinsData.col(ib3) = outerBose.data.col(ind3);
        }

        if (smearing->getType() == DeltaFunction::gaussian) {
          // v3ps = states3Plus.getGroupVelocities();
          // v3ms = states3Mins.getGroupVelocities();
          v3ps = v3pss[iiq1];
          v3ms = v3mss[iiq1];
        }

      } else {
        bose3PlusData = Eigen::MatrixXd::Zero(numCalcs, nb3Plus);
        bose3MinsData = Eigen::MatrixXd::Zero(numCalcs, nb3Plus);

        auto q3PlusC = q3PlusCs[iiq1];
        auto q3MinsC = q3MinsCs[iiq1];

        for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
          double temperature =
              statisticsSweep.getCalcStatistics(iCalc).temperature;
          for (long ib3 = 0; ib3 < nb3Plus; ib3++) {
            bose3PlusData(iCalc, ib3) =
                particle.getPopulation(state3PlusEnergies(ib3), temperature);
          }
          for (long ib3 = 0; ib3 < nb3Mins; ib3++) {
            bose3MinsData(iCalc, ib3) =
                particle.getPopulation(state3MinsEnergies(ib3), temperature);
          }
        }

        if (smearing->getType() == DeltaFunction::adaptiveGaussian) {
          Eigen::Tensor<std::complex<double>, 3> v3psTmp =
              h0->diagonalizeVelocityFromCoords(q3PlusC);
          Eigen::Tensor<std::complex<double>, 3> v3msTmp =
              h0->diagonalizeVelocityFromCoords(q3MinsC);

          // we only need the diagonal elements of the velocity operator
          // i.e. the group velocity
          v3ps = Eigen::MatrixXd::Zero(nb3Plus, 3);
          v3ms = Eigen::MatrixXd::Zero(nb3Mins, 3);
          for (int i : {0, 1, 2}) {
            for (int ib3 = 0; ib3 < nb3Plus; ib3++) {
              v3ps(ib3, i) = v3psTmp(ib3, ib3, i).real();
            }
            for (int ib3 = 0; ib3 < nb3Mins; ib3++) {
              v3ms(ib3, i) = v3msTmp(ib3, ib3, i).real();
            }
          }
        }
      }

      for (long ib1 = 0; ib1 < nb1; ib1++) {
        en1 = state1Energies(ib1);
        ind1 =
            outerBandStructure.getIndex(WavevectorIndex(iq1), BandIndex(ib1));
        if (en1 < energyCutoff) {
          continue;
        }

        for (long ib2 = 0; ib2 < nb2; ib2++) {
          en2 = state2Energies(ib2);
          ind2 =
              innerBandStructure.getIndex(WavevectorIndex(iq2), BandIndex(ib2));
          if (en2 < energyCutoff) {
            continue;
          }

          if (switchCase == 0) {
            // note: above we are parallelizing over wavevectors.
            // (for convenience of computing coupling3ph)
            // Not the same way as Matrix() is parallelized.
            // here we check that we don't duplicate efforts
            if (!matrix.indecesAreLocal(ind1, ind2)) {
              continue;
            }
          }

          // Isotope scattering
          if (doIsotopes) {
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

            termIso = 0.;
            for (int iat = 0; iat < numAtoms; iat++) {
              zzIso = complexZero;
              for (int kdim : {0, 1, 2}) { // cartesian indices
                zzIso += std::conj(ev1(kdim, iat, ib1)) * ev2(kdim, iat, ib2);
              }
              termIso += std::norm(zzIso) * massVariance(iat);
            }
            termIso *= pi * 0.5 / innerNumFullPoints * en1 * en2 * deltaIso;

            for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
              bose1 = outerBose.data(iCalc, ind1);
              bose2 = innerBose.data(iCalc, ind2);

              rateIso = termIso * (bose1 * bose2 + 0.5 * (bose1 + bose2));

              switch (switchCase) {
              case (0):
                // case of matrix construction
                // we build the scattering matrix S
                matrix(ind1, ind2) += rateIso;
                linewidth->data(iCalc, ind1) += rateIso;
                break;
              case (1):
                // case of matrix-vector multiplication
                // we build the scattering matrix A = S*n(n+1)
                for (long i : {0, 1, 2}) {
                  outPopulation->data(3 * iCalc + i, ind1) +=
                      rateIso * inPopulation->data(3 * iCalc + i, ind2);
                  outPopulation->data(3 * iCalc + i, ind1) +=
                      rateIso * inPopulation->data(3 * iCalc + i, ind1);
                }
                break;
              case (2):
                // case of linewidth construction
                // there's still a missing norm done later
                linewidth->data(iCalc, ind1) += rateIso;
                break;
              }
            }
          }

          // split into two cases since there may be different bands
          for (long ib3 = 0; ib3 < nb3Plus; ib3++) {
            en3Plus = state3PlusEnergies(ib3);
            if (en3Plus < energyCutoff) {
              continue;
            }

            double enProd = en1 * en2 * en3Plus;

            switch (smearing->getType()) {
            case (DeltaFunction::gaussian):
              deltaPlus = smearing->getSmearing(en1 + en2 - en3Plus);
              break;
            case (DeltaFunction::adaptiveGaussian):
              v = v2s.row(ib2) - v3ps.row(ib3);
              deltaPlus = smearing->getSmearing(en1 + en2 - en3Plus, v);
              break;
            default:
              deltaPlus = smearing->getSmearing(en3Plus - en1, iq2, ib2);
              break;
            }

            if (deltaPlus < 0)
              continue;

            // loop on temperature
            for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
              bose1 = outerBose.data(iCalc, ind1);
              bose2 = innerBose.data(iCalc, ind2);
              bose3Plus = bose3PlusData(iCalc, ib3);

              // Calculate transition probability W+
              ratePlus = pi * 0.25 * bose1 * bose2 * (bose3Plus + 1.) *
                         couplingPlus(ib1, ib2, ib3) * deltaPlus /
                         innerNumFullPoints / enProd;

              switch (switchCase) {
              case (0):
                // case of matrix construction
                // we build the scattering matrix S
                matrix(ind1, ind2) += ratePlus;
                linewidth->data(iCalc, ind1) += ratePlus;
                break;
              case (1):
                // case of matrix-vector multiplication
                // we build the scattering matrix A = S*n(n+1)
                for (long i : {0, 1, 2}) {
                  outPopulation->data(3 * iCalc + i, ind1) +=
                      ratePlus * inPopulation->data(3 * iCalc + i, ind2);
                  outPopulation->data(3 * iCalc + i, ind1) +=
                      ratePlus * inPopulation->data(3 * iCalc + i, ind1);
                }
                break;
              case (2):
                // case of linewidth construction
                linewidth->data(iCalc, ind1) += ratePlus;
                break;
              }
            }
          }

          for (long ib3 = 0; ib3 < nb3Mins; ib3++) {
            en3Mins = state3MinsEnergies(ib3);
            if (en3Mins < energyCutoff) {
              continue;
            }

            double enProd = en1 * en2 * en3Mins;

            switch (smearing->getType()) {
            case (DeltaFunction::gaussian):
              deltaMins1 = smearing->getSmearing(en1 + en3Mins - en2);
              deltaMins2 = smearing->getSmearing(en2 + en3Mins - en1);
              break;
            case (DeltaFunction::adaptiveGaussian):
              v = v2s.row(ib2) - v3ms.row(ib3);
              deltaMins1 = smearing->getSmearing(en1 + en3Mins - en2, v);
              deltaMins2 = smearing->getSmearing(en2 + en3Mins - en1, v);
              break;
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

            for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
              bose1 = outerBose.data(iCalc, ind1);
              bose2 = innerBose.data(iCalc, ind2);
              bose3Mins = bose3MinsData(iCalc, ib3);

              // Calculatate transition probability W-
              rateMins1 = pi * 0.25 * bose3Mins * bose1 * (bose2 + 1.) *
                          couplingMins(ib1, ib2, ib3) * deltaMins1 /
                          innerNumFullPoints / enProd;
              rateMins2 = pi * 0.25 * bose2 * bose3Mins * (bose1 + 1.) *
                          couplingMins(ib1, ib2, ib3) * deltaMins2 /
                          innerNumFullPoints / enProd;

              switch (switchCase) {
              case (0):
                // case of matrix construction
                matrix(ind1, ind2) -= rateMins1 + rateMins2;
                linewidth->data(iCalc, ind1) += 0.5 * rateMins2;
                break;
              case (1):
                // case of matrix-vector multiplication
                for (long i : {0, 1, 2}) {
                  outPopulation->data(3 * iCalc + i, ind1) -=
                      (rateMins1 + rateMins2) *
                      inPopulation->data(3 * iCalc + i, ind2);
                  outPopulation->data(3 * iCalc + i, ind1) +=
                      0.5 * rateMins2 * inPopulation->data(3 * iCalc + i, ind1);
                }
                break;
              case (2):
                // case of linewidth construction
                linewidth->data(iCalc, ind1) += 0.5 * rateMins2;
                break;
              }
            }
          }
        }
      }
    }
  }
  loopPrint.close();

  if (switchCase == 1) {
    mpi->allReduceSum(&outPopulation->data);
  } else {
    mpi->allReduceSum(&linewidth->data);
  }

  // Add boundary scattering

  if (doBoundary) {
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
          matrix(is1, is1) += rate; // this is redundant, see below
          linewidth->data(iCalc, is1) += rate;
          break;
        case (1):
          // case of matrix-vector multiplication
          for (long i : {0, 1, 2}) {
            outPopulation->data(3 * iCalc + i, is1) +=
                rate * inPopulation->data(3 * iCalc + i, is1);
          }
          break;
        case (2):
          // case of linewidth construction
          linewidth->data(iCalc, is1) += rate;
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
        matrix(is1, is2) = 0.;
      }
    }

  } else if (switchCase == 1) {
    // case of matrix-vector multiplication
    for (auto is1 : excludeIndeces) {
      outPopulation->data.col(is1).setZero();
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
      matrix(i, i) = linewidth->data(iCalc, i);
    }
  }
}
