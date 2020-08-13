#include "helper_3rd_state.h"
#include <set>
#include "constants.h"
#include "mpiHelper.h"
#include "delta_function.h"

Helper3rdState::Helper3rdState(BaseBandStructure &innerBandStructure_,
                               BaseBandStructure &outerBandStructure_,
                               VectorBTE &outerBose_, const int &smearingType_,
                               PhononH0 *h0_)
    : innerBandStructure(innerBandStructure_),
      outerBandStructure(outerBandStructure_),
      outerBose(outerBose_),
      smearingType(smearingType_),
      h0(h0_) {
  // three conditions must be met to avoid recomputing q3
  // 1 - q1 and q2 mesh must be the same
  // 2 - the mesh is gamma-centered
  // 3 - the mesh is complete (if q1 and q2 are only around 0, q3 might be
  //     at the border)
  auto tup = outerBandStructure.getPoints().getMesh();
  auto mesh = std::get<0>(tup);
  auto offset = std::get<1>(tup);
  if ((&innerBandStructure == &outerBandStructure) && (offset.norm() == 0.) &&
      innerBandStructure.hasWindow() == 0) {
    storedAllQ3 = true;
    storedAllQ3Case = storedAllQ3Case1;

  } else if ((&innerBandStructure == &outerBandStructure) &&
             (offset.norm() == 0.) && innerBandStructure.hasWindow() != 0) {
    storedAllQ3 = true;
    storedAllQ3Case = storedAllQ3Case2;

    // in this case, we filtered some wavevectors out of inner and outer
    // bandstructures. We must find all q3 that conserve momentum

    // first, we build the full grid that the 3rd point would fall into
    auto innerPoints = innerBandStructure.getPoints();
    auto tup = innerPoints.getMesh();
    auto mesh = std::get<0>(tup);
    auto offset = std::get<1>(tup);
    fullPoints3 = std::make_unique<FullPoints>(
        innerBandStructure.getPoints().getCrystal(), mesh, offset);

    // now, we loop over the pairs of wavevectors
    std::set<int> listOfIndexes;
    int numPoints = innerPoints.getNumPoints();
    for (int iq2 : mpi->divideWorkIter(numPoints)) {
      for (int iq1 = 0; iq1 < numPoints; iq1++) {
        Eigen::Vector3d q1Coords =
            innerPoints.getPointCoords(iq1, Points::crystalCoords);
        Eigen::Vector3d q2Coords =
            innerPoints.getPointCoords(iq2, Points::crystalCoords);

        Eigen::Vector3d q3PlusCoords = q1Coords + q2Coords;
        Eigen::Vector3d q3MinsCoords = q1Coords - q2Coords;

        int iq3Plus = fullPoints3->getIndex(q3PlusCoords);
        int iq3Mins = fullPoints3->getIndex(q3MinsCoords);

        listOfIndexes.insert(iq3Plus);
        listOfIndexes.insert(iq3Mins);
      }
    }

    std::vector<int> localCounts(mpi->getSize(), 0);
    localCounts[mpi->getRank()] = listOfIndexes.size();
    mpi->allReduceSum(&localCounts);
    int totalCounts = 0;
    for ( int iRank=0; iRank<mpi->getSize(); iRank++ ) {
      totalCounts += localCounts[iRank];
    }

    // these are the offset to be used for MPI gather
    std::vector<int> displacements(mpi->getSize(),0);
    for (int i = 1; i < mpi->getSize(); i++) {
      displacements[i] = displacements[i - 1] + localCounts[i-1];
    }

    std::vector<int> globalListOfIndexes(totalCounts,0);
    int i=0;
    for (auto it : listOfIndexes) {
      int index = i + displacements[mpi->getRank()];
      i++; // after computing index!
      globalListOfIndexes[index] = it;
    }
    mpi->allReduceSum(&globalListOfIndexes);

    // now we must avoid duplicates between different MPI processes
    std::set<int> setOfIndexes;
    for ( auto x : globalListOfIndexes ) {
      setOfIndexes.insert(x);
    }

    // create the filtered list of points
    Eigen::VectorXi filter(setOfIndexes.size());
    i = 0;
    for (long iq : setOfIndexes) {
      filter(i) = iq;
      i++;
    }
    ActivePoints activePoints3(*fullPoints3, filter);

    // build band structure
    bool withEigenvectors = true;
    bool withVelocities = true;
    std::unique_ptr<ActiveBandStructure> bs(new ActiveBandStructure(
        activePoints3, h0, withEigenvectors, withVelocities));
    bandStructure3 = std::make_unique<ActiveBandStructure>(
        activePoints3, h0, withEigenvectors, withVelocities);
  }
}

// auto [eigvals3Mins, nb3Mins, eigvecs3Mins, v3Mins, bose3]
/** This function receives in input the cartesian coordinates of a vector,
 * and returns the harmonic info for that vector.
 * This is to be used for the third wavevector of the 3-phonon scattering.
 */
std::tuple<Eigen::Vector3d, Eigen::VectorXd, long, Eigen::MatrixXcd,
           Eigen::MatrixXd, Eigen::MatrixXd>
Helper3rdState ::get(Point &point1, Point &point2, const int &thisCase) {
  Eigen::Vector3d q3;
  if (thisCase == casePlus) {
    q3 = point1.getCoords(Points::cartesianCoords) +
         point2.getCoords(Points::cartesianCoords);
  } else {
    q3 = point1.getCoords(Points::cartesianCoords) -
         point2.getCoords(Points::cartesianCoords);
  }

  if (storedAllQ3) {
    // if the meshes are the same (and gamma centered)
    // q3 will fall into the same grid, and it's easy to get

    // note: 3rdBandStructure might still be different from inner/outer bs.
    // so, we must use the points from 3rdBandStructure to get the values

    long iq3;
    long nb3;
    Eigen::VectorXd energies3;
    Eigen::MatrixXcd eigvecs3;
    Eigen::MatrixXd v3s;
    Eigen::MatrixXd bose3Data;

    if (storedAllQ3Case == storedAllQ3Case1) {  // we use innerBandStruc
      auto points = innerBandStructure.getPoints();
      Eigen::Vector3d crystalPoints = points.cartesianToCrystal(q3);
      iq3 = points.getIndex(crystalPoints);
      auto iq3Index = WavevectorIndex(iq3);
      energies3 = innerBandStructure.getEnergies(iq3Index);
      eigvecs3 = innerBandStructure.getEigenvectors(iq3Index);

      if (smearingType == DeltaFunction::adaptiveGaussian) {
        v3s = innerBandStructure.getGroupVelocities(iq3Index);
      }
      nb3 = energies3.size();
      bose3Data = Eigen::MatrixXd(outerBose.numCalcs, nb3);
      for (long ib3 = 0; ib3 < nb3; ib3++) {
        long ind3 =
            outerBandStructure.getIndex(WavevectorIndex(iq3), BandIndex(ib3));
        bose3Data.col(ib3) = outerBose.data.col(ind3);
      }
    } else {
      auto points = bandStructure3->getPoints();
      Eigen::Vector3d crystalPoints = points.cartesianToCrystal(q3);
      iq3 = points.getIndex(crystalPoints);
      auto iq3Index = WavevectorIndex(iq3);
      energies3 = bandStructure3->getEnergies(iq3Index);
      eigvecs3 = bandStructure3->getEigenvectors(iq3Index);
      if (smearingType == DeltaFunction::adaptiveGaussian) {
        v3s = bandStructure3->getGroupVelocities(iq3Index);
      }
      nb3 = energies3.size();
      bose3Data = Eigen::MatrixXd(outerBose.numCalcs, nb3);
      for (long iCalc = 0; iCalc < outerBose.numCalcs; iCalc++) {
        double temperature =
            outerBose.statisticsSweep.getCalcStatistics(iCalc).temperature;
        for (long ib3 = 0; ib3 < nb3; ib3++) {
          bose3Data(iCalc, ib3) = h0->getParticle().getPopulation(
              energies3(ib3), temperature);
        }
      }
    }

    return {q3, energies3, nb3, eigvecs3, v3s, bose3Data};

  } else {
    // otherwise, q3 doesn't fall into the same grid
    // and we must therefore compute it from the hamiltonian

    long iq1 = point1.getIndex();

    Eigen::VectorXd energies3;
    Eigen::MatrixXcd eigvecs3;
    Eigen::MatrixXd v3s;
    Eigen::MatrixXd bose3Data;

    if (thisCase == casePlus) {
      energies3 = cachePlusEnergies[iq1];
      eigvecs3 = cachePlusEigvecs[iq1];
      v3s = cachePlusVelocity[iq1];
      bose3Data = cachePlusBose[iq1];
    } else {
      energies3 = cacheMinsEnergies[iq1];
      eigvecs3 = cacheMinsEigvecs[iq1];
      v3s = cacheMinsVelocity[iq1];
      bose3Data = cacheMinsBose[iq1];
    }

    long nb3 = energies3.size();
    return {q3, energies3, nb3, eigvecs3, v3s, bose3Data};
  }
}

void Helper3rdState::prepare(const std::vector<long> q1Indexes,
                             const long &iq2) {
  if (!storedAllQ3) {
    int numPoints = q1Indexes.size();

    cachePlusEnergies.resize(numPoints);
    cachePlusEigvecs.resize(numPoints);
    cachePlusBose.resize(numPoints);
    cachePlusVelocity.resize(numPoints);

    cacheMinsEnergies.resize(numPoints);
    cacheMinsEigvecs.resize(numPoints);
    cacheMinsBose.resize(numPoints);
    cacheMinsVelocity.resize(numPoints);

    Particle particle = h0->getParticle();

    for (long iq1 : q1Indexes) {
      Eigen::Vector3d q1 =
          outerBandStructure.getPoint(iq1).getCoords(Points::cartesianCoords);
      Eigen::Vector3d q2 =
          innerBandStructure.getPoint(iq2).getCoords(Points::cartesianCoords);

      Eigen::Vector3d q3Plus = q1 + q2;
      Eigen::Vector3d q3Mins = q1 - q2;

      auto tup = h0->diagonalizeFromCoords(q3Plus);
      auto energies3Plus = std::get<0>(tup);
      auto eigvecs3Plus = std::get<1>(tup);
      auto tup1 = h0->diagonalizeFromCoords(q3Mins);
      auto energies3Mins = std::get<0>(tup1);
      auto eigvecs3Mins = std::get<1>(tup1);

      int nb3Plus = energies3Plus.size();
      int nb3Mins = energies3Mins.size();

      Eigen::MatrixXd bose3DataPlus(outerBose.numCalcs, nb3Plus);
      Eigen::MatrixXd bose3DataMins(outerBose.numCalcs, nb3Mins);
      bose3DataPlus.setZero();
      bose3DataMins.setZero();

      for (long iCalc = 0; iCalc < outerBose.numCalcs; iCalc++) {
        double temperature =
            outerBose.statisticsSweep.getCalcStatistics(iCalc).temperature;
        for (long ib3 = 0; ib3 < nb3Plus; ib3++) {
          bose3DataPlus(iCalc, ib3) =
              particle.getPopulation(energies3Plus(ib3), temperature);
        }
        for (long ib3 = 0; ib3 < nb3Mins; ib3++) {
          bose3DataMins(iCalc, ib3) =
              particle.getPopulation(energies3Mins(ib3), temperature);
        }
      }

      Eigen::MatrixXd v3sPlus(nb3Plus, 3);
      Eigen::MatrixXd v3sMins(nb3Mins, 3);
      v3sPlus.setZero();
      v3sMins.setZero();
      if (smearingType == DeltaFunction::adaptiveGaussian) {
        Eigen::Tensor<std::complex<double>, 3> v3sTmpPlus =
            h0->diagonalizeVelocityFromCoords(q3Plus);
        Eigen::Tensor<std::complex<double>, 3> v3sTmpMins =
            h0->diagonalizeVelocityFromCoords(q3Mins);

        // we only need the diagonal elements of the velocity operator
        // i.e. the group velocity
        for (int i : {0, 1, 2}) {
          for (int ib3 = 0; ib3 < nb3Plus; ib3++) {
            v3sPlus(ib3, i) = v3sTmpPlus(ib3, ib3, i).real();
          }
          for (int ib3 = 0; ib3 < nb3Mins; ib3++) {
            v3sMins(ib3, i) = v3sTmpMins(ib3, ib3, i).real();
          }
        }
      }

      cachePlusEnergies[iq1] = energies3Plus;
      cachePlusEigvecs[iq1] = eigvecs3Plus;
      cachePlusBose[iq1] = bose3DataPlus;
      cachePlusVelocity[iq1] = v3sPlus;

      cacheMinsEnergies[iq1] = energies3Mins;
      cacheMinsEigvecs[iq1] = eigvecs3Mins;
      cacheMinsBose[iq1] = bose3DataMins;
      cacheMinsVelocity[iq1] = v3sMins;
    }
  }
}

const int Helper3rdState::casePlus = 0;

const int Helper3rdState::caseMins = 1;
