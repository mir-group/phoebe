#include "helper_el_scattering.h"
#include <set>
#include "constants.h"
#include "mpiHelper.h"
#include "delta_function.h"

HelperElScattering::HelperElScattering(BaseBandStructure &innerBandStructure_,
                               BaseBandStructure &outerBandStructure_,
                               StatisticsSweep &statisticsSweep_,
                               const int &smearingType_, PhononH0 &h0_)
    : innerBandStructure(innerBandStructure_),
      outerBandStructure(outerBandStructure_),
      statisticsSweep(statisticsSweep_),
      smearingType(smearingType_),
      h0(h0_) {
  // three conditions must be met to avoid recomputing q3
  // 1 - q1 and q2 mesh must be the same
  // 2 - the mesh is gamma-centered
  // 3 - the mesh is complete (if q1 and q2 are only around 0, q3 might be
  //     at the border)
  auto t1 = outerBandStructure.getPoints().getMesh();
//  auto mesh = std::get<0>(t1);
  auto offset = std::get<1>(t1);
  storedAllQ3 = false;

  if ((&innerBandStructure == &outerBandStructure) && (offset.norm() == 0.) &&
      innerBandStructure.hasWindow() == 0) {

    storedAllQ3 = true;
    storedAllQ3Case = storedAllQ3Case1;

    auto t2 = innerBandStructure.getPoints().getMesh();
    auto mesh = std::get<0>(t2);
    auto offset = std::get<1>(t2);

    fullPoints3 = std::make_unique<FullPoints>(
        innerBandStructure.getPoints().getCrystal(), mesh, offset);
    bool withVelocities = true;
    bool withEigenvectors = true;
    FullBandStructure bs = h0.populate(*fullPoints3, withVelocities,
                                               withEigenvectors);
    bandStructure3 = std::make_unique<FullBandStructure>(bs);

  } else if ((&innerBandStructure == &outerBandStructure) &&
             (offset.norm() == 0.) && innerBandStructure.hasWindow() != 0) {

    storedAllQ3 = true;
    storedAllQ3Case = storedAllQ3Case2;

    // in this case, we filtered some wavevectors out of inner and outer
    // bandstructures. We must find all q3 that conserve momentum

    // first, we build the full grid that the 3rd point would fall into
    auto innerPoints = innerBandStructure.getPoints();
    auto t2 = innerPoints.getMesh();
    auto mesh = std::get<0>(t2);
    auto offset = std::get<1>(t2);

    fullPoints3 = std::make_unique<FullPoints>(
        innerBandStructure.getPoints().getCrystal(), mesh, offset);

    // now, we loop over the pairs of wavevectors
    std::set<int> listOfIndexes;
    int numPoints = innerBandStructure.getNumPoints();
    for (int ik2 : mpi->divideWorkIter(numPoints)) {
      auto ik2Index = WavevectorIndex(ik2);
      Eigen::Vector3d k2Coords_ = innerBandStructure.getWavevector(ik2Index);

      auto rotations = innerPoints.getRotationsStar(ik2);
      for ( Eigen::Matrix3d rotation : rotations ) {
        Eigen::Vector3d k2Coords = rotation * k2Coords_;
        // Eigen::Vector3d k2Coords = k2Coords_;

        for (int ik1 = 0; ik1 < outerBandStructure.getNumPoints(); ik1++) {
          auto ik1Index = WavevectorIndex(ik1);
          Eigen::Vector3d k1Coords = outerBandStructure.getWavevector(ik1Index);

          // k' = k + q : phonon absorption
          Eigen::Vector3d q3Coords = k2Coords - k1Coords;
          Eigen::Vector3d q3Cart = fullPoints3->cartesianToCrystal(q3Coords);

          int iq3 = fullPoints3->getIndex(q3Cart);
          listOfIndexes.insert(iq3);
        }
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
    for (int iq : setOfIndexes) {
      filter(i) = iq;
      i++;
    }

    ActivePoints ap3 = ActivePoints(*fullPoints3, filter);
    activePoints3 = std::make_unique<ActivePoints>(ap3);

    // build band structure
    bool withEigenvectors = true;
    bool withVelocities = true;
    // note: bandStructure3 stores a copy of ap3: can't pass unique_ptr
    bandStructure3 = std::make_unique<ActiveBandStructure>(
        ap3, &h0, withEigenvectors, withVelocities);
  }
}

// auto [eigvals3Mins, nb3Mins, eigvecs3Mins, v3Mins, bose3]
/** This function receives in input the cartesian coordinates of a vector,
 * and returns the harmonic info for that vector.
 * This is to be used for the third wavevector of the 3-phonon scattering.
 */
std::tuple<Eigen::Vector3d, Eigen::VectorXd, int, Eigen::MatrixXcd,
           Eigen::MatrixXd, Eigen::MatrixXd>
    HelperElScattering::get(Eigen::Vector3d &k1, const int &ik2) {

  auto ik2Idx = WavevectorIndex(ik2);
  Eigen::Vector3d k2 = innerBandStructure.getWavevector(ik2Idx);
  Eigen::Vector3d q3 = k2 - k1;

  if (storedAllQ3) {
    // if the meshes are the same (and gamma centered)
    // q3 will fall into the same grid, and it's easy to get

    // note: 3rdBandStructure might still be different from inner/outer bs.
    // so, we must use the points from 3rdBandStructure to get the values

    int iq3;
    if (storedAllQ3Case == storedAllQ3Case1) {  // we use innerBandStruc
      Eigen::Vector3d crystalPoints = fullPoints3->cartesianToCrystal(q3);
      iq3 = fullPoints3->getIndex(crystalPoints);
    } else {
      Eigen::Vector3d crystalPoints = activePoints3->cartesianToCrystal(q3);
      iq3 = activePoints3->getIndex(crystalPoints);
    }
    auto iq3Index = WavevectorIndex(iq3);

    Eigen::VectorXd energies3 = bandStructure3->getEnergies(iq3Index);
    Eigen::MatrixXcd eigvecs3 = bandStructure3->getEigenvectors(iq3Index);
    Eigen::MatrixXd v3s;
    if (smearingType == DeltaFunction::adaptiveGaussian) {
      v3s = bandStructure3->getGroupVelocities(iq3Index);
    }
    int nb3 = energies3.size();
    auto particle = h0.getParticle();
    Eigen::MatrixXd bose3Data = Eigen::MatrixXd(statisticsSweep.getNumCalcs(), nb3);
    for (int iCalc = 0; iCalc < statisticsSweep.getNumCalcs(); iCalc++) {
      double temp = statisticsSweep.getCalcStatistics(iCalc).temperature;
      for (int ib3 = 0; ib3 < nb3; ib3++) {
         bose3Data(iCalc, ib3) = particle.getPopulation(energies3(ib3), temp);
      }
    }

    return {q3, energies3, nb3, eigvecs3, v3s, bose3Data};

  } else {
    // otherwise, q3 doesn't fall into the same grid
    // and we must therefore compute it from the hamiltonian

    // iq2 in input is an index over wavevectors
    // we need to find the index over the local cache
    int ik2Counter = ik2 - cacheOffset;
    Eigen::VectorXd energies3 = cacheEnergies[ik2Counter];
    Eigen::MatrixXcd eigvecs3 = cacheEigvecs[ik2Counter];
    Eigen::MatrixXd v3s = cacheVelocity[ik2Counter];
    Eigen::MatrixXd bose3Data = cacheBose[ik2Counter];
    int nb3 = energies3.size();

    return {q3, energies3, nb3, eigvecs3, v3s, bose3Data};
  }
}

void HelperElScattering::prepare(const Eigen::Vector3d &k1,
                                 const std::vector<int> k2Indexes) {
  if (!storedAllQ3) {
    int numPoints = k2Indexes.size();
    cacheEnergies.resize(numPoints);
    cacheEigvecs.resize(numPoints);
    cacheBose.resize(numPoints);
    cacheVelocity.resize(numPoints);
    cacheOffset = k2Indexes[0];

    Particle particle = h0.getParticle();

    int ik2Counter = -1;
    for (int ik2 : k2Indexes) {
      ik2Counter++;
      auto ik2Idx = WavevectorIndex(ik2);
      Eigen::Vector3d k2 = innerBandStructure.getWavevector(ik2Idx);

      Eigen::Vector3d q3 = k2 - k1;

      auto t1 = h0.diagonalizeFromCoords(q3);
      auto energies3 = std::get<0>(t1);
      auto eigvecs3 = std::get<1>(t1);

      int nb3 = energies3.size();

      Eigen::MatrixXd bose3Data(statisticsSweep.getNumCalcs(), nb3);
      bose3Data.setZero();

      for (int iCalc = 0; iCalc < statisticsSweep.getNumCalcs(); iCalc++) {
        double temp = statisticsSweep.getCalcStatistics(iCalc).temperature;
        for (int ib3 = 0; ib3 < nb3; ib3++) {
          bose3Data(iCalc, ib3) = particle.getPopulation(energies3(ib3), temp);
        }
      }

      Eigen::MatrixXd v3s(nb3, 3);
      v3s.setZero();
      if (smearingType == DeltaFunction::adaptiveGaussian) {
        Eigen::Tensor<std::complex<double>, 3> v3sTmp =
            h0.diagonalizeVelocityFromCoords(q3);

        // we only need the diagonal elements of the velocity operator
        // i.e. the group velocity
        for (int i : {0, 1, 2}) {
          for (int ib3 = 0; ib3 < nb3; ib3++) {
            v3s(ib3, i) = v3sTmp(ib3, ib3, i).real();
          }
        }
      }

      cacheEnergies[ik2Counter] = energies3;
      cacheEigvecs[ik2Counter] = eigvecs3;
      cacheBose[ik2Counter] = bose3Data;
      cacheVelocity[ik2Counter] = v3s;

    }
  }
}
