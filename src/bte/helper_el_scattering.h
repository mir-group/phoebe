#ifndef HELPER_EL_SCATTERING_H
#define HELPER_EL_SCATTERING_H

#include "phonon_h0.h"
#include "vector_bte.h"

/** This class manages the calculation of phonon properties during the
 * intragration of the scattering matrix.
 *
 * If possible, we precompute all phonon properties. This happens if inner=outer
 * bandstructure, see ElScatteringMatrix for explanation of what they are.
 * Otherwise, we compute phonon properties in batches (all q=k'-k for a given k)
 */
class HelperElScattering {
 public:

  /** Default constructor
   * In the constructor, if possible (outer==inner bandstructure), we precompute
   * all phonon properties for all valid q=k'-k points
   *
   * @param innerBandStructure_: bandstructure for scattering integration
   * @param outerBandStructure_: bandstructure for defining where we want
   * scattering properties computed.
   * @param statisticsSweep_: object with temperatures and chemical potentials
   * @param smearingType_: kind of smearing used for approximating the dirac
   * delta on energies in the scattering matrix. If the smearing is adaptive, we
   * also need to compute phonon group velocities.
   * @param h0_: Phonon hamiltonian needed to compute the phonon properties.
   */
  HelperElScattering(BaseBandStructure &innerBandStructure_,
                 BaseBandStructure &outerBandStructure_,
                 StatisticsSweep &statisticsSweep_,
                 const int &smearingType_,
                 PhononH0 &h0_);

  /** This function creates a "cache" of phonon properties.
   * To be called in the loop over k2, before the loop on k2.
   * If needed and have not been precomputed in the constructor, we compute the
   * phonon properties for all q for this list of k' at fixed k.
   *
   * @param k1: k1 wavevector in cartesian coordinates
   * @param k2Indexes: list of k-point indices running in the innerbandstructure
   */
  void prepare(const Eigen::Vector3d &k1, const std::vector<int> k2Indexes);

  /** This function returns the phonon properties at the point q such that
   * q = k1 - innerBandStructure.points(ik2)
   *
   * @paramm k1: first wavevector in outerbandstructure
   * @param ik2: index of second wavevector that lives in innerbandstructure
   * @return: a tuple with [qPointCartesianCoordinate,
   * qPointEnergies, qPointEnergies.size(), qPointEigenvectors,
   * qPointVelocities, boseEinsteinPopulation] of phonons for this q-point.
   */
  std::tuple<Eigen::Vector3d, Eigen::VectorXd, int, Eigen::MatrixXcd,
             Eigen::MatrixXd, Eigen::MatrixXd> get(Eigen::Vector3d &k1,
                                                   const int &ik2);

 private:
  BaseBandStructure &innerBandStructure;
  BaseBandStructure &outerBandStructure;
  StatisticsSweep &statisticsSweep;
  int smearingType;
  PhononH0 &h0;

  std::unique_ptr<BaseBandStructure> bandStructure3;
  std::unique_ptr<FullPoints> fullPoints3;
  std::unique_ptr<ActivePoints> activePoints3;

  bool storedAllQ3; // if true, q3 falls on a grid and we store the full bands
  const int storedAllQ3Case1=1;
  const int storedAllQ3Case2=2;
  int storedAllQ3Case;
  int cacheOffset=0;

  std::vector<Eigen::VectorXd> cacheEnergies;
  std::vector<Eigen::MatrixXcd> cacheEigvecs;
  std::vector<Eigen::MatrixXd> cacheBose;
  std::vector<Eigen::MatrixXd> cacheVelocity;
};

#endif
