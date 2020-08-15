#ifndef HELPER_EL_SCATTERING_H
#define HELPER_EL_SCATTERING_H

#include "phonon_h0.h"
#include "vector_bte.h"

class HelperElScattering {
 public:

  /** setup the class, and if possible, compute all 3rd points once.
   */
  HelperElScattering(BaseBandStructure &innerBandStructure_,
                 BaseBandStructure &outerBandStructure_,
                 StatisticsSweep &statisticsSweep_,
                 const int &smearingType_,
                 PhononH0 *h0_ = nullptr);

  /** to be called in the loop over q2, before the loop on q1.
   * If we have 3rd bandstructure, we don't do anything.
   * Otherwise, we compute all q3 for this list of q1
   */
  void prepare(const std::vector<long> k1Indexes, const Eigen::Vector3d &k2);

  /** to be called to get the info on q3.
   */
  std::tuple<Eigen::Vector3d, Eigen::VectorXd, int, Eigen::MatrixXcd,
             Eigen::MatrixXd, Eigen::MatrixXd> get(const long &ik1,
                                                   Eigen::Vector3d &k2);

 private:
  BaseBandStructure &innerBandStructure;
  BaseBandStructure &outerBandStructure;
  StatisticsSweep &statisticsSweep;
  int smearingType;
  PhononH0 *h0 = nullptr;

  std::unique_ptr<BaseBandStructure> bandStructure3;
  std::unique_ptr<FullPoints> fullPoints3;

  bool storedAllQ3; // if true, q3 falls on a grid and we store the full bands
  const int storedAllQ3Case1=1;
  const int storedAllQ3Case2=2;
  int storedAllQ3Case;

  std::vector<Eigen::VectorXd> cacheEnergies;
  std::vector<Eigen::MatrixXcd> cacheEigvecs;
  std::vector<Eigen::MatrixXd> cacheBose;
  std::vector<Eigen::MatrixXd> cacheVelocity;
};

#endif
