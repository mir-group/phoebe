#ifndef HELPER_3RD_STATE_H
#define HELPER_3RD_STATE_H

#include "phonon_h0.h"
#include "vector_bte.h"

class Helper3rdState {
 public:

  /** setup the class, and if possible, compute all 3rd points once.
   */
  Helper3rdState(BaseBandStructure &innerBandStructure_,
                 BaseBandStructure &outerBandStructure_,
                 VectorBTE &outerBose_,
                 const int &smearingType_,
                 PhononH0 *h0_ = nullptr);

  /** to be called in the loop over q2, before the loop on q1.
   * If we have 3rd bandstructure, we don't do anything.
   * Otherwise, we compute all q3 for this list of q1
   */
  void prepare(const std::vector<long> q1Indexes, const long &iq2);

  /** to be called to get the info on q3.
   */
  std::tuple<Eigen::Vector3d, Eigen::VectorXd, long, Eigen::MatrixXcd,
             Eigen::MatrixXd, Eigen::MatrixXd> get(Point &point1,
                                                   Point &point2,
                                                   const int &thisCase);

  static const int casePlus;
  static const int caseMins;

 private:
  BaseBandStructure &innerBandStructure;
  BaseBandStructure &outerBandStructure;
  VectorBTE &outerBose;
  int smearingType;
  PhononH0 *h0 = nullptr;

  std::unique_ptr<BaseBandStructure> bandStructure3;
  std::unique_ptr<FullPoints> fullPoints3;

  bool storedAllQ3; // if true, q3 falls on a grid and we store the full bands

  const int storedAllQ3Case1=1;
  const int storedAllQ3Case2=2;
  int storedAllQ3Case;

  std::vector<Eigen::VectorXd> cachePlusEnergies, cacheMinsEnergies;
  std::vector<Eigen::MatrixXcd> cachePlusEigvecs, cacheMinsEigvecs;
  std::vector<Eigen::MatrixXd> cachePlusBose, cacheMinsBose;
  std::vector<Eigen::MatrixXd> cachePlusVelocity, cacheMinsVelocity;

};

#endif
