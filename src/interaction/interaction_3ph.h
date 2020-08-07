#ifndef PHINTERACTION_H
#define PHINTERACTION_H

#include "constants.h"
#include "crystal.h"
#include "eigen.h"
#include "points.h"
#include "state.h"
#include "utilities.h"
#include <Kokkos_Core.hpp>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>

// * Class to calculate the coupling of a 3-phonon interaction given three
// * phonon modes.
class Interaction3Ph {
private:
  Crystal &crystal;

  Kokkos::View<double *****> D3_k;
  Kokkos::View<Kokkos::complex<double> ****> D3PlusCached_k, D3MinsCached_k;
  Kokkos::View<double **> cellPositions2_k, cellPositions3_k;

  typedef std::chrono::steady_clock::time_point time_point;
  typedef std::chrono::steady_clock::duration time_delta;

  std::vector<time_delta> dts;
  std::vector<time_delta> newdts;
  double tosec(time_delta dt) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(dt).count() /
           1e9;
  };

public:
  int nr2, nr3, numAtoms, numBands;
  Interaction3Ph(
      Crystal &crystal_, long &numTriplets,
      Eigen::Tensor<double, 4> &ifc3Tensor,
      Eigen::Tensor<double, 3> &cellPositions,
      Eigen::Tensor<long, 2> &displacedAtoms);          // default constructor
  Interaction3Ph(const Interaction3Ph &that);            // copy constructor
  Interaction3Ph &operator=(const Interaction3Ph &that); // assignment op

  std::tuple<std::vector<Eigen::Tensor<double, 3>>,
             std::vector<Eigen::Tensor<double, 3>>>
  getCouplingsSquared(std::vector<Eigen::Vector3d> q1s_e, Eigen::Vector3d q2_e,
                      std::vector<Eigen::MatrixXcd> ev1s_e,
                      Eigen::MatrixXcd ev2_e,
                      std::vector<Eigen::MatrixXcd> ev3Pluss_e,
                      std::vector<Eigen::MatrixXcd> ev3Minss_e,
                      std::vector<int> nb1s_e, int nb2,
                      std::vector<int> nb3Pluss_e, std::vector<int> nb3Minss_e);
  void cacheD3(Eigen::Vector3d q2_e);

  ~Interaction3Ph();
};

#endif
