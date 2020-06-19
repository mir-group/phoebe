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
  long numTriplets;
  Eigen::Tensor<double, 4> ifc3Tensor;
  Eigen::Tensor<double, 3> cellPositions;
  Eigen::Tensor<long, 2> displacedAtoms;
  Eigen::MatrixXi tableAtCIndex1;
  Eigen::MatrixXi tableAtCIndex2;
  Eigen::MatrixXi tableAtCIndex3;

  Kokkos::View<double *****> D3_k;
  Kokkos::View<Kokkos::complex<double> ****> D3PlusCached_k, D3MinsCached_k;
  Kokkos::View<double **> cellPositions2_k, cellPositions3_k;

  template <typename A, typename B, typename C>
  std::tuple<Eigen::Tensor<double, 3>, Eigen::Tensor<double, 3>>
  calcCouplingSquared(A &state1, B &state2, C &state3Plus, C &state3Mins);

  // left here for future reference, but will be obsolete
  template <typename A, typename B, typename C>
  std::tuple<Eigen::Tensor<double, 3>, Eigen::Tensor<double, 3>>
  slowestCalcCouplingSquared(A &state1, B &state2, C &state3Plus,
                             C &state3Mins);

  // variables for the caching mechanism.
  // TODO: sparsify D3[+/-]Cached over the indices ir2 and ir3
  bool useD3Caching = true;
  Eigen::MatrixXd cellPositions2;
  Eigen::MatrixXd cellPositions3;
  Eigen::Tensor<double, 5> D3;
  long nr2, nr3, numAtoms, numBands;
  Eigen::Vector3d cachedCoords;
  Eigen::Tensor<std::complex<double>, 4> D3PlusCached;
  Eigen::Tensor<std::complex<double>, 4> D3MinsCached;

  typedef std::chrono::steady_clock::time_point time_point;
  typedef std::chrono::steady_clock::duration time_delta;

  std::vector<time_delta> dts;
  double tosec(time_delta dt) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(dt).count() /
           1e9;
  };

public:
  Interaction3Ph(
      Crystal &crystal_, long &numTriplets_,
      Eigen::Tensor<double, 4> &ifc3Tensor_,
      Eigen::Tensor<double, 3> &cellPositions_,
      Eigen::Tensor<long, 2> &displacedAtoms_);          // default constructor
  Interaction3Ph(const Interaction3Ph &that);            // copy constructor
  Interaction3Ph &operator=(const Interaction3Ph &that); // assignment op

  // this is an interface: we can compute it on the fly or read the cache.
  template <typename A, typename B, typename C>
  std::tuple<Eigen::Tensor<double, 3>, Eigen::Tensor<double, 3>>
  getCouplingSquared(A &state1, B &state2, C &state3Plus, C &state3Mins);

  ~Interaction3Ph() {
    std::cout << "calcCouplingSquared timing breakdown:"
              << "\n";
    std::cout << "nr2, nr3 phase loop: " << tosec(dts[0]) << "\n";
    std::cout << "D3Cached loop: " << tosec(dts[1]) << "\n";
    std::cout << "nr3 phase loop: " << tosec(dts[2]) << "\n";
    std::cout << "tmp loop: " << tosec(dts[3]) << "\n";
    std::cout << "tmp1 loop: " << tosec(dts[4]) << "\n";
    std::cout << "tmp2 loop: " << tosec(dts[5]) << "\n";
    std::cout << "vp loop: " << tosec(dts[6]) << "\n";
    std::cout << "vm loop: " << tosec(dts[7]) << "\n";
    std::cout << "cp loop: " << tosec(dts[8]) << "\n";
    std::cout << "cm loop: " << tosec(dts[9]) << "\n";
  }

  /**
   * Calculate single three-phonon matrix element (V^{+}/V^{-1}).
   *
   * Method for calculating the matrix element of a single, plus
   * or minus, three-phonon scattering process.
   *
   * @param[in] interactingPhonons phononTriplet data structure
   *  containing the three phonon modes (in the full Brillouin zone)
   *  taking part in an interaction event.
   * @param[in] phononEigenvectors Eigenvectors of all phonons in the
   *  full Brillouin zone.
   * @param[in] numTriplets Number of triplets considered in the supercell
   *  for the third order force constants calculation.
   * @param[in] ifc3Tensor Third order force constants.
   * @param[in] cellPositions Cartesian coordinates of the 2nd and 3rd
   * unitcells.
   * @param[in] displacedAtoms Index of the displaced atom for every triplet
   *  of unitcells.
   * @param[in] crystal Contains of crystal information
   * @param[in] procType Character '+' or '-' to choose the +/- type process.
   * @return The complex matrix element.
   */

  //	double calculateSingleV(const PhononTriplet &interactingPhonons, const
  // Eigen::MatrixXd &q, 			const int numTriplets, const
  // Eigen::Tensor<double,4> &ifc3Tensor, 			const
  // Eigen::Tensor<double,3> &cellPositions, 			const
  // Eigen::Tensor<int,2> &displacedAtoms,const CrystalInfo &crysInfo,
  // const char procType);

  //  void calculateAllVminus(const int *grid, const PhononMode &mode,
  //			    const Eigen::MatrixXd &qFBZ,
  //			    const Eigen::Tensor<complex<double>,3> &ev, const
  // int numTriplets, 			    const Eigen::Tensor<double,4>
  // &ifc3Tensor, 			    const Eigen::Tensor<double,3>
  // &cellPositions, 			    const Eigen::Tensor<int,2>
  //&displacedAtoms,const CrystalInfo &crysInfo);
  //
  //  void calculateAllW(const double T,const int *grid, const PhononMode &mode,
  //		     const Eigen::MatrixXi &indexMesh, const CrystalInfo
  //&crysInfo, 		     const Eigen::MatrixXd omega, const TetraData
  // tetra);
};

// template <typename A, typename B, typename C>
// std::tuple<Eigen::Tensor<double, 3>, Eigen::Tensor<double, 3>>
// Interaction3Ph::getCouplingSquared(A &state1, B &state2, C &state3Plus,
//                                    C &state3Mins) {
//   return calcCouplingSquared(state1, state2, state3Plus, state3Mins);
// }
//

/*
template <typename A, typename B, typename C>
std::tuple<Eigen::Tensor<double, 3>, Eigen::Tensor<double, 3>>
Interaction3Ph::slowestCalcCouplingSquared(A &state1, B &state2, C &state3Plus,
                                           C &state3Mins) {

  Eigen::Vector3d cell2Pos, cell3Pos;

  // Cartesian phonon wave vectors: q1,q2,q3
  auto q2 = state2.getCoords(Points::cartesianCoords);
  auto q3Plus = state3Plus.getCoords(Points::cartesianCoords);
  auto q3Mins = state3Mins.getCoords(Points::cartesianCoords);

  // phonon branches:
  // we allow the number of bands to be different in each direction
  long nb1 = state1.getNumBands(); // <- numBands
  long nb2 = state2.getNumBands();
  long nb3Plus = state3Plus.getNumBands();
  long nb3Mins = state3Mins.getNumBands();

  Eigen::Tensor<std::complex<double>, 3> vPlus(nb1, nb2, nb3Plus);
  Eigen::Tensor<std::complex<double>, 3> vMins(nb1, nb2, nb3Mins);
  vPlus.setZero();
  vMins.setZero();

  // this is the first implementation of this function
  // it's reasonable, but it's about 30 times slower than the
  // implementation below

  Eigen::Tensor<std::complex<double>, 3> ev1, ev2, ev3Plus, ev3Mins;

  state1.getEigenvectors(ev1); // size (3,numAtoms,numBands)
  state2.getEigenvectors(ev2);
  state3Plus.getEigenvectors(ev3Plus);
  state3Mins.getEigenvectors(ev3Mins);

  for (long it = 0; it < numTriplets; it++) { // sum over all triplets

    Eigen::Tensor<std::complex<double>, 3> v0Plus(nb1, nb2, nb3Plus);
    Eigen::Tensor<std::complex<double>, 3> v0Mins(nb1, nb2, nb3Mins);
    v0Plus.setZero();
    v0Mins.setZero();

    for (int ib1 = 0; ib1 < nb1; ib1++) {
      for (int ib2 = 0; ib2 < nb2; ib2++) {
        for (int ic3 : {0, 1, 2}) {
          for (int ic2 : {0, 1, 2}) {
            for (int ic1 : {0, 1, 2}) {
              // note: we should use ev2 without conjugation
              // but computed at -q2. Instead we use it at
              // +q2 and note that z^*(q)=z(-q)
              for (int ib3 = 0; ib3 < nb3Plus; ib3++) {
                v0Plus(ib1, ib2, ib3) +=
                    ifc3Tensor(it, ic1, ic2, ic3) *
                    ev1(ic1, displacedAtoms(it, 0), ib1) *
                    ev2(ic2, displacedAtoms(it, 1), ib2) *
                    std::conj(ev3Plus(ic3, displacedAtoms(it, 2), ib3));
              }

              for (int ib3 = 0; ib3 < nb3Mins; ib3++) {
                v0Mins(ib1, ib2, ib3) +=
                    ifc3Tensor(it, ic1, ic2, ic3) *
                    ev1(ic1, displacedAtoms(it, 0), ib1) *
                    std::conj(ev2(ic2, displacedAtoms(it, 1), ib2)) *
                    std::conj(ev3Mins(ic3, displacedAtoms(it, 2), ib3));
              }
            }
          }
        }
      }
    }

    for (int ic : {0, 1, 2}) {
      cell2Pos(ic) = cellPositions(it, 0, ic);
      cell3Pos(ic) = cellPositions(it, 1, ic);
    }

    // As a convention, the first primitive cell in the triplet is
    // restricted to the origin, so the phase for that cell is unity.

    double arg = +q2.dot(cell2Pos) - q3Plus.dot(cell3Pos);
    std::complex<double> phasePlus = exp(complexI * arg);
    arg = -q2.dot(cell2Pos) - q3Mins.dot(cell3Pos);
    std::complex<double> phaseMins = exp(complexI * arg);

    for (int ib1 = 0; ib1 < nb1; ib1++) {
      for (int ib2 = 0; ib2 < nb2; ib2++) {
        for (int ib3 = 0; ib3 < nb3Plus; ib3++) {
          // case +
          vPlus(ib1, ib2, ib3) += v0Plus(ib1, ib2, ib3) * phasePlus;
        }
        for (int ib3 = 0; ib3 < nb3Mins; ib3++) {
          // case -
          vMins(ib1, ib2, ib3) += v0Mins(ib1, ib2, ib3) * phaseMins;
        }
      }
    }
  }

  Eigen::Tensor<double, 3> couplingPlus(nb1, nb2, nb3Plus);
  Eigen::Tensor<double, 3> couplingMins(nb1, nb2, nb3Mins);
  // case +
  for (int ib3 = 0; ib3 < nb3Plus; ib3++) {
    for (int ib2 = 0; ib2 < nb2; ib2++) {
      for (int ib1 = 0; ib1 < nb1; ib1++) {
        couplingPlus(ib1, ib2, ib3) = std::norm(vPlus(ib1, ib2, ib3));
      }
    }
  }
  // case -
  for (int ib3 = 0; ib3 < nb3Mins; ib3++) {
    for (int ib2 = 0; ib2 < nb2; ib2++) {
      for (int ib1 = 0; ib1 < nb1; ib1++) {
        couplingMins(ib1, ib2, ib3) = std::norm(vMins(ib1, ib2, ib3));
      }
    }
  }
  return {couplingPlus, couplingMins};
}
*/

#endif
