#ifndef PHINTERACTION_H
#define PHINTERACTION_H

#include "constants.h"
#include "crystal.h"
#include "eigen.h"
#include "points.h"
#include "state.h"
#include "utilities.h"
#include <Kokkos_Core.hpp>
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

template <typename A, typename B, typename C>
std::tuple<Eigen::Tensor<double, 3>, Eigen::Tensor<double, 3>>
Interaction3Ph::getCouplingSquared(A &state1, B &state2, C &state3Plus,
                                   C &state3Mins) {
  return calcCouplingSquared(state1, state2, state3Plus, state3Mins);
}

template <typename A, typename B, typename C>
std::tuple<Eigen::Tensor<double, 3>, Eigen::Tensor<double, 3>>
Interaction3Ph::calcCouplingSquared(A &state1, B &state2, C &state3Plus,
                                    C &state3Mins) {

  Eigen::Vector3d cell2Pos, cell3Pos;

  // Cartesian phonon wave vectors: q1,q2,q3
  auto q1 = state1.getCoords(Points::cartesianCoords);
  auto q2 = state2.getCoords(Points::cartesianCoords);

  Kokkos::View<double *> q1_k("q1", 3), q2_k("q2", 3);
  auto q1_h = Kokkos::create_mirror_view(q1_k),
       q2_h = Kokkos::create_mirror_view(q2_k);
  for (int i = 0; i < 3; i++) {
    q1_h(i) = q1(i);
    q2_h(i) = q2(i);
  }
  Kokkos::deep_copy(q1_k, q1_h);
  Kokkos::deep_copy(q2_k, q2_h);

  // phonon branches:
  // we allow the number of bands to be different in each direction
  long nb1 = state1.getNumBands(); // <- numBands
  long nb2 = state2.getNumBands();
  long nb3Plus = state3Plus.getNumBands();
  long nb3Mins = state3Mins.getNumBands();

  // Eigen::Tensor<std::complex<double>, 3> vPlus(nb1, nb2, nb3Plus);
  // Eigen::Tensor<std::complex<double>, 3> vMins(nb1, nb2, nb3Mins);
  // vPlus.setZero();
  // vMins.setZero();
  Kokkos::View<Kokkos::complex<double> ***> vPlus("vp", nb1, nb2, nb3Plus),
      vMins("vm", nb1, nb2, nb3Mins);

  // this is a bit more convoluted than the implementation above,
  // but it should be faster, as we break loops in two sections

  // Eigen::MatrixXcd ev1, ev2, ev3Plus, ev3Mins;

  Kokkos::View<Kokkos::complex<double> **> ev1("ev1", nb1, nb1),
      ev2("ev2", nb2, nb2), ev3Plus("ev3p", nb3Plus, nb3Plus),
      ev3Mins("ev3m", nb3Mins, nb3Mins);

  state1.getEigenvectors(ev1);
  state2.getEigenvectors(ev2);
  state3Plus.getEigenvectors(ev3Plus);
  state3Mins.getEigenvectors(ev3Mins);

  if (state2.getCoords(Points::cartesianCoords) != cachedCoords) {
    cachedCoords = state2.getCoords(Points::cartesianCoords);
    Kokkos::View<Kokkos::complex<double> **> phasePlus("pp", nr2, nr3),
        phaseMins("pm", nr2, nr3);

    Kokkos::parallel_for(
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {nr2, nr3}),
        KOKKOS_LAMBDA(int ir2, int ir3) {
          double argP = 0, argM = 0;
          for (int ic : {0, 1, 2}) {
            argP += +q2_k(ic) *
                    (cellPositions2_k(ir2, ic) - cellPositions3_k(ir3, ic));
            argM += -q2_k(ic) *
                    (cellPositions2_k(ir2, ic) - cellPositions3_k(ir3, ic));
          }
          phasePlus(ir2, ir3) = exp(complexI * argP);
          phaseMins(ir2, ir3) = exp(complexI * argM);
        });
    // std::cout << phasePlus(5, 8) << ", " << phaseMins(5, 8) << std::endl;

    Kokkos::parallel_for(
        Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
            {0, 0, 0, 0}, {numBands, numBands, numBands, nr3}),
        KOKKOS_LAMBDA(int ind1, int ind2, int ind3, int ir3) {
          Kokkos::complex<double> tmpp = 0, tmpm = 0;
          for (int ir2 = 0; ir2 < nr2; ir2++) { // sum over all triplets

            // As a convention, the first primitive cell in the triplet is
            // restricted to the origin, so the phase for that cell is unity.

            tmpp += D3_k(ind1, ind2, ind3, ir2, ir3) * phasePlus(ir2, ir3);
            tmpm += D3_k(ind1, ind2, ind3, ir2, ir3) * phaseMins(ir2, ir3);
          }
          D3PlusCached_k(ind1, ind2, ind3, ir3) = tmpp;
          D3MinsCached_k(ind1, ind2, ind3, ir3) = tmpm;
        });
  }

  //  std::cout << D3PlusCached_k(1, 1, 1, 1) << ", " << D3MinsCached_k(1, 1, 1,
  //  1)
  //            << std::endl;
  Kokkos::View<Kokkos::complex<double> *> phasePlus("pp", nr3),
      phaseMins("pm", nr3);

  Kokkos::parallel_for(
      nr3, KOKKOS_LAMBDA(int ir3) {
        double argP = 0, argM = 0;
        for (int ic : {0, 1, 2}) {
          argP += -q1_k(ic) * cellPositions3_k(ir3, ic);
          argM += -q1_k(ic) * cellPositions3_k(ir3, ic);
        }
        phasePlus(ir3) = exp(complexI * argP);
        phaseMins(ir3) = exp(complexI * argM);
      });

  // As a convention, the first primitive cell in the triplet is
  // restricted to the origin, so the phase for that cell is unity.

  // note: tmp* is a tensor over cartesian and atomic indices
  // (whose size coincides with the band numbers)
  // Eigen::Tensor<std::complex<double>, 3> tmpPlus(numBands, numBands,
  // numBands); Eigen::Tensor<std::complex<double>, 3> tmpMins(numBands,
  // numBands, numBands); tmpPlus.setZero(); tmpMins.setZero();
  Kokkos::View<Kokkos::complex<double> ***> tmpPlus("tmpp", numBands, numBands,
                                                    numBands),
      tmpMins("tmpm", numBands, numBands, numBands);
  Kokkos::parallel_for(
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0},
                                             {numBands, numBands, numBands}),
      KOKKOS_LAMBDA(int iac1, int iac2, int iac3) {
        Kokkos::complex<double> tmpp = 0, tmpm = 0;
        for (int ir3 = 0; ir3 < nr3; ir3++) { // sum over all triplets
          tmpp += D3PlusCached_k(iac1, iac2, iac3, ir3) * phasePlus(ir3);
          tmpm += D3MinsCached_k(iac1, iac2, iac3, ir3) * phaseMins(ir3);
        }
        tmpPlus(iac1, iac2, iac3) = tmpp;
        tmpMins(iac1, iac2, iac3) = tmpm;
      });

  // now we want to multiply
  // vPlus(ib1,ib2,ib3) = tmpPlus(iac1,iac2,iac3) * ev1(iac1,ib1)
  //          * std::conj(ev2(iac2,ib2)) * std::conj(ev3Plus(iac3,ib3));
  // vMins(ib1,ib2,ib3) += tmpMins(iac1,iac2,iac3) * ev1(iac1,ib1)
  //			* std::conj(ev2(iac2,ib2)) *
  // std::conj(ev3Mins(iac3,ib3));
  // a single loop over all 6 indices at once is easier to read, but slow
  // instead, we break it into 3 loops with 4 nested loops each
  // ( ~numBands^2 faster)

  // Eigen::Tensor<std::complex<double>, 3> tmp1Plus(nb1, numBands, numBands);
  // Eigen::Tensor<std::complex<double>, 3> tmp1Mins(nb1, numBands, numBands);
  // tmp1Plus.setZero();
  // tmp1Mins.setZero();
  Kokkos::View<Kokkos::complex<double> ***> tmp1Plus("t1p", nb1, numBands,
                                                     numBands),
      tmp1Mins("t1m", nb1, numBands, numBands);

  Kokkos::parallel_for(
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0},
                                             {nb1, numBands, numBands}),
      KOKKOS_LAMBDA(int ib1, int iac2, int iac3) {
        Kokkos::complex<double> tmpp = 0, tmpm = 0;
        for (int iac1 = 0; iac1 < numBands; iac1++) {
          tmpp += tmpPlus(iac1, iac2, iac3) * ev1(iac1, ib1);
          tmpm += tmpMins(iac1, iac2, iac3) * ev1(iac1, ib1);
        }
        tmp1Plus(ib1, iac2, iac3) = tmpp;
        tmp1Mins(ib1, iac2, iac3) = tmpm;
      });
  // Eigen::Tensor<std::complex<double>, 3> tmp2Plus(nb1, nb2, numBands);
  // Eigen::Tensor<std::complex<double>, 3> tmp2Mins(nb1, nb2, numBands);
  // tmp2Plus.setZero();
  // tmp2Mins.setZero();
  Kokkos::View<Kokkos::complex<double> ***> tmp2Plus("t2p", nb1, nb2, numBands),
      tmp2Mins("t2m", nb1, nb2, numBands);

  Kokkos::parallel_for(
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {nb1, nb2, numBands}),
      KOKKOS_LAMBDA(int ib1, int ib2, int iac3) {
        Kokkos::complex<double> tmpp = 0, tmpm = 0;
        for (int iac2 = 0; iac2 < numBands; iac2++) {
          tmpp += tmp1Plus(ib1, iac2, iac3) * ev2(iac2, ib2);
          tmpm += tmp1Mins(ib1, iac2, iac3) * Kokkos::conj(ev2(iac2, ib2));
        }
        tmp2Plus(ib1, ib2, iac3) = tmpp;
        tmp2Mins(ib1, ib2, iac3) = tmpm;
      });

  // the last loop is split in two because nb3 is not the same for + and -
  Kokkos::parallel_for(
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {nb1, nb2, nb3Plus}),
      KOKKOS_LAMBDA(int ib1, int ib2, int ib3) {
        Kokkos::complex<double> tmpp = 0;
        for (int iac3 = 0; iac3 < numBands; iac3++) {
          tmpp += tmp2Plus(ib1, ib2, iac3) * Kokkos::conj(ev3Plus(iac3, ib3));
        }
        vPlus(ib1, ib2, ib3) = tmpp;
      });
  Kokkos::parallel_for(
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {nb1, nb2, nb3Mins}),
      KOKKOS_LAMBDA(int ib1, int ib2, int ib3) {
        Kokkos::complex<double> tmpm = 0;
        for (int iac3 = 0; iac3 < numBands; iac3++) {
          tmpm += tmp2Mins(ib1, ib2, iac3) * Kokkos::conj(ev3Mins(iac3, ib3));
        }
        vMins(ib1, ib2, ib3) = tmpm;
      });

  Eigen::Tensor<double, 3> couplingPlus(nb1, nb2, nb3Plus);
  auto vPlus_h = Kokkos::create_mirror_view(vPlus);
  Kokkos::deep_copy(vPlus_h, vPlus);
  auto norm = [](auto x) { return x.real() * x.real() + x.imag() * x.imag(); };
// case +
#pragma omp parallel for collapse(3)
  for (int ib1 = 0; ib1 < nb1; ib1++) {
    for (int ib2 = 0; ib2 < nb2; ib2++) {
      for (int ib3 = 0; ib3 < nb3Plus; ib3++) {
        couplingPlus(ib1, ib2, ib3) = norm(vPlus_h(ib1, ib2, ib3));
      }
    }
  }
  // case -
  Eigen::Tensor<double, 3> couplingMins(nb1, nb2, nb3Mins);
  auto vMins_h = Kokkos::create_mirror_view(vMins);
  Kokkos::deep_copy(vMins_h, vMins);
#pragma omp parallel for collapse(3)
  for (int ib1 = 0; ib1 < nb1; ib1++) {
    for (int ib2 = 0; ib2 < nb2; ib2++) {
      for (int ib3 = 0; ib3 < nb3Mins; ib3++) {
        couplingMins(ib1, ib2, ib3) = norm(vMins_h(ib1, ib2, ib3));
      }
    }
  }
  return {couplingPlus, couplingMins};
}

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

#endif
