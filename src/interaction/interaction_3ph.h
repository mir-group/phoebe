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

/** Class to calculate the probability rate for one 3-phonon scattering event.
 * In physical notation, this corresponds to the calculation of the term |V3|^2.
 * This class must be schematically used as follow:
 * for (iq2) {
 *   cacheD3()
 *   for (iq1) {
 *     getCouplingsSquared()
 *
 * The cacheD3() does a partial Fourier transform over the R2 Bravais lattice
 * vector, while getCouplingSquared() does the Fourier transform over the R1
 * Bravais lattice vectors for a bunch of q1 and (q3 = q1 +- q2) wavevectors.
 *
 * This calculation is GPU optimized. We found that the R1 fourier-transform
 * should be done for a bunch of q1/q3 wavevectors in order to optimize the
 * usage of a GPU (calculations would be too fast for a single q1 wavevector,
 * and most of the time would be spent in data transfer latency). We use Kokkos
 * as the library for GPU acceleration, so that this class can also revert back
 * to a CPU-only calculation if the GPU is not found at compilation time.
 *
 * The class assumes that the coupling matrix elements in real space passed to
 * input to phoebe are in the form D(R1,R2,R3) -> D(0,R2,R3), i.e. we take
 * advantage of the crystal periodicity to set the phases of R1 to zero, so that
 * a 6-dimensional Fourier transform is enough (and avoid a 9D FT).
 */
class Interaction3Ph {
private:
  Crystal &crystal_;

  // variables to be saved on the GPU
  Kokkos::View<double *****> D3_k;
  Kokkos::View<Kokkos::complex<double> ****> D3PlusCached_k, D3MinsCached_k;
  Kokkos::View<double **> cellPositions2_k, cellPositions3_k;

  // variables for timing
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

  /** Default constructor.
   * This method mostly reshapes the input into something more manageable and
   * moves data to the GPU if necessary.
   *
   * @param crystal: crystal object
   * @param ifc3Tensor: tensor with the third-derivative anharmonic force
   * constants in a format similar to that of ShengBTE (may be changed in the
   * future), which is provided by the parser function. This content is
   * processed in the constructor and reshaped into something more convenient
   * for execution. In input is shaped as (3*numAtoms,3*numAtoms,3*numAtoms,
   * numtriplets), where numTriplets is an index over the triplets of atoms
   * displaced to compute the 3rd derivative, and the first three indices run
   * over 3 cartesian coordinates and the atomic basis, indicating which atoms
   * have been moved to compute the tensor element. Values stored are the value
   * of energy change upon such displacements of triplets, in Rydbergs.
   * @param cellPositions: Bravais lattice vectors associated to each element
   * (triplet) of the ifc3tensor.
   * @param displacedAtoms: indices of the atoms that are displaced for each
   * element (triplet) of the Ifc3Tensor.
   */
  Interaction3Ph(
      Crystal &crystal_, long &numTriplets,
      Eigen::Tensor<double, 4> &ifc3Tensor,
      Eigen::Tensor<double, 3> &cellPositions,
      Eigen::Tensor<long, 2> &displacedAtoms);

  /** Copy constructor
   */
  Interaction3Ph(const Interaction3Ph &that);

  /** Assignment operator
   */
  Interaction3Ph &operator=(const Interaction3Ph &that);

  /** Computes the |V3|^2 matrix elements for a bunch of q1 wavevectors at fixed
   * q2 wavevector.
   *
   * @param q1s_e: list of wavevectors at q1 in cartesian coordinates
   * @param q2_e: wavevector at q2 in cartesian coordinates
   * @param ev1s_e: list of eigenvectors at all the q1 points
   * @param ev2_e: eigenvector at q2
   * @param ev3Pluss_e: list of eigenvectors at all the q3 = q1 + q2 points
   * @param ev3Minss_e: list of eigenvectors at all the q3 = q1 - q2 points
   * @param nb1s_e: list of number of bands at q1
   * @parma nb2: number of bands at q2
   * @param nb3Pluss_e: list of number of bands at q3 = q1 + q2
   * @param nb3Minss_e: list of number of bands at q3 = q1 - q2
   *
   * Note: we pass the number of bands as a mechanism to cope with the
   * reduction of bands to only the "active" ones, ,rather than always using all
   * 3*numAtoms phonon bands.
   */
  std::tuple<std::vector<Eigen::Tensor<double, 3>>,
             std::vector<Eigen::Tensor<double, 3>>>
  getCouplingsSquared(std::vector<Eigen::Vector3d> q1s_e, Eigen::Vector3d q2_e,
                      std::vector<Eigen::MatrixXcd> ev1s_e,
                      Eigen::MatrixXcd ev2_e,
                      std::vector<Eigen::MatrixXcd> ev3Pluss_e,
                      std::vector<Eigen::MatrixXcd> ev3Minss_e,
                      std::vector<int> nb1s_e, int nb2,
                      std::vector<int> nb3Pluss_e, std::vector<int> nb3Minss_e);

  /** Computes a partial Fourier transform over the q2/R2 variables.
   * @param q2_e: values of the q2 cartesian coordinates over which the Fourier
   * transform is computed.
   */
  void cacheD3(Eigen::Vector3d q2_e);

  /** Destructor method, which prints a summary of execution times.
   */
  ~Interaction3Ph();
};

#endif
