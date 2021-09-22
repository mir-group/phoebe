#ifndef PH_INTERACTION_H
#define PH_INTERACTION_H

#include <Kokkos_Core.hpp>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>

#include "constants.h"
#include "crystal.h"
#include "eigen.h"
#include "points.h"
#include "utilities.h"

/** Class to calculate the probability rate for one 3-phonon scattering event.
 * In physical notation, this corresponds to the calculation of the term |V3|^2.
 * This class must be schematically used as follow:
 * for (iq2) {
 *   cacheD3()
 *   int numBatches = coupling3Ph->estimateNumBatches(nq1, nb2);
 *   for (batch : batches) {
 *     getCouplingsSquared()
 *     for ( iq1 : batch ) {
 *       ...
 *
 * The cacheD3() does a partial Fourier transform over the R2 Bravais lattice
 * vector. The list of all q1 wavevectors is split in batches, and the
 * getCouplingSquared() does the Fourier transform over the R1  Bravais lattice
 * vectors for all the wavevectors q1 and (q3 = q1 +- q2) in the batch.
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
 *
 * Use the environmental variable MAXMEM to set the amount of VRAM in gigabytes
 * available on the GPU/node (depending on the Kokkos installation).
 */
class Interaction3Ph {
private:
  Crystal &crystal_;

  // variables to be saved on the GPU
  Kokkos::View<double *****> D3_k;
  Kokkos::View<Kokkos::complex<double> ****> D3PlusCached_k, D3MinsCached_k;
  Kokkos::View<double **> cellPositions2_k, cellPositions3_k;
  Kokkos::View<double ***> weights2_k, weights3_k;

  double maxmem = 16.0e9; // default 16 Gb memory space for computation

  // dimensions
  int nr2, nr3, numAtoms, numBands;

public:

  /** Default constructor.
   * This method mostly moves data to the GPU if necessary.
   *
   * @param crystal: crystal object
   * @param D3: is a 5-dimensional real tensor with the third-derivative
   * anharmonic force constants. In details, it contains the tensor
   * D3(i,j,k,0,R',R"), where R are Bravais lattice vectors, and i, j, k are
   * indices running over both the Cartesian directions and the basis of atoms
   * in the primitive unit cell. See the theory section for a longer
   * description. The dimensions are
   * (3*numAtoms,3*numAtom,3*numAtom,numBravais,numBravais), where numBravais
   * must be consistent with the cellPositions chosen below.
   * @param cellPositions2: Bravais lattice vectors associated to each element
   * in the 4th index of D3 tensor.
   * @param cellPositions3: Bravais lattice vectors associated to each element
   * in the 5th index of D3 tensor.
   * @param weights2: weight of the Bravais lattice vector, i.e. a counter of
   * the degeneracy of symmetry-equivalent Bravais lattice vectors in
   * cellPositions2. Set to unity (e.g. in ShengBTE) if this is not used.
   * @param weights3: weight of the Bravais lattice vector, i.e. a counter of
   * the degeneracy of symmetry-equivalent Bravais lattice vectors in
   * cellPositions3. Set to unity (e.g. in ShengBTE) if this is not used.
   */
  Interaction3Ph(Crystal &crystal, Eigen::Tensor<double, 5> &D3,
                 Eigen::MatrixXd &cellPositions2,
                 Eigen::MatrixXd &cellPositions3, Eigen::VectorXd weights2,
                 Eigen::VectorXd weights3);

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
  getCouplingsSquared(const std::vector<Eigen::Vector3d> &q1s_e,
                      const Eigen::Vector3d &q2_e,
                      const std::vector<Eigen::MatrixXcd> &ev1s_e,
                      const Eigen::MatrixXcd &ev2_e,
                      const std::vector<Eigen::MatrixXcd> &ev3Pluss_e,
                      const std::vector<Eigen::MatrixXcd> &ev3Minss_e,
                      const std::vector<int> &nb1s_e, const int nb2,
                      const std::vector<int> &nb3Pluss_e,
                      std::vector<int> &nb3Minss_e);

  /** Computes a partial Fourier transform over the q2/R2 variables.
   * @param q2_e: values of the q2 cartesian coordinates over which the Fourier
   * transform is computed.
   */
  void cacheD3(const Eigen::Vector3d &q2_e);

  /** Estimate the number of batches that the list of q1 wavevectors must be
   * split into, in order to fit in memory.
   *
   * @param nq1: total number of q1 wavevectors to be split in batches
   * @param nb2: number of bands at the q2 wavevector.
   */
  int estimateNumBatches(const int &nq1, const int &nb2);
};

#endif
