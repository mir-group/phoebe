#ifndef HARMONIC_H
#define HARMONIC_H

#include "bandstructure.h"
#include "particle.h"
#include "points.h"
#include "common_kokkos.h"

/** Virtual base class for Harmonic Hamiltonian.
 * The subclasses of this base class are the objects responsible for storing
 * the DFT harmonic properties computed on a coarse grid, and containing the
 * methods necessary for interpolating (diagonalizing) such harmonic
 * hamiltonian on a fine grid.
 */
class HarmonicHamiltonian {
 public:
  /** Returns the total number of phonon branches / electron bands that are
   * available in the interpolator.
   * @return numBands: integer number of bands.
   */
  virtual int getNumBands() = 0;

  /** Returns the Particle object to distinguish between electrons and phonons
   * @return particle: a Particle object
   */
  virtual Particle getParticle() = 0;

  /** Diagonalize the Harmonic Hamiltonian at an arbitrary wavevector
   * @param point: the Point object describing the wavevector.
   * @return eigenvalues: the values of quasiparticle energies, a double
   * vector of size numBands.
   * @return eigenvectors: a complex matrix of size (numBands,numBands) with
   * the eigenvectors, ordered in columns, computed at the input wavevector.
   */
  virtual std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> diagonalize(
      Point &point) = 0;

  /** Diagonalize the Harmonic Hamiltonian at an arbitrary wavevector.
   * Same as diagonalize(), but the wavevector coordinates are explicitly
   * passed in input.
   * @param point: the wavevector in cartesian coordinates.
   * @return eigenvalues: the values of quasiparticle energies, a double
   * vector of size numBands.
   * @return eigenvectors: a complex matrix of size (numBands,numBands) with
   * the eigenvectors, ordered in columns, computed at the input wavevector.
   */
  virtual std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
  diagonalizeFromCoordinates(
      Eigen::Vector3d &k) = 0;

  /** Computes the velocity operator (if possible, otherwise just its
   * diagonal matrix elements i.e. the group velocity) of the quasiparticle
   * at the input wavevector. v = \partial energy / \partial k.
   * @param point: Point object with the input wavevector.
   * @return velocity: a tensor of size(numBands,numBands,3) with the matrix
   * elements of the operator at the given wavevector. (third index is the
   * cartesian direction).
   */
  virtual Eigen::Tensor<std::complex<double>, 3> diagonalizeVelocity(
      Point &point) = 0;
  virtual Eigen::Tensor<std::complex<double>, 3>
  diagonalizeVelocityFromCoordinates(
      Eigen::Vector3d &coordinates) = 0;

  /** Method for the construction of the band structure on a grid of
   * of wavevectors of the Brillouin zone. In particular, we save the
   * information for all the bands available to the harmonic hamiltonian.
   * @param fullPoints: a Points object class over which we will compute the
   * band structure (e.g. FullPoints, or PathPoints)
   * @param withVelocities: bool. If set to true, we store the velocity
   * operator.
   * @param withEigenvectors: bool. If set to true, we store the eigenvectors
   * @return fullBandStructure: a FullBandStructure object that stores at
   * least the quasiparticle energies and, optionally, velocities and
   * eigenvectors.
   */
  virtual FullBandStructure populate(Points &fullPoints, const bool &withVelocities,
                                     const bool &withEigenvectors,
                                     const bool isDistributed = false) = 0;
  virtual std::tuple<DoubleView2D, StridedComplexView3D, ComplexView4D>
      kokkosBatchedDiagonalizeWithVelocities(
      const DoubleView2D &cartesianCoordinates) = 0;
  void kokkosBatchedTreatDegenerateVelocities(
      const DoubleView2D& cartesianCoordinates,
      const DoubleView2D& resultEnergies, ComplexView4D& resultVelocities,
      const double& threshold);
  virtual StridedComplexView3D kokkosBatchedBuildBlochHamiltonian(
      const DoubleView2D &cartesianCoordinates) = 0;
  virtual std::tuple<DoubleView2D, StridedComplexView3D>
  kokkosBatchedDiagonalizeFromCoordinates(
      const DoubleView2D &cartesianCoordinates) = 0;
};

#endif
