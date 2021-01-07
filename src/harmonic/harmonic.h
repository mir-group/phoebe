#ifndef HARMONIC_H
#define HARMONIC_H

#include "bandstructure.h"
#include "particle.h"
#include "points.h"

/** Base class for the Harmonic Hamiltonian.
 * The subclasses of this base class are the objects responsible for storing
 * the DFT harmonic properties computed on a coarse grid, and containing the
 * methods necessary for interpolating (diagonalizing) such harmonic
 * hamiltonian on a fine grid.
 */
class HarmonicHamiltonian {
 public:
  /** Default constructor
   */
  HarmonicHamiltonian();

  /** Returns the total number of phonon branches / electron bands that are
   * available in the interpolator.
   * @return numBands: integer number of bands.
   */
  virtual int getNumBands();

  /** Returns the Particle object to distinguish between electrons and phonons
   * @return particle: a Particle object
   */
  virtual Particle getParticle();

  /** Diagonalize the Harmonic Hamiltonian at an arbitrary wavevector
   * @param point: the Point object describing the wavevector.
   * @return eigenvalues: the values of quasiparticle energies, a double
   * vector of size numBands.
   * @return eigenvectors: a complex matrix of size (numBands,numBands) with
   * the eigenvectors, ordered in columns, computed at the input wavevector.
   */
  virtual std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> diagonalize(
      Point &point);

  /** Diagonalize the Harmonic Hamiltonian at an arbitrary wavevector.
   * Same as diagonalize(), but the wavevector coordinates are explicitly
   * passed in input.
   * @param point: the wavevector in cartesian coordinates.
   * @return eigenvalues: the values of quasiparticle energies, a double
   * vector of size numBands.
   * @return eigenvectors: a complex matrix of size (numBands,numBands) with
   * the eigenvectors, ordered in columns, computed at the input wavevector.
   */
  virtual std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> diagonalizeFromCoords(
      Eigen::Vector3d &k);

  /** Computes the velocity operator (if possible, otherwise just its
   * diagonal matrix elements i.e. the group velocity) of the quasiparticle
   * at the input wavevector. v = \partial energy / \partial k.
   * @param point: Point object with the input wavevector.
   * @return velocity: a tensor of size(numBands,numBands,3) with the matrix
   * elements of the operator at the given wavevector. (third index is the
   * cartesian direction).
   */
  virtual Eigen::Tensor<std::complex<double>, 3> diagonalizeVelocity(
      Point &point);
  virtual Eigen::Tensor<std::complex<double>, 3> diagonalizeVelocityFromCoords(
      Eigen::Vector3d &coords);

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
  virtual FullBandStructure populate(Points &fullPoints, bool &withVelocities,
                                     bool &withEigenvectors,
                                     bool isDistributed = false);

 protected:
  Particle particle;

  int numBands = 0;
};

#endif
