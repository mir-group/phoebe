#ifndef STATE_H
#define STATE_H

#include "exceptions.h"
#include "points.h"
#include "utilities.h"

/** Class containing harmonic information for all bands at a given k(q) point.
 * Detached refers to the fact that this object can be created just by knowing
 * the properties of the Bloch state, and doesn't come from any BandStructure
 * object. A BandStructureObject instead returns a State object.
 * Detached state copies and stores input energies, eigenvectors, velocities and
 * wavevectors coordinates, and is thus a little slower than State below.
 */
class DetachedState {
public:
  /** Class constructor
   * @param point: the 3d-vector with coordinates of the wavevector.
   * Prefer cartesian coordinates.
   * Nota bene: "point" is the only value that can be returned by getCoords(),
   * so, make sure you only need point in cartesian coordinates.
   * @param energies: vector with the energies of all bands at the given
   * wavevector.
   * @param numBands1: number of rows of the eigenvectors. e.g. for phonons
   * this is equal to 3*numAtoms, i.e. the full number of bands in the
   * harmonic hamiltonian.
   * @param numBands2: number of columns of the eigenvectors, and/or the
   * number of bands after applying the energy window.
   * @param velocities: complex tensor of size (numBands2, numBands2, 3) with
   * the values of the velocity operator of the given point.
   * @param eigenvectors: complex matrix of size (numBands1,numBands2) with
   * the values of the eigenvector of the given point.
   */
  DetachedState(Eigen::Vector3d &point_, Eigen::VectorXd &energies_,
                long numBands1_, long numBands2_,
                Eigen::MatrixXcd &eigenvectors_,
                Eigen::Tensor<std::complex<double>, 3> *velocities_ = nullptr);

  /** Copy constructor
   */
  DetachedState(const DetachedState &that);

  /** Copy assignment operator
   */
  DetachedState &operator=(const DetachedState &that);

  /** get the coordinates of a wavevector (Point object)
   * Note: even though we specify the basis, Detached state is only going to
   * return the same value as the "point" variable used in the constructor.
   * @param basis: Points::cartesianCoordinates.
   * @return point: an Eigen::Vector3d object with the point coordinates.
   */
  Eigen::Vector3d getCoords(const int &basis);

  /** get the energy of a single band
   * @param bandIndex: integer from 0 to numBands2-1
   * @return energy: Bloch state energy in rydbergs.
   */
  double getEnergy(const long &bandIndex, double chemicalPotential = 0.);

  /** get all energies at a given point
   * @param chemicalPotential: return energies relative to the chemical
   * potential, i.e. energy-chemicalPotential
   * @return energies: a vector of energies in rydbergs for all bands present
   */
  Eigen::VectorXd getEnergies(double chemicalPotential = 0.);

  /** get the group velocity of a single Bloch state.
   * @param bandIndex: integer from 0 to numBands2-1.
   * @return velocity: the 3d-vector of the group velocity.
   */
  Eigen::Vector3d getVelocity(const long &bandIndex);

  /** get the off-diagonal velocity operator of two Bloch states.
   * @param bandIndex1: bloch index of the bra state,
   * integer from 0 to numBands2-1.
   * @param bandIndex2: bloch index of the ket state,
   * integer from 0 to numBands2-1.
   * @return velocity: the complex 3d-vector of the velocity.
   */
  Eigen::Vector3cd getVelocity(const long &bandIndex1, const long &bandIndex2);

  /** get the group velocities of all bands for given k/q point
   * @return groupVelocities: double matrix of size (numBands2,3) with the
   * group velocities.
   */
  Eigen::MatrixXd getGroupVelocities();

  /** get the velocities of all bands for given k/q point
   * @return Velocities: a complex tensor of dimensions
   * (numBands2,numBands2,3) with the velocity operator.
   */
  Eigen::Tensor<std::complex<double>, 3> getVelocities();

  /** get the eigenvectors for the current Point
   * @input/output eigenvectors: a complex tensor of size
   * (3,numAtoms,numBands2) with the phonon eigenvectors. Error if
   * the eigenvectors are not set.
   */
  void getEigenvectors(Eigen::Tensor<std::complex<double>, 3> &eigs);

  /** get the eigenvectors for the current Point
   * @input/output eigenvectors: a complex matrix of size
   * (numBands1,numBands2) with the electron eigenvectors. Error if
   * eigenvectors are not set. numBands1 is the full number of bands in the
   * hamiltonian, numBands2 is the number of filtered bands.
   */
  void getEigenvectors(Eigen::MatrixXcd &eigs);

protected:
  // pointers to the bandstructure, I don't want to duplicate storage here
  Eigen::Vector3d point;
  Eigen::VectorXd energies;
  long numBands1;
  long numBands2;
  Eigen::Tensor<std::complex<double>, 3> velocities;
  Eigen::MatrixXcd eigenvectors;
};

/** Class containing harmonic information for all bands at a given k(q) point,
 * such as energies, wavevectors, velocities and eigenvectors.
 * This object doesn't store any vector, it is a collection of pointers and
 * references to the values stored, typically, in a BandStructure object.
 * Thus, it must be used without having references (the bandstructure) going
 * out of scope.
 *
 * Note: given that we use pointers, the State must know the order in which
 * velocity and eigenvectors are stored in the raw buffer (to reconstruct
 * tensors and matrices). Thus, if the storage order in BandStructure is
 * modified. so must be done for the State class get* methods.
 */
class State {
public:
  /** Class constructor
   * @param point: a Point instance with the coordinates in the brillouin
   * zone.
   * @param *energies: pointer to the beginning of the data with the
   * energies of the given point.
   * @param numBands1: number of rows of the eigenvectors. e.g. for phonons
   * this is equal to 3*numAtoms, i.e. the full number of bands in the harmonic
   * hamiltonian.
   * @param numBands2: number of columns of the eigenvectors, and/or the
   * number of bands after applying the energy window. For FullBandStructure
   * numBands1=numBands2, for ActiveBandStructure numBands1>=numBands2.
   * @param *velocities: pointer to the start of the vector containing the
   * values of the velocity operator of the given point.
   * @param *eigenvectors: pointer to the start of the vector containing the
   * values of the eigenvector of the given point.
   */
  State(Point &point_, double *energies_, long numBands1_, long numBands2_,
        std::complex<double> *velocities_ = nullptr,
        std::complex<double> *eigenvectors_ = nullptr);

  /** Copy constructor
   */
  State(const State &that);

  /** Copy assignment operator
   */
  State &operator=(const State &that);

  /** get the wavevector (as a Point object)
   * @return point: a Point object.
   */
  Point getPoint();

  /** get the coordinates of a wavevector (Point object)
   * @param basis: either Points::cartesianCoordinates or
   * Points::crystalCoordinates
   * @return point: an Eigen::Vector3d object with the point coordinates.
   */
  Eigen::Vector3d getCoords(const int &basis);

  /** get the energy of a single band
   * @param bandIndex: integer from 0 to numBands2-1
   * @return energy: Bloch state energy in rydbergs.
   */
  // TODO: remove this method?
  double getEnergy(const long &bandIndex, double chemicalPotential = 0.);

  /** get all energies at a given point
   * @param chemicalPotential: return energies relative to the chemical
   * potential, i.e. energy-chemicalPotential
   * @return energies: a vector of energies in rydbergs for all bands present
   */
  Eigen::VectorXd getEnergies(double chemicalPotential = 0.);

  /** get the group velocity of a single Bloch state.
   * @param bandIndex: integer from 0 to numBands2-1.
   * @return velocity: the 3d-vector of the group velocity.
   */
  Eigen::Vector3d getVelocity(const long &bandIndex);

  /** get the off-diagonal velocity operator of two Bloch states.
   * @param bandIndex1: bloch index of the bra state,
   * integer from 0 to numBands2-1.
   * @param bandIndex2: bloch index of the ket state,
   * integer from 0 to numBands2-1.
   * @return velocity: the complex 3d-vector of the velocity.
   */
  Eigen::Vector3cd getVelocity(const long &bandIndex1, const long &bandIndex2);

  /** get the group velocities of all bands for given k/q point
   * @return groupVelocities: double matrix of size (numBands2,3) with the
   * group velocities.
   */
  Eigen::MatrixXd getGroupVelocities();

  /** get the velocities of all bands for given k/q point
   * @return Velocities: a complex tensor of dimensions
   * (numBands2,numBands2,3) with the velocity operator.
   */
  Eigen::Tensor<std::complex<double>, 3> getVelocities();

  /** get the weight of the k/q point. Used for integrations over the
   * brillouin zone with an irreducible mesh of points.
   * @return weight: a vector of size (numBands2).
   */
  double getWeight();

  /** get the eigenvectors for the current Point
   * @input/output eigenvectors: a complex tensor of size
   * (3,numAtoms,numBands2) with the phonon eigenvectors. Error if
   * the eigenvectors are not set.
   */
  void getEigenvectors(Eigen::Tensor<std::complex<double>, 3> &eigs);

  /** get the eigenvectors for the current Point
   * @input/output eigenvectors: a complex matrix of size
   * (numBands1,numBands2) with the electron eigenvectors. Error if
   * eigenvectors are not set. numBands1 is the full number of bands in the
   * hamiltonian, numBands2 is the number of filtered bands.
   */
  void getEigenvectors(Eigen::MatrixXcd &eigs);

protected:
  // pointers to the bandstructure, I don't want to duplicate storage here
  Point point;
  double *energies;
  long numBands1; // this is full
  long numBands2; // this is reduced by the window
  std::complex<double> *velocities = nullptr;
  std::complex<double> *eigenvectors = nullptr;
  bool hasVelocities = false;
  bool hasEigenvectors = false;
};

#endif
