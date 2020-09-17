#ifndef BANDSTRUCTURE_H
#define BANDSTRUCTURE_H

#include "exceptions.h"
#include "particle.h"
#include "points.h"
#include "utilities.h"
#include "Matrix.h"

/** Base class for describing objects containing the band structure, i.e.
 * the harmonic properties of a quasiparticle, as a function of wavevectors
 * (Points).
 */
class BaseBandStructure {
 public:
  /** Get the Particle object associated with this class
   * @return particle: a Particle object, describing e.g. whether this
   * is a phonon or electron bandStructure
   */
  virtual Particle getParticle() = 0;

  /** Returns the wavevectors on which the bandstructure is computed.
   * @return Points: the object representing the Brillouin zone wavevectors.
   */
  virtual Points getPoints() = 0;

  /** Returns a wavevector, given a wavevector index.
   * The wavevector index runs from 0 to numPoints-1
   */
  virtual Point getPoint(const long &pointIndex) = 0;

  /** Returns the total number of k/q-points.
   * @param useFullGrid: default = false. If true, returns the number of
   * points of the full monkhorst-pack grid of wavevectors, otherwise,
   * returns the number of wavevectors stored in the BandStructure.
   * @return numPoints: the total number of wavevectors of the bandStructure.
   */
  virtual long getNumPoints(const bool &useFullGrid = false) = 0;

  /** Returns the number of bands.
   * @return numPoints: the total number of wavevectors of the bandStructure.
   * If the number of bands is not constant, calls an error.
   */
  virtual long getNumBands() = 0;

  /** Checks whether the bandStructure has been built discarding some Bloch
   * states from those available.
   * @return windowMethod: one of the values of Window::filterMethod. 0 for
   * no filter.
   */
  virtual long hasWindow() = 0;
  /** Returns the boolean determining if this bandstructure is distributed
   * or not.
   */
  virtual bool getIsDistributed() = 0;

  // needed in the BTE
  /** Builds a Bloch state index, which combines both wavevector index and
   * band index.
   * It's used to view the various matrices such as energy as a 1D vector,
   * and can be used in combination with get() methods.
   * @param wavevectorIndex: strong-typed index on wavevector
   * @return stateIndex: integer from 0 to numStates-1=numBands*numPoints-1
   */
  virtual long getIndex(const WavevectorIndex &ik, const BandIndex &ib) = 0;

  /** Given a Bloch state index, finds the corresponding wavevector and band
   * index.
   * @param stateIndex: integer from 0 to numStates-1=numBands*numPoints-1
   * @return WavevectorIndex: strong-typed index on wavevector
   * @return BandIndex: strong-typed index on bands
   */
  virtual std::tuple<WavevectorIndex, BandIndex> getIndex(const long &is) = 0;
  virtual std::tuple<WavevectorIndex, BandIndex> getIndex(StateIndex &is) = 0;

  /** Returns the total number of Bloch states.
   * @return numStates: the integer number of Bloch states.
   */
  virtual long getNumStates() = 0;

  /** Returns an iterator to be used for loops over the Bloch state index.
   * The values of the iterator are distributed in N blocks over N MPI ranks.
   */
  std::vector<int> parallelStateIterator();

  /** Returns the energy of a quasiparticle from its Bloch index
   * Used for accessing the bandstructure in the BTE.
   * @param stateIndex: an integer index in range [0,numStates[
   * @return energy: the value of the QP energy for that given Bloch index.
   * Phonon energies are referred to zero, with negative energies being
   * actually complex phonon frequencies. Electronic energies are not saved
   * with any particular reference, and should be used together with the
   * chemical potential computed by StatisticsSweep. By policy, it's in
   * rydbergs units.
   */
  virtual const double &getEnergy(const long &stateIndex) = 0;
  virtual const double &getEnergy(StateIndex &is) = 0;
  virtual Eigen::VectorXd getEnergies(WavevectorIndex &ik) = 0;

  /** Returns the energy of a quasiparticle from its Bloch index
   * Used for accessing the bandstructure in the BTE.
   * @param stateIndex: an integer index in range [0,numStates[
   * @return velocity: a 3d vector with velocity. By policy, we save it in
   * the cartesian basis and in atomic rydberg units.
   */
  virtual Eigen::Vector3d getGroupVelocity(const long &stateIndex) = 0;
  virtual Eigen::Vector3d getGroupVelocity(StateIndex &is) = 0;
  virtual Eigen::MatrixXd getGroupVelocities(WavevectorIndex &ik) = 0;
  virtual Eigen::Tensor<std::complex<double>, 3> getVelocities(
      WavevectorIndex &ik) = 0;

  virtual Eigen::MatrixXcd getEigenvectors(WavevectorIndex &ik) = 0;
  virtual Eigen::Tensor<std::complex<double>, 3> getPhEigenvectors(
      WavevectorIndex &ik) = 0;

  /** Returns the energy of a quasiparticle from its Bloch index
   * Used for accessing the bandstructure in the BTE.
   * @param stateIndex: an integer index in range [0,numStates[
   * @return wavevector: a 3d vector with the wavevector in cartesian
   * coordinates in units of Bohr^-1.
   */
  virtual Eigen::Vector3d getWavevector(const long &stateIndex) = 0;
  virtual Eigen::Vector3d getWavevector(StateIndex &is) = 0;
  virtual Eigen::Vector3d getWavevector(WavevectorIndex &ik) = 0;

  virtual double getWeight(const long &stateIndex) = 0;
  virtual double getWeight(StateIndex &is) = 0;
  virtual double getWeight(WavevectorIndex &ik) = 0;

  /** Method to save quasiparticle eigenvectors inside FullBandStructure().
   * @param point: a vector of 3 crystal coordinates. The method will look
   * for the wavevector index.
   * @param energies: a vector of size (numBands) with the quasiparticle
   * energies
   */
  virtual void setEnergies(Point &point, Eigen::VectorXd &energies_) = 0;

  /** Method to save quasiparticle eigenvectors inside FullBandStructure().
   * Note that in this case, eigenvectors are passed as a matrix, which is
   * the case e.g. for the Wannier interpolation, where the eigenvectors
   * represent the unitary transformation matrix U for Wannier localization.
   * @param point: a Point object with the coordinates of the wavevector,
   * which should come from the same Point class stored in FullBandStructure
   * @param eigenvectors: a complex matrix of size (numBands,numBands)
   */
  virtual void setEigenvectors(Point &point,
                               Eigen::MatrixXcd &eigenvectors_) = 0;

  /** Saves in the class the velocities computed at a particular point.
   * @param point: a Point object representing the wavevector where these
   * velocities have been computed.
   * @param velocities: a rank-3 tensor of size (numBands,numBands,3)
   * containing the matrix elements of the velocity operator. Diagonal
   * elements are the quasiparticle group velocities.
   */
  virtual void setVelocities(
      Point &point, Eigen::Tensor<std::complex<double>, 3> &velocities_) = 0;
};

class ActiveBandStructure;
// forward declaration of friend class
// we could remove this by making ActiveBS a subclass of FullBS

/** FullBandStructure is the class that stores the energies, velocities and
 * eigenvectors of a quasiparticle computed on a set of wavevectors (as defined
 * by Points() ) in the Brillouin zone.
 * By default, each MPI process holds a full copy of the bandstructure.
 * However, the bandstructure can be distributed over the wavevectors, if so
 * specified in the constructor.
 *
 * An important note for developers: When using a distributed bandstructure, 
 * looping over numStates of the bandstructure will not work -- you need to 
 * loop over the iterator of indices provided by getStateIndices or 
 * getWavevectorIndices. The class will throw errors when nonlocal values are 
 * being requested. All functions in bandstructure are written to take the 
 * global wavevector indices associated with the Points object internal to the 
 * bandstructure (because we use the Point class to find wavevector indices in
 * get and set functions of bandstructure).
 */
class FullBandStructure : public BaseBandStructure {
 public:
  /** Constructor of the FullBandStructure
   * @param numBands: an integer with the number of bands in the system
   * @param particle: a Particle object that contains the type of
   * quasiparticle
   * @param withVelocities: a boolean to decide whether to store velocities
   * @param withEigenvectors: a boolean to decide whether to store eigenvecs
   * @param points: the underlying mesh of wavevectors.
   * @param isDistributed: if true, we distribute in memory the
   * storage of quantities, parallelizing over wavevectors (points). By default
   * we don't distribute data in parallel.
   */
  FullBandStructure(long numBands_, Particle &particle_, bool withVelocities,
                    bool withEigenvectors, Points &points_,
                    bool isDistributed_ = false);

  FullBandStructure();

  /** Copy constructor
   */
  FullBandStructure(const FullBandStructure &that);

  /** Assignment operator
   */
  FullBandStructure &operator=(const FullBandStructure &that);

  /** Get the Particle object associated with this class
   * @return particle: a Particle object, describing e.g. whether this
   * is a phonon or electron bandStructure
   */
  Particle getParticle();

  /** Returns the wavevectors on which the bandstructure is computed.
   * @return Points: the object representing the Brillouin zone wavevectors.
   */
  Points getPoints();

  /** Returns a wavevector, given a wavevector index.
   * The wavevector index runs from 0 to numPoints-1
   */
  Point getPoint(const long &pointIndex);

  /** Returns the total number of k/q-points.
   * @param useFullGrid: used for compatibility with the ActiveBandStructure,
   * it doesn't affect the returned value.
   * @return numPoints: the total number of wavevectors of the bandStructure.
   */
  long getNumPoints(const bool &useFullGrid = false);

  /** Returns the number of bands.
   * @return numPoints: the total number of wavevectors of the bandStructure.
   */
  long getNumBands();

  long hasWindow();

  bool getIsDistributed();

  /** Builds a Bloch state index, which runs on both wavevector index and
   * band index. ik runs from 0 to numPoints-1, ib from 0 to numBands-1.
   * It's used to view the various matrices such as energy as a 1D vector,
   * and can be used in combination with get() methods.
   * @param wavevectorIndex: strong-typed index on wavevector
   * @return stateIndex: integer from 0 to numStates-1=numBands*numPoints-1
   */
  long getIndex(const WavevectorIndex &ik, const BandIndex &ib);

  /** Given a Bloch state index, finds the corresponding wavevector and band
   * index.
   * @param stateIndex: integer from 0 to numStates-1=numBands*numPoints-1
   * @return WavevectorIndex: strong-typed index on wavevector
   * @return BandIndex: strong-typed index on bands
   */
  std::tuple<WavevectorIndex, BandIndex> getIndex(const long &is);

  /** Given a Bloch state index, finds the corresponding wavevector and band
   * index.
   * @param stateIndex: StateIndex(is) object where is is an integer from 0 to
   * numStates-1=numBands*numPoints-1
   * @return WavevectorIndex: strong-typed index on wavevector
   * @return BandIndex: strong-typed index on bands
   */
  std::tuple<WavevectorIndex, BandIndex> getIndex(StateIndex &is);

  /** Returns the total number of Bloch states, equal to numPoints*numBands.
   * @return numStates: the total number of Bloch states in the class.
   */
  long getNumStates();

  /** Returns the indices of all wavevector indices on this process, or in
  * an undistributed case, returns all wavevector indices.
  * @return wavevectorIndices: a vector of wavevector indices for use as
  * an iterator.
  */
  std::vector<long> getWavevectorIndices();

  /** Returns the indices of all state indices on this process, or in
  * an undistributed case, returns all state indices.
  * @return stateIndices: a vector of tuples of (ib,ik) indices for use as
  * an iterator.
  */
  std::vector<std::tuple<long,long>> getStateIndices();

  /** Returns the indices of all bands on this process, or in an
  * undistributed case, returns all band indices.
  * @return bandIndices: a vector of band indices for use as an iterator.
  */
  std::vector<long> getBandIndices();

  // TODO in the distributed case, this might not be true!
  // then, is will not run from [0,numStates].
  /** Returns the energy of a quasiparticle from its Bloch index.
   * Used for accessing the bandstructure in the BTE.
   * @param stateIndex: an integer index in range [0,numStates[
   * @return energy: the value of the QP energy for that given Bloch index.
   * Phonon energies are referred to zero, with negative energies being
   * actually complex phonon frequencies. Electronic energies are not saved
   * with any particular reference, and should be used together with the
   * chemical potential computed by StatisticsSweep. By policy, it's in
   * rydbergs units.
   */
  const double &getEnergy(const long &stateIndex);

  // TODO in the distributed case, this might not be true!
  // then, is will not run from [0,numStates-1].
  /** Returns the energy of a quasiparticle from its Bloch index.
   * Same as getEnergy(const long &stateIndex), but using a StateIndex input
   * @param stateIndex: a StateIndex(is) object where 'is' is an integer
   * running over the number of states [0,numStates-1].
   * @return energy: the value of the QP energy for that given Bloch index.
   * Phonon energies are referred to zero, with negative energies being
   * actually complex phonon frequencies. Electronic energies are not saved
   * with any particular reference, and should be used together with the
   * chemical potential computed by StatisticsSweep. By policy, it's in
   * rydbergs units.
   */
  const double &getEnergy(StateIndex &is);

  /** Returns the energy of a quasiparticle from its band and wavevector index.
   * Same as getEnergy(const long &stateIndex), but using ib,ik instead.
   * ik should always be the global wavevector index, or this will be wrong!
   * @param ik: the wavevector index of the particle state
   * @param ib: the band index of the particle state
   * @return energy: the value of the QP energy for that given Bloch index.
   * Phonon energies are referred to zero, with negative energies being
   * actually complex phonon frequencies. Electronic energies are not saved
   * with any particular reference, and should be used together with the
   * chemical potential computed by StatisticsSweep. By policy, it's in
   * rydbergs units.
   */
  const double &getEnergy(WavevectorIndex &ik, BandIndex &ib);

  /** Returns the energies of all quasiparticle computed at a specified
   * wavevector.
   * @param wavevectorIndex: a WavevectorIndex(ik) where ik is an integer
   * index running over the wavevectors [0,numPoints-1]
   * @return energies: an Eigen vector of quasiparticle energies for that
   * given wavevector.
   * Phonon energies are referred to zero, with negative energies being
   * actually complex phonon frequencies. Electronic energies are not saved
   * with any particular reference, and should be used together with the
   * chemical potential computed by StatisticsSweep. In rydbergs units.
   */
  Eigen::VectorXd getEnergies(WavevectorIndex &ik);

  /** Returns the group velocity of a quasiparticle from its Bloch index.
   * Used for accessing the bandstructure in the BTE.
   * @param stateIndex: an integer index in range [0,numStates-1]
   * @return velocity: a 3d vector with velocity. By policy, we save it in
   * the cartesian basis and in atomic rydberg units.
   */
  Eigen::Vector3d getGroupVelocity(const long &stateIndex);

  /** Returns the group velocity of a quasiparticle from its Bloch index.
   * Used for accessing the bandstructure in the BTE.
   * @param stateIndex: a StateIndex(is) object where 'is' is an integer index
   * in the range [0,numStates[
   * @return velocity: a 3d vector with velocity. By policy, we save it in
   * the cartesian basis and in atomic rydberg units.
   */
  Eigen::Vector3d getGroupVelocity(StateIndex &is);

  /** Returns the group velocity of a quasiparticle for all bands at a
   * specified wavevector index.
   * Used for accessing the bandstructure in the BTE.
   * @param wavevectorIndex: a WavevectorIndex(ik) object where 'ik' is an
   * integer index running over the wavevectors in range [0,numPoints-1]
   * @return velocity: a matrix(numBands,3) with the group velocity, in
   * cartesian basis and in atomic rydberg units.
   */
  Eigen::MatrixXd getGroupVelocities(WavevectorIndex &ik);

  /** Returns the velocity operator (including off-diagonal matrix elements)
   * of the quasiparticles at the specified wavevector index.
   * Used for accessing the bandstructure in the BTE.
   * @param wavevectorIndex: a WavevectorIndex(ik) object where 'ik' is an
   * integer index running over the wavevectors in range [0,numPoints-1]
   * @return velocity: a tensor (numBands,numBands,3) with the velocity
   * operator matrix elements, in cartesian basis and in atomic rydberg units.
   */
  Eigen::Tensor<std::complex<double>, 3> getVelocities(WavevectorIndex &ik);

  /** Obtain the eigenvectors of the quasiparticles at a specified wavevector.
   * @param wavevectorIndex: a WavevectorIndex(ik) object where ik is the
   * integer wavevector index running over [0,numPoints-1].
   * @return eigenvectors: a complex matrix(numBands,numBands) where numBands
   * is the number of Bloch states present at the specified wavevector.
   * Eigenvectors are ordered along columns.
   * Note that all bandstructure interpolators may give eigenvectors.
   */
  Eigen::MatrixXcd getEigenvectors(WavevectorIndex &ik);

  /** Obtain the eigenvectors of the quasiparticles at a specified wavevector.
   * It's only meaningful for the phonon bandstructure, where eigenvectors
   * are more naturally represented in this shape!
   * @param wavevectorIndex: a WavevectorIndex(ik) object where ik is the
   * integer wavevector index running over [0,numPoints-1].
   * @return eigenvectors: a complex tensor (3,numAtoms,numBands) where
   * numBands is the number of Bloch states present at the specified
   * wavevector, numAtoms is the number of atoms in the crystal unit cell and
   * 3 is a cartesian directions.
   */
  Eigen::Tensor<std::complex<double>, 3> getPhEigenvectors(WavevectorIndex &ik);

  /** Returns the energy of a quasiparticle from its Bloch index
   * Used for accessing the bandstructure in the BTE.
   * @param stateIndex: an integer index in range [0,numStates[
   * @return wavevector: a 3d vector with the wavevector in cartesian
   * coordinates in units of Bohr^-1.
   */
  Eigen::Vector3d getWavevector(const long &stateIndex);

  /** Returns the energy of a quasiparticle from its Bloch index.
   * @param stateIndex: a StateIndex(is) object where 'is' is an integer
   * index in range [0,numStates-1].
   * @return wavevector: a 3d vector with the wavevector in cartesian
   * coordinates in units of Bohr^-1.
   */
  Eigen::Vector3d getWavevector(StateIndex &is);

  /** Returns the energy of a quasiparticle from its Bloch index.
   * @param wavevectorIndex: a WavevectorIndex(ik) object where 'ik' is an
   * integer index in range [0,numPoints-1].
   * @return wavevector: a 3d vector with the wavevector in cartesian
   * coordinates in units of Bohr^-1.
   */
  Eigen::Vector3d getWavevector(WavevectorIndex &ik);

  /** Returns the weight of a quasiparticle from its Bloch index, to be used
   * when integrating the Brillouin zone.
   * @param stateIndex: an integer index in range [0,numStates-1].
   * @return weight: a double value normalized such that
   * \f$\sum_{ik} weight(ik) = 1\f$.
   */
  double getWeight(const long &stateIndex);

  /** Returns the weight of a quasiparticle from its Bloch index, to be used
   * when integrating the Brillouin zone.
   * @param stateIndex: a StateIndex(is) object where 'is' is an integer
   * index in range [0,numStates-1].
   * @return weight: a double value normalized such that
   * \f$\sum_{ik} weight(ik) = 1\f$.
   */
  double getWeight(StateIndex &is);

  /** Returns the weight of a quasiparticle from its wavevector index, to be
   * used when integrating the Brillouin zone.
   * @param wavevectorIndex: a WavevectorIndex(ik) object where 'ik' is an
   * integer index in range [0,numPoints-1].
   * @return weight: a double value normalized such that
   * \f$\sum_{ik} weight(ik) = 1\f$.
   */
  double getWeight(WavevectorIndex &ik);

  /** Method to save quasiparticle energies inside FullBandStructure().
   * @param point: a point object, which also provides the wavevector index.
   * @param energies: a vector of size (numBands) with the quasiparticle
   * energies
   */
  void setEnergies(Point &point, Eigen::VectorXd &energies_);

  /** Method to save quasiparticle eigenvectors inside FullBandStructure().
   * Note that in this case, eigenvectors are passed as a matrix, which is
   * the case e.g. for the Wannier interpolation, where the eigenvectors
   * represent the unitary transformation matrix U for Wannier localization.
   * @param point: a Point object with the coordinates of the wavevector,
   * which should come from the same Point class stored in FullBandStructure
   * @param eigenvectors: a complex matrix of size (numBands,numBands)
   */
  void setEigenvectors(Point &point, Eigen::MatrixXcd &eigenvectors_);

  /** Saves in the class the velocities computed at a particular point.
   * @param point: a Point object representing the wavevector where these
   * velocities have been computed.
   * @param velocities: a rank-3 tensor of size (numBands,numBands,3)
   * containing the matrix elements of the velocity operator. Diagonal
   * elements are the quasiparticle group velocities.
   */
  void setVelocities(Point &point,
                     Eigen::Tensor<std::complex<double>, 3> &velocities_);

  /** Method to save quasiparticle eigenvectors inside FullBandStructure().
   * @param point: a vector of 3 crystal coordinates. The method will look
   * for the wavevector index.
   * @param energies: a vector of size (numBands) with the quasiparticle
   * energies
   */
  void setEnergies(Eigen::Vector3d &point, Eigen::VectorXd &energies_);

  /** Returns all electronic energies for all wavevectors at fixed band index
   * Used by the Fourier interpolation of the band structure.
   * @param bandIndex: index in [0,numBands[ for the quasiparticle band
   * @return energies: a vector of size numPoints with the QP energies.
   * Energies are ordered as the underlying wavevectors in Points.
   */
  Eigen::VectorXd getBandEnergies(long &bandIndex);

 protected:
  // stores the quasiparticle kind
  Particle particle;

  // link to the underlying mesh of points
  Points &points;  // these may be FullPoints or PathPoints

  // if the bands are distributed or not
  bool isDistributed = false;

  // matrices storing the raw data
  Matrix<double> energies;       // size(bands,points)
  Matrix<std::complex<double>> velocities;    // size(3*bands^2,points)
  Matrix<std::complex<double>> eigenvectors;  // size(bands^2,points)

  // auxiliary variables
  int numBands = 0;
  int numAtoms = 0;
  long numPoints = 0;
  long numLocalPoints = 0;

  bool hasEigenvectors = false;
  bool hasVelocities = false;

  // method to find the index of the kpoint, from its crystal coordinates
  long getIndex(Eigen::Vector3d &pointCoords);

  // ActiveBandStructure, which restricts the FullBandStructure to a subset
  // of wavevectors, needs low-level access to the raw data (for now)
  friend class ActiveBandStructure;
};

#endif
