#ifndef ACT_BANDSTRUCTURE_H
#define ACT_BANDSTRUCTURE_H

#include "active_points.h"
#include "bandstructure.h"
#include "constants.h"
#include "full_points.h"
#include "harmonic.h"
#include "particle.h"
#include "points.h"
#include "statistics_sweep.h"
#include "window.h"

/** Class container of the quasiparticle bandstructure, i.e. energies,
 * velocities, eigenvectors and wavevectors.
 * In contrast to FullBandStructure, which stores this information for all
 * wavevectors and bands, ActiveBandStructure stores this data for a specified
 * subset of wavevectors, and a subset of bands which may change from one
 * wavevector to another.
 * Should be initialized through the static method provided below, which makes
 * use of the Window that is selected by user input.
 * The subset of wavevectors is specified using an ActivePoints class.
 */
class ActiveBandStructure : public BaseBandStructure {
public:
  /** Almost empty constructor, to be used internally.
   */
  ActiveBandStructure(Particle &particle_, ActivePoints &activePoints);

  /** Copy constructor
   */
  ActiveBandStructure(const ActiveBandStructure &that);

  /** Assignment operator
   */
  ActiveBandStructure &operator=(const ActiveBandStructure &that);

  /** Default constructor.
   * Populate the bandstructure for the subset of wavevectors specified by
   * activePoints, but saving all bands for such wavevector.
   * Using this constructor, we populate the bandstructure on the wavevectors
   * specified by activePoints.
   * Currently used for the caching of the 3rd wavevectors info in the
   * 3-phonon scattering matrix.
   */
  ActiveBandStructure(const ActivePoints &activePoints_,
                      HarmonicHamiltonian *h0, const bool &withEigenvectors,
                      const bool &withVelocities);

  /** Get the Particle object associated with this class
   * @return particle: a Particle object, describing e.g. whether this
   * is a phonon or electron bandStructure
   */
  Particle getParticle();

  /** Returns the wavevectors on which the bandstructure is computed.
   * @return Points: the object representing the Brillouin zone wavevectors.
   * Should be an ActivePoints class.
   */
  Points getPoints();

  /** Returns a wavevector, given a wavevector index.
   * The wavevector index runs from 0 to numPoints-1, where numPoints is the
   * number of active wavevectors to be used in the calculations.
   */
  Point getPoint(const long &pointIndex);

  /** Returns the total number of k/q-points.
   * @param useFullGrid: default = false. If true, returns the number of
   * points of the full monkhorst-pack grid of wavevectors, otherwise,
   * returns the number of wavevectors stored in the BandStructure.
   * In particular for this class, the number of wavevectors stored in the list
   * is obtained with false, if true, and if activePoints.parentPoints is a
   * FullPoints obhject, then we will return the number of points in the
   * original full grid of points (useful when normalizing equations).
   * @return numPoints: the total number of wavevectors of the bandStructure.
   */
  long getNumPoints(const bool &useFullGrid = false);

  /** Returns the number of bands.
   * If the number of bands is not constant, it raises an error.
   * In that case, it's best to check the size of the objects returned by
   * the get*() methods or use getNumBands(ik).
   * @return numBandsFull: the total number of bands of the bandStructure
   *    (only returned in the case where the bands have not been filtered.)
   */
  long getNumBands();

  /** Returns the number of bands at a given wavevector.
   * @return numBands: the number of bands at the requested ik.
   */
  long getNumBands(WavevectorIndex &ik);

  /** Checks whether the bandStructure has been built discarding some Bloch
   * states from those available.
   * @return windowMethod: one of the values of Window::filterMethod. 0 for
   * no filter.
   */
  long hasWindow();

  /** Returns if this bandstructure is distributed. In the case of
   * activeBandstructure, currently always returns false.
   */
  bool getIsDistributed();

  /** Builds a Bloch state index, which runs on both wavevector index and
   * band index.
   * @param wavevectorIndex: strong-typed index on wavevector index
   * @param bandIndex: strong-typed index on band index. Note that bandIndex
   * has different range for every wavevector
   * @return stateIndex: integer from 0 to numStates-1
   */
  long getIndex(const WavevectorIndex &ik, const BandIndex &ib);

  /** Given a Bloch state index, finds the corresponding wavevector and band
   * index.
   * @param stateIndex: integer from 0 to numStates-1
   * @return WavevectorIndex: a strong-typed index on wavevector
   * @return BandIndex: a strong-typed index on bands
   */
  std::tuple<WavevectorIndex, BandIndex> getIndex(const long &is);

  /** Given a Bloch state index, finds the corresponding wavevector and band
   * index.
   * @param stateIndex: a StateIndex(is) object where is is an integer from 0
   * to numStates-1
   * @return WavevectorIndex: strong-typed index on wavevector
   * @return BandIndex: strong-typed index on bands
   */
  std::tuple<WavevectorIndex, BandIndex> getIndex(StateIndex &is);

  /** Returns the total number of Bloch states, equal to numPoints*numBands.
   * @return numStates: the total number of Bloch states in the class.
   */
  long getNumStates();

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
   * @return velocity: a matrix(numActiveBands,3) with the group velocity, in
   * cartesian basis and in atomic rydberg units, where numActiveBands is the
   * number of active bands present at the specified wavevector.
   */
  Eigen::MatrixXd getGroupVelocities(WavevectorIndex &ik);

  /** Returns the velocity operator (including off-diagonal matrix elements)
   * of the quasiparticles at the specified wavevector index.
   * Used for accessing the bandstructure in the BTE.
   * @param wavevectorIndex: a WavevectorIndex(ik) object where 'ik' is an
   * integer index running over the wavevectors in range [0,numPoints-1]
   * @return velocity: a tensor (numActiveBands,numActiveBands,3) with the
   * velocity operator matrix elements, in cartesian basis and in atomic rydberg
   * units, where numActiveBands is the number of active bands present at the
   * specified wavevector.
   */
  Eigen::Tensor<std::complex<double>, 3> getVelocities(WavevectorIndex &ik);

  /** Obtain the eigenvectors of the quasiparticles at a specified wavevector.
   * @param wavevectorIndex: a WavevectorIndex(ik) object where ik is the
   * integer wavevector index running over [0,numPoints-1].
   * @return eigenvectors: a complex matrix(numBands,numActiveBands) where
   * numBands is the full (before applying the window) number of Bloch states
   * present at the specified wavevector and numActiveBands is the number of
   * filtered bands at the desired wavevector.
   * Eigenvectors are ordered along columns.
   * Note that all bandstructure interpolators may give eigenvectors.
   */
  Eigen::MatrixXcd getEigenvectors(WavevectorIndex &ik);

  /** Obtain the eigenvectors of the quasiparticles at a specified wavevector.
   * It's only meaningful for the phonon bandstructure, where eigenvectors
   * are more naturally represented in this shape!
   * @param wavevectorIndex: a WavevectorIndex(ik) object where ik is the
   * integer wavevector index running over [0,numPoints-1].
   * @return eigenvectors: a complex tensor (3,numAtoms,numActiveBands) where
   * numActiveBands is the number of Bloch states present at the specified
   * wavevector, numAtoms is the number of atoms in the crystal unit cell and
   * 3 are the cartesian directions.
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
   * @param stateIndex: a StateIndex(is) object where 'is' is an integer
   * index in range [0,numStates-1].
   * @return weight: a double value normalized such that the summation
   * \f$\sum_{ik} weight(ik) = 1\f$ if we were summing over the complete
   * (non-filtered) list of wavevectors (but we probably aren't).
   */
  double getWeight(StateIndex &is);

  /** Returns the weight of a quasiparticle from its wavevector index, to be
   * used when integrating the Brillouin zone.
   * @param wavevectorIndex: a WavevectorIndex(ik) object where 'ik' is an
   * integer index in range [0,numPoints-1].
   * @return weight: a double value normalized such that the summation
   * \f$\sum_{ik} weight(ik) = 1\f$ if we were summing over the complete
   * (non-filtered) list of wavevectors (but we probably aren't).
   */
  double getWeight(WavevectorIndex &ik);

  /** Method to save quasiparticle energies inside ActiveBandStructure().
   * @param point: a point object, which will be used to lookup the wavevector
   * index ik.
   * @param energies: a vector of size (numBands(ik)) with the quasiparticle
   * energies for that point.
   */
  void setEnergies(Point &point, Eigen::VectorXd &energies_);

  /** Method to save quasiparticle energies inside ActiveBandStructure().
   * @param point: a point object, which will be used to lookup the wavevector
   * index ik.
   * @param energies: a vector of size (numBands(ik)) with the quasiparticle
   * energies for that point.
   */
  void setEnergies(Point &point, std::vector<double> &energies_);

  /** Method to save quasiparticle eigenvectors inside ActiveBandStructure.
   * @param point: a point object, which will be used to lookup the wavevector
   * index ik.
   * @param eigenvectors: a matrix of size (numFullBands, numBands(ik)) with
   * the quasiparticle eigenvectors for this point.
   */
  void setEigenvectors(Point &point, Eigen::MatrixXcd &eigenvectors_);

  /** Method to save quasiparticle velocity operator in ActiveBandStructure.
   * @param point: a point object, which will be used to lookup the wavevector
   * index ik.
   * @param velocities: a tensor of size (numBands(ik),numBands(ik),3) with
   * the quasiparticle velocity matrix elements for this wavevector.
   */
  void setVelocities(Point &point,
                     Eigen::Tensor<std::complex<double>, 3> &velocities_);

  /** Preferred method to initialize the ActiveBandStructure class.
   * @param context: object with user input variables.
   * @param h0: hamiltonian object to be diagonalized.
   * @param points: initial unfiltered list of wavevectors, which will be
   * filtered into an ActivePoints object.
   * @param withEigenvectors: compute and store the eigenvectors
   * @param withVelocities: compute and store the velocity matrix elements
   * @param forceBuildAPP: forces activeBandStructure to be built
   * using the internal buildAsPostprocessing method, even if the input
   * H0 is for phonons.
   */
  static std::tuple<ActiveBandStructure, StatisticsSweep>
  builder(Context &context, HarmonicHamiltonian &h0, Points &points,
          const bool &withEigenvectors = true,
          const bool &withVelocities = true, const bool &forceBuildAPP = false);

  std::vector<Eigen::Matrix3d> getRotationsStar(WavevectorIndex &ikIndex);
  std::vector<Eigen::Matrix3d> getRotationsStar(StateIndex &isIndex);

  std::tuple<long, Eigen::Matrix3d> getRotationToIrreducible(
      const Eigen::Vector3d &x, const int &basis = Points::crystalCoords);

  BteIndex stateToBte(StateIndex &isIndex);
  StateIndex bteToState(BteIndex &ibteIndex);

  std::vector<long> irrStateIterator();
  std::vector<long> parallelIrrStateIterator();
  std::vector<long> irrPointsIterator();
  std::vector<long> parallelIrrPointsIterator();

  long getPointIndex(const Eigen::Vector3d &crystalCoords,
                     const bool &suppressError = false);
  std::vector<long> getReduciblesFromIrreducible(const long &ik);
 protected:
  // stores the quasiparticle kind
  Particle particle;

  ActivePoints activePoints;

  // note: we don't store a matrix: we are storing an object (Nk,Nb),
  // with a variable number of bands Nb per point
  std::vector<double> energies;
  std::vector<std::complex<double>> velocities;
  std::vector<std::complex<double>> eigenvectors;

  bool hasEigenvectors = false;
  long numStates = 0;
  long numIrrStates;
  long numIrrPoints;
  long numPoints;
  bool hasPoints();

  Eigen::VectorXi numBands;
  long numFullBands;
  long windowMethod;

  // index management
  // these are two auxiliary vectors to store indices
  Eigen::MatrixXi auxBloch2Comb;
  Eigen::VectorXi cumulativeKbOffset;
  Eigen::MatrixXi bteAuxBloch2Comb;
  Eigen::VectorXi bteCumulativeKbOffset;
  Eigen::VectorXi cumulativeKbbOffset;
  // this is the functionality to build the indices
  void buildIndeces(); // to be called after building the band structure
  // and these are the tools to convert indices

  // utilities to convert Bloch indices into internal indices
  long velBloch2Comb(const long &ik, const long &ib1, const long &ib2,
                     const long &i);
  long eigBloch2Comb(const long &ik, const long &ibFull, const long &ibRed);
  long bloch2Comb(const long &k, const long &b);
  std::tuple<long, long> comb2Bloch(const long &is);

  long bteBloch2Comb(const long &k, const long &b);
  std::tuple<long, long> bteComb2Bloch(const long &is);

  void buildOnTheFly(Window &window, Points &points, HarmonicHamiltonian &h0,
                     const bool &withEigenvectors = true,
                     const bool &withVelocities = true);

  StatisticsSweep buildAsPostprocessing(Context &context, Points &points,
                                        HarmonicHamiltonian &h0,
                                        const bool &withEigenvector = true,
                                        const bool &withVelocities = true);
};

#endif
