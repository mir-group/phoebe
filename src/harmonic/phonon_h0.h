#ifndef PHONON_H0_H
#define PHONON_H0_H

#include <cmath>

#include "bandstructure.h"
#include "constants.h"
#include "crystal.h"
#include "eigen.h"
#include "harmonic.h"
#include "particle.h"
#include "points.h"

/** class that computes phonon energies, velocities and eigenvectors.
 * First, it contains the force constants, i.e. the second derivative of the
 * total energy w.r.t. ionic displacements. Additionally, it has all the
 * infrastructure to Fourier transform this matrix, and diagonalize it to get
 * the harmonic phonon properties.
 */
class PhononH0 : public HarmonicHamiltonian {
 public:
  /** Constructor, which stores all input data.
   * @param crystal: the object with the information on the crystal structure
   * @param dielectricMatrix: 3x3 matrix with the dielectric matrix
   * @param bornCharges: real tensor of size (numAtoms,3,3) with the Born
   * effective charges
   * @param forceConstants: a tensor of doubles with the force constants
   * size is (meshX, meshY, meshZ, 3, 3, numAtoms, numAtoms)
   */
  PhononH0(Crystal &crystal, const Eigen::Matrix3d &dielectricMatrix_,
           const Eigen::Tensor<double, 3> &bornCharges_,
           const Eigen::Tensor<double, 7> &forceConstants_,
           const std::string &sumRule);

  /** Copy constructor
   */
  PhononH0(const PhononH0 &that);

  /** Copy assignment operator
   */
  PhononH0 &operator=(const PhononH0 &that);

  /** Returns the number of phonon bands for the crystal in consideration.
   */
  int getNumBands() override;

  /** Returns the underlying phonon-boson particle.
   */
  Particle getParticle() override;

  /** get the phonon energies (in Ry) at a single q-point.
   * @param q: a point object with the wavevector. Must know the cartesian
   * coordinates of the wavevector.
   * @return tuple(energies, eigenvectors): the energies are a double vector
   * of size (numBands=3*numAtoms). Eigenvectors are a complex tensor of
   * size (3,numAtoms,numBands). The eigenvector is rescaled by the
   * sqrt(masses) (masses in rydberg)
   */
  std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> diagonalize(Point &point) override;

  /** Equivalent to diagonalize() computes phonon eigenValues/Vectors given the
   * wavevector, but the wavevector is passed by coordinates.
   * @param q: a 3d eigen vector with the cartesian coordinates of the
   * phonon wavevector.
   * @param withMassScaling: if true, rescales the eigenvectors by the
   * mass z -> z/sqrt(m)
   * @return eigenvalues: all values of phonon energies for this point.
   * @return eigenvectors: the phonon eigenvectors, in matrix form, for this
   * point.
   */
  std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> diagonalizeFromCoordinates(
      Eigen::Vector3d &q, const bool &withMassScaling);
  std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> diagonalizeFromCoordinates(
      Eigen::Vector3d &q) override;

  /** get the phonon velocities (in atomic units) at a single q-point.
   * @param q: a Point object with the wavevector coordinates.
   * @return velocity(numBands,numBands,3): values of the velocity operator
   * for this state, in atomic units.
   */
  Eigen::Tensor<std::complex<double>, 3> diagonalizeVelocity(Point &point) override;
  Eigen::Tensor<std::complex<double>, 3> diagonalizeVelocityFromCoordinates(
      Eigen::Vector3d &coordinates) override;

  /** This method constructs a phonon band structure.
   * @param points: the object with the list/mesh of wavevectors
   * @param withVelocities: if true, compute the phonon velocity operator.
   * @param withEigenvectors: if true, stores the phonon eigenvectors.
   * @return FullBandStructure: the band structure object containing the
   * complete phonon band structure.
   */
  FullBandStructure populate(Points &points, bool &withVelocities,
                             bool &withEigenvectors, bool isDistributed = false) override;

  std::tuple<std::vector<Eigen::VectorXd>, std::vector<Eigen::MatrixXcd>,
             std::vector<Eigen::Tensor<std::complex<double>,3>>> populate(
      const std::vector<Eigen::Vector3d>& cartesianCoordinates,
      const bool& withVelocities=false, const bool &withMassScaling=true);

  /** Returns the size of the q-point coarse grid on which the force constants
   * have been computed.
   * @return qCoarseGrid: an Eigen vector of 3 integers.
   */
  Eigen::Vector3i getCoarseGrid();

  /** Utility to convert the indices of atom basis and polarization into
   *  the index that has to be used with the eigenvector in matrix form.
   * @param iAt: atomic basis index.
   * @param iPol: polarization index (0,1,2).
   * @return k: index to be used in the phonon eigenvector.
   */
  int getIndexEigenvector(const int &iAt, const int &iPol) const;

  /** same as getIndexEigenvector, but as a static member
   * @param nAtoms: the number of atoms in the unit cell
   */
  static int getIndexEigenvector(const int &iAt, const int &iPol, const int &nAtoms);

  /** Get the static dielectric constant matrix.
   * @return dielectricMatrix: a 3x3 eigen matrix.
   */
  Eigen::Matrix3d getDielectricMatrix();

  /** Get the Born effective charges for the ions in the unit cell
   * @return bornCharges: a real tensor of shape (numAtoms,3,3) with the charges
   */
  Eigen::Tensor<double, 3> getBornCharges();

protected:
  /** Impose the acoustic sum rule on force constants and Born charges
   * @param sumRule: name of the sum rule to be used
   * Currently supported values are akin to those from Quantum ESPRESSO
   * i.e. "simple" (for a rescaling of the diagonal elements) or "crystal"
   * (to find the closest matrix which satisfies the sum rule)
   */
  void setAcousticSumRule(const std::string &sumRule);

  void reorderDynamicalMatrix();

  Particle particle;

  bool hasDielectric = false;
  int numAtoms;
  int numBands;
  Eigen::MatrixXd directUnitCell;
  Eigen::MatrixXd reciprocalUnitCell;
  double volumeUnitCell;
  Eigen::MatrixXi atomicSpecies;
  Eigen::VectorXd speciesMasses;
  Eigen::MatrixXd atomicPositions;
  Eigen::Matrix3d dielectricMatrix;
  Eigen::Tensor<double, 3> bornCharges;
  Eigen::Vector3i qCoarseGrid;
  Eigen::Tensor<double, 7> forceConstants;
  Eigen::Tensor<double, 5> wsCache;
  int nr1Big, nr2Big, nr3Big;

  int numBravaisVectors = 0;
  Eigen::MatrixXd bravaisVectors;
  Eigen::VectorXd weights;
  Eigen::Tensor<double,3> mat2R;

  Eigen::MatrixXd gVectors;
  Eigen::Tensor<double,3> longRangeCorrection1;
  double gMax = 14.; // cutoff for ewald summation

  // private methods, used to diagonalize the Dyn matrix

  /** In wsInit, starting from the primitive crystal unit cell, we build the
   * list of bravais lattice vectors used for the phonon Fourier transform.
   */
  void wsInit(const Eigen::MatrixXd &unitCell);

  /** wsWeight computes the `weights`, i.e. the number of symmetry-equivalent
   * Bravais lattice vectors, that are used in the phonon Fourier transform.
   */
  static double wsWeight(const Eigen::VectorXd &r, const Eigen::MatrixXd &rws);

  /** Auxiliary methods for sum rule on Born charges
   */
  void sp_zeu(Eigen::Tensor<double, 3> &zeu_u, Eigen::Tensor<double, 3> &zeu_v,
              double &scalar) const;

  /** This method does the Fourier transform of the force constants to construct
   * the dynamical matrix, and returns the phonon frequencies and eigenvectors.
   *
   * @param cartesianCoordinates: vector of q-point coordinates.
   * @param withMassScaling: controls the normalization of the phonon
   * eigenvectors. If false, phonon eigenvectors are ortho-normalized to 1. If
   * true, the normalization of the eigenvector is proportional to the masses.
   * @return tuple with phonon frequencies and eigenvectors. Remember that the
   * square of the frequencies are the eigenvalues. Also, if phonon frequencies
   * should be imaginary, they are returned with negative values.
   */
  std::tuple<std::vector<Eigen::VectorXd>,
             std::vector<Eigen::MatrixXcd>> internalPopulate(
      const std::vector<Eigen::Vector3d>& cartesianCoordinates,
      const bool& withMassScaling);

  /** Adds the long range correction to the dynamical matrix due to dipole-ion
   * interaction. Called by internalPopulate()
   */
  void addLongRangeTerm(std::vector<Eigen::MatrixXcd> &dyns,
                        const std::vector<Eigen::Vector3d>& qs);

};

#endif
