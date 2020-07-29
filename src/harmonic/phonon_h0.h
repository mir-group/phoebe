#ifndef PHONONH0_H
#define PHONONH0_H

#include <math.h>

#include "bandstructure.h"
#include "constants.h"
#include "crystal.h"
#include "eigen.h"
#include "harmonic.h"
#include "particle.h"
#include "points.h"

/** class that computes phonon energies, velocities and eigenvectors.
 * FIrst, it contains the force constants, i.e. the second derivative of the
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
   * size is (meshx, meshy, meshz, 3, 3, numAtoms, numAtoms)
   */
  PhononH0(Crystal &crystal, const Eigen::MatrixXd &dielectricMatrix_,
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
  long getNumBands();

  /** Returns the underlying phonon-boson particle.
   */
  Particle getParticle();

  /** get the phonon energies (in Ry) at a single q-point.
   * @param q: a point object with the wavevector. Must know the cartesian
   * coordinates of the wavevector.
   * @return tuple(energies, eigenvectors): the energies are a double vector
   * of size (numBands=3*numAtoms). Eigenvectors are a complex tensor of
   * size (3,numAtoms,numBands). The eigenvector is rescaled by the
   * sqrt(masses) (masses in rydbergs)
   */
  std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> diagonalize(Point &point);

  /** Equivalent to diagonalize() computes phonon eigenvals/vecs given the
   * wavevector, but the wavevector is passed by coordinates.
   * @param q: a 3d eigen vector with the cartesian coordinates of the
   * phonon wavevector.
   * @param withMassScaling: if true, rescales the eigenvectors by the
   * mass z -> z/sqrt(m)
   * @return eigenvalues: all values of phonon energies for this point.
   * @return eigenvectors: the phonon eigenvectors, in matrix form, for this
   * point.
   */
  std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> diagonalizeFromCoords(
      Eigen::Vector3d &q, const bool &withMassScaling);
  std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> diagonalizeFromCoords(
      Eigen::Vector3d &q);

  /** get the phonon velocities (in atomic units) at a single q-point.
   * @param q: a Point object with the wavevector coordinates.
   * @return velocity(numBands,numBands,3): values of the velocity operator
   * for this state, in atomic units.
   */
  Eigen::Tensor<std::complex<double>, 3> diagonalizeVelocity(Point &point);
  Eigen::Tensor<std::complex<double>, 3> diagonalizeVelocityFromCoords(
      Eigen::Vector3d &coords);

  /** This method constructs a phonon bandstructure.
   * @param points: the object with the list/mesh of wavevectors
   * @param withVelocities: if true, compute the phonon velocity operator.
   * @param withEigenvectors: if true, stores the phonon eigenvectors.
   * @return FullBandStructure: the bandstructure object containing the
   * complete phonon band structure.
   */
  FullBandStructure populate(Points &points, bool &withVelocities,
                             bool &withEigenvectors);

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

  Eigen::Vector3i getCoarseGrid();
  // internal variables

  // these 3 variables might be used for extending future functionalities.
  // for the first tests, they can be left at these default values
  // in the future, we might expose them to the user input
  bool na_ifc = false;
  bool loto_2d = false;
  bool frozenPhonon = false;

  bool hasDielectric;
  long numAtoms;
  long numBands;
  Eigen::MatrixXd directUnitCell;
  Eigen::MatrixXd reciprocalUnitCell;
  double latticeParameter;
  double volumeUnitCell;
  Eigen::MatrixXi atomicSpecies;
  Eigen::VectorXd speciesMasses;
  Eigen::MatrixXd atomicPositions;
  Eigen::MatrixXd dielectricMatrix;
  Eigen::Tensor<double, 3> bornCharges;
  Eigen::Vector3i qCoarseGrid;
  Eigen::Tensor<double, 7> forceConstants;
  Eigen::Tensor<double, 5> wscache;
  int nr1Big, nr2Big, nr3Big;

  int numBravaisVectors;
  Eigen::MatrixXd bravaisVectors;
  Eigen::VectorXd weights;
  Eigen::Tensor<double,5> mat2R;

  // private methods, used to diagonalize the Dyn matrix
  void wsinit(const Eigen::MatrixXd &unitCell);
  double wsweight(const Eigen::VectorXd &r, const Eigen::MatrixXd &rws);
  void longRangeTerm(Eigen::Tensor<std::complex<double>, 4> &dyn,
                     const Eigen::VectorXd &q, const long sign);
  void nonAnaliticTerm(const Eigen::VectorXd &q,
                       Eigen::Tensor<std::complex<double>, 4> &dyn);
  void nonAnalIFC(const Eigen::VectorXd &q,
                  Eigen::Tensor<std::complex<double>, 4> &f_of_q);
  void shortRangeTerm(Eigen::Tensor<std::complex<double>, 4> &dyn,
                      const Eigen::VectorXd &q,
                      Eigen::Tensor<std::complex<double>, 4> &f_of_q);
  std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> dyndiag(
      Eigen::Tensor<std::complex<double>, 4> &dyn);

  // methods for sum rule on Born charges
  void sp_zeu(Eigen::Tensor<double, 3> &zeu_u, Eigen::Tensor<double, 3> &zeu_v,
              double &scal);
};

#endif
