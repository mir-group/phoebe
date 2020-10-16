#ifndef ELPHINTERACTION_H
#define ELPHINTERACTION_H

#include <complex>

#include "constants.h"
#include "crystal.h"
#include "eigen.h"
#include "phonon_h0.h"
#include "points.h"
#include "utilities.h"
#include "context.h"

/** Class to handle the coupling between electron and phonons.
 * Currently implements the calculation of the diagram for the interaction
 * k+q -> k'.
 * Use the static method to initialize an instance of this class.
 * Then, use calc + get to compute and retrieve the values of the
 * electron-phonon interaction strength in Bloch space |g|^2.
 *
 * This class starts from the interaction matrix elements in real space Wannier
 * representation and mostly does:
 * 1) a double Fourier transform on phonon and electron coordinates;
 * 2) multiply by phonon/electron eigenvectors/rotation matrices;
 * 3) for polar materials, adds the long-range Frohlich interaction.
 */
// TODO: add flag to let user decide whether to use or not polar corrections
class InteractionElPhWan {
private:
  Crystal &crystal;
  PhononH0 *phononH0 = nullptr;

  Eigen::Tensor<std::complex<double>, 5> couplingWannier;
  // numElBands,numElBands,numPhBands,numPhBravaisVectors,numElBravaisVectors);

  Eigen::MatrixXd elBravaisVectors;
  Eigen::VectorXd elBravaisVectorsWeights;
  Eigen::MatrixXd phBravaisVectors;
  Eigen::VectorXd phBravaisVectorsWeights;

  int numPhBands, numElBands, numElBravaisVectors, numPhBravaisVectors;

  std::vector<Eigen::Tensor<double, 3>> cacheCoupling;

  Eigen::Tensor<std::complex<double>, 4> elPhCached;
  Eigen::Vector3d cachedK1;
  bool usePolarCorrection = false;

  /** Add polar correction to the electron-phonon coupling.
   * @param q3: phonon wavevector, in cartesian coordinates
   * @param ev1: eigenvector (rotation matrix U) at k
   * @param ev2: eigenvector (rotation matrix U) at k'
   * @param ev3: phonon eigenvector at q = k'-k
   * @return g^L: the long-range (Frohlich) component of the el-ph interaction,
   * as a tensor of shape (nb1,nb2,numPhBands)
   */
  Eigen::Tensor<std::complex<double>, 3>
  getPolarCorrection(const Eigen::Vector3d &q3, const Eigen::MatrixXcd &ev1,
                     const Eigen::MatrixXcd &ev2, const Eigen::MatrixXcd &ev3);

public:


  /** Default constructor
   * @param crystal_: object describing the crystal unit cell.
   * @param couplingWannier_: matrix elements of the electron phonon
   * interaction. A tensor of shape (iw1,iw2,imode,rPh,rEl), where iw1 iw2 are
   * indices on Wannier functions, imode is a phonon mode index in real space,
   * rPh is an index on phonon Bravais lattice vectors, and rEl is an index on
   * electronic Bravais Lattice vectors. Built such that the iw2 Wannier
   * functions are set in the origin (k2 doesn't contribute to the Fourier
   * transform).
   * @param elBravaisVectors_: list of Bravais lattice vectors used in the
   * electronic Fourier transform of the coupling.
   * @param elBravaisVectorsWeights_: weights (degeneracies) of the
   * lattice vectors used in the electronic Fourier transform of the coupling.
   * @param phBravaisVectors_: list of Bravais lattice vectors used in the
   * phonon Fourier transform of the coupling.
   * @param phBravaisVectorsWeights_: weights (degeneracies) of the
   * lattice vectors used in the phonon Fourier transform of the coupling.
   * @param phononH0_: pointer to the phonon dynamical matrix object. Used for
   * adding the polar interaction.
   */
  InteractionElPhWan(
      Crystal &crystal_,
      const Eigen::Tensor<std::complex<double>, 5> &couplingWannier_,
      const Eigen::MatrixXd &elBravaisVectors_,
      const Eigen::VectorXd &elBravaisVectorsWeights_,
      const Eigen::MatrixXd &phBravaisVectors_,
      const Eigen::VectorXd &phBravaisVectorsWeights_,
      PhononH0 *phononH0_ = nullptr);

  /** Copy constructor
   */
  InteractionElPhWan(const InteractionElPhWan &that);

  /** Assignment operator
   */
  InteractionElPhWan &operator=(const InteractionElPhWan &that);

  /** Computes the values of the el-ph coupling strength for transitions of
   * type k1,q3 -> k2, where k1 is one fixed wavevector, and k2,q3 are
   * wavevectors running in lists of wavevectors.
   * It is assumed that the relation (k2 = k1 + q3) holds.
   * Note: this method must be used in conjunction with getCouplingSquared,
   * which is used to return the values computed here.
   * @param el1Eigenvec: electron eigenvector matrix U_{mb}(k1), where U is
   * obtained by diagonalizing the Wannier Hamiltonian.
   * @param el2Eigenvecs: vector of electron eigenvectors matrix U_{mb}(k2),
   * where U is obtained by diagonalizing the Wannier Hamiltonian at a bunch
   * of k2 wavevectors.
   * @param phEigvecs: phonon eigenvectors, in matrix form, for a bunch of
   * wavevectors q3
   * @param k1: value of first wavevector.
   * @param k2s: list of k2 wavevectors.
   * @param q3s: list of phonon wavevectors.
   */
  void calcCouplingSquared(
      const Eigen::MatrixXcd &eigvec1,
      const std::vector<Eigen::MatrixXcd> &eigvecs2,
      const std::vector<Eigen::MatrixXcd> &eigvecs3,
      const Eigen::Vector3d &k1C,
      const std::vector<Eigen::Vector3d> &k2Cs,
      const std::vector<Eigen::Vector3d> &q3Cs);

  /** Get the coupling for the values of the wavevectors triplet (k1,k2,q3),
   * where k1 is the wavevector used at calcCoupling Squared(),
   * k2 (at index ik2) is the wavevector of the scattered electron in the
   * final state, and q3 = k2 - k1 is the phonon wavevector.
   * Note: this method only works AFTER calcCouplingSquared has been called.
   * @param ik2: index of the 2nd wavevector, aligned with the list of
   * wavevectors passed to calcCouplingSquared().
   * @return g2: a tensor of shape (nb1,nb2,numPhBands=3*numAtoms) with the
   * values of the coupling squared |g(ik1,ik2,iq3)|^2 for the el-ph transition
   * k1,q3 -> k2
   */
  Eigen::Tensor<double, 3> getCouplingSquared(const int &ik2);

  /** Static method to initialize the class by parsing a file.
   * @param fileName: name of the file containing the coupling matrix elements
   * in real space, and the information on the lattice vectors and degeneracies.
   * @param crystal: object describing the crystal unit cell.
   * @return intElPh: an instance of InteractionElPh.
   */
  static InteractionElPhWan parse(Context &context, Crystal &crystal,
                                  PhononH0 *phononH0_ = nullptr);

  static Eigen::Tensor<std::complex<double>, 3> getPolarCorrectionStatic(
      const Eigen::Vector3d &q3, const Eigen::MatrixXcd &ev1,
      const Eigen::MatrixXcd &ev2, const Eigen::MatrixXcd &ev3,
      const double &volume, const Eigen::Matrix3d &reciprocalUnitCell,
      const Eigen::Matrix3d &epsilon,
      const Eigen::Tensor<double, 3> &bornCharges,
      const Eigen::MatrixXd &atomicPositions,
      const Eigen::Vector3i &qCoarseMesh);
};

#endif
