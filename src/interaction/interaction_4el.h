#ifndef EL4_INTERACTION_H
#define EL4_INTERACTION_H

#include <complex>

#include "constants.h"
#include "crystal.h"
#include "eigen.h"
#include "phonon_h0.h"
#include "points.h"
#include "utilities.h"
#include "context.h"
#include <Kokkos_Core.hpp>

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
class Interaction4El {
  Crystal &crystal;

  int numWannier, numElBravaisVectors;

  std::vector<Eigen::Tensor<double, 4>> cacheCoupling1;
  std::vector<Eigen::Tensor<double, 4>> cacheCoupling2;
  std::vector<Eigen::Tensor<double, 4>> cacheCoupling3;

  // Kokkos View types
  using ComplexView1D = Kokkos::View<Kokkos::complex<double> *, Kokkos::LayoutRight>;
  using ComplexView2D = Kokkos::View<Kokkos::complex<double> **, Kokkos::LayoutRight>;
  using ComplexView3D = Kokkos::View<Kokkos::complex<double> ***, Kokkos::LayoutRight>;
  using ComplexView4D = Kokkos::View<Kokkos::complex<double> ****, Kokkos::LayoutRight>;
  using ComplexView5D = Kokkos::View<Kokkos::complex<double> *****, Kokkos::LayoutRight>;
  using ComplexView6D = Kokkos::View<Kokkos::complex<double> ******, Kokkos::LayoutRight>;
  using ComplexView7D = Kokkos::View<Kokkos::complex<double> *******, Kokkos::LayoutRight>;
  using IntView1D = Kokkos::View<int *, Kokkos::LayoutRight>;
  using DoubleView1D = Kokkos::View<double *, Kokkos::LayoutRight>;
  using DoubleView2D = Kokkos::View<double **, Kokkos::LayoutRight>;
  using DoubleView4D = Kokkos::View<double ****, Kokkos::LayoutRight>;
  using DoubleView5D = Kokkos::View<double *****, Kokkos::LayoutRight>;

  using HostComplexView1D = Kokkos::View<Kokkos::complex<double>*, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using HostComplexView2D = Kokkos::View<Kokkos::complex<double>**, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using HostComplexView3D = Kokkos::View<Kokkos::complex<double>***, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using HostComplexView4D = Kokkos::View<Kokkos::complex<double>****, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using HostComplexView5D = Kokkos::View<Kokkos::complex<double>*****, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using HostComplexView6D = Kokkos::View<Kokkos::complex<double>******, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using HostComplexView7D = Kokkos::View<Kokkos::complex<double>*******, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using HostDoubleView1D = Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using HostDoubleView2D = Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  // Kokkos Range types
  using Range2D = Kokkos::MDRangePolicy<Kokkos::Rank<2,Kokkos::Iterate::Right,Kokkos::Iterate::Right>>;
  using Range3D = Kokkos::MDRangePolicy<Kokkos::Rank<3,Kokkos::Iterate::Right,Kokkos::Iterate::Right>>;
  using Range4D = Kokkos::MDRangePolicy<Kokkos::Rank<4,Kokkos::Iterate::Right,Kokkos::Iterate::Right>>;
  using Range5D = Kokkos::MDRangePolicy<Kokkos::Rank<5,Kokkos::Iterate::Right,Kokkos::Iterate::Right>>;
  using Range6D = Kokkos::MDRangePolicy<Kokkos::Rank<6,Kokkos::Iterate::Right,Kokkos::Iterate::Right>>;
  using Range7D = Kokkos::MDRangePolicy<Kokkos::Rank<7,Kokkos::Iterate::Right,Kokkos::Iterate::Right>>;

  ComplexView6D elPhCached1a, elPhCached1b;
  ComplexView5D elPhCached2a, elPhCached2b, elPhCached2c;
  ComplexView7D couplingWannier_d;
  DoubleView2D elBravaisVectors_d;
  DoubleView1D elBravaisVectorsDegeneracies_d;

public:

  /** Base constructor
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
  Interaction4El(
      Crystal &crystal_,
      const Eigen::Tensor<std::complex<double>, 7> &couplingWannier_,
      const Eigen::MatrixXd &elBravaisVectors_,
      const Eigen::VectorXd &elBravaisVectorsDegeneracies_);

  /** Default constructor */
  Interaction4El(Crystal &crystal_);

  /** Copy constructor
   */
  Interaction4El(const Interaction4El &that);

  /** Assignment operator
   */
  Interaction4El &operator=(const Interaction4El &that);

  /** Computes the values of the el-ph coupling strength for transitions of
   * type k1,q3 -> k2, where k1 is one fixed wavevector, and k2,q3 are
   * wavevectors running in lists of wavevectors.
   * It is assumed that the relation (k2 = k1 + q3) holds.
   *
   * The call to this function must be preceded by a call to cacheElPh(),
   * which does a precomputation at fixed value of k1.
   * If not, results will be wrong.
   * Hence, structure a code calling this functions as:
   * for k1:
   *   cacheElPh(k1)
   *   for k2:
   *     k3 = k2 - k1
   *     calcCouplingSquared(k2,k3)
   *
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
      const std::vector<Eigen::MatrixXcd> &eigvecs3,
      const std::vector<Eigen::MatrixXcd> &eigvecs4,
      const std::vector<Eigen::Vector3d> &k3Cs,
      const std::vector<Eigen::Vector3d> &k4Cs);

  /** Computes a partial Fourier transform over the k1/R_el variables.
   * @param k1C: values of the k1 cartesian coordinates over which the Fourier
   * transform is computed.
   * @param eigvec1: Wannier rotation matrix U at point k1.
   */
  void cache1stEl(const Eigen::MatrixXcd &eigvec1, const Eigen::Vector3d &k1C);
  void cache2ndEl(const Eigen::MatrixXcd &eigvec2, const Eigen::Vector3d &k2C);

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
  std::tuple<Eigen::Tensor<double, 4>,Eigen::Tensor<double, 4>,
             Eigen::Tensor<double, 4>> getCouplingSquared(const int &ik3);

  /** Static method to initialize the class by parsing a file.
   * @param fileName: name of the file containing the coupling matrix elements
   * in real space, and the information on the lattice vectors and degeneracies.
   * @param crystal: object describing the crystal unit cell.
   * @return intElPh: an instance of InteractionElPh.
   */
  static Interaction4El parse(Context &context, Crystal &crystal);
};

#endif
