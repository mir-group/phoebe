#ifndef H0FOURIER_H
#define H0FOURIER_H

#include <string>

#include "bandstructure.h"
#include "harmonic.h"
#include "points.h"
#include "full_points.h"
#include "context.h"

/** Class for a Fourier-like interpolation of an electronic band structure.
 * Takes the information on the band structure computed on a uniform coarse
 * grid of k-points, and interpolates it to finer grids with a plane wave based
 * method.
 */
class ElectronH0Fourier : public HarmonicHamiltonian {
 public:
  const bool hasEigenvectors = false;

  /** Constructor of the Fourier interpolation
   * This class stores a copy of the electronic band structure on the coarse
   * grid.
   * @param crystal: the crystal used in the band structure calculation
   * @param coarseBandStructure: values of the electronic bands over a full
   * grid of kpoints, provided by an external (DFT) code.
   * @param cutoff: a parameter used to define the number of coefficients
   * in the plane-wave interpolation. It should be an integer >1.
   * This parameter controls that the interpolation is generated using a grid
   * of lattice vectors which is of size 2*Grid(i)*cutoff, where i is the
   * direction index and grid(i) is the size of the wavevector coarse grid.
   */
  ElectronH0Fourier(Crystal &crystal_, FullPoints coarsePoints_,
                    FullBandStructure coarseBandStructure_, double cutoff);

  /** Copy constructor
   */
  ElectronH0Fourier(const ElectronH0Fourier &that);

  /** Copy assignment
   */
  ElectronH0Fourier &operator=(const ElectronH0Fourier &that);

  /** Method to return that the underlying is that of an electronic Fermion.
   */
  Particle getParticle();

  /** Get the total number of bands available at ech wavevector.
   *
   */
  long getNumBands();

  /** get the electronic energies (in Ry) at a single k-point.
   * Energies don't have any reference value, and must be used in connection
   * with a chemical potential.
   * @param k: a point object with the wavevector. Must have the cartesian
   * coordinates of the wavevector.
   * @return tuple(energies, eigenvectors): the energies are a double vector
   * of size (numBands). Eigenvectors of size (numBands,numBands), but are
   * simply set to zero, since there is no diagonalization happening here.
   */
  std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> diagonalize(Point &point);

  std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> diagonalizeFromCoords(
      Eigen::Vector3d &wavevector);

  /** get the electron velocities (in atomic units) at a single k-point.
   * @param k: a Point object with the wavevector coordinates.
   * @return velocity(numBands,numBands,3): values of the velocity operator
   * for this state, in atomic units. Note that the off-diagonal matrix
   * elements are set to zero, because this kind of interpolation, at the
   * moment, doesn't have any information on the off-diagonal elements.
   */
  Eigen::Tensor<std::complex<double>, 3> diagonalizeVelocity(Point &point);
  Eigen::Tensor<std::complex<double>, 3> diagonalizeVelocityFromCoords(
      Eigen::Vector3d &coords);

  /** This method constructs an electron bandstructure.
   * @param points: the object with the list/mesh of wavevectors
   * @param withVelocities: if true, compute the electron velocity operator.
   * @param withEigenvectors: can only be false, as there are no eigenvectors
   * with this kind of interpolation.
   * @return FullBandStructure: the bandstructure object containing the
   * complete electronic band structure.
   */
  FullBandStructure populate(Points &fullPoints, bool &withVelocities,
                             bool &withEigenvectors, bool isDistributed=false);

  void trimBands(Context &context, const double &minEn, const double &maxEn);

 protected:
  Crystal &crystal;
  FullBandStructure coarseBandStructure;
  FullPoints coarsePoints;
  Particle particle;

  Eigen::MatrixXcd expansionCoefficients;

  long numBands;
  double cutoff;
  long numDataPoints;
  long numPositionVectors;
  double minDistance;
  Eigen::VectorXd positionDegeneracies;
  Eigen::MatrixXd positionVectors;
  Eigen::Vector3d refWavevector;

  void setPositionVectors();
  Eigen::VectorXcd getLagrangeMultipliers(Eigen::VectorXd energies);
  Eigen::VectorXcd getCoefficients(Eigen::VectorXd energies);
  std::complex<double> getStarFunction(Eigen::Vector3d &wavevector, long &iR);
  Eigen::Vector3cd getDerivativeStarFunction(Eigen::Vector3d &wavevector,
                                             long &iR);
  double getRoughnessFunction(Eigen::Vector3d position);
  const double coeff1 = 0.75;  // 3/4
  const double coeff2 = 0.75;
  double getEnergyFromCoords(Eigen::Vector3d &wavevector, long &bandIndex);
  Eigen::Vector3d getGroupVelocityFromCoords(Eigen::Vector3d &wavevector,
                                             long &bandIndex);
};

#endif
