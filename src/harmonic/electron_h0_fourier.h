#include <string>
#include "points.h"
#include "bandstructure.h"
#include "harmonic.h"

/** Class for a Fourier-like interpolation of an electronic band
 *
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
	 * in the plane-wave interpolation.
	 */
	ElectronH0Fourier(Crystal & crystal_, FullPoints coarsePoints_,
			FullBandStructure coarseBandStructure_, double cutoff);

	/** calculate the energy of a single electronic Bloch state.
	 * @param point: the Point object representing the wavevector.
	 * @param bandIndex: the number identifying the band from 0 to numBands-1
	 * @return energy: in Rydbergs
	 */
	double getEnergy(Point & point, long & bandIndex);

	/** get all electronic energies at a given wavevector
	 * @param point: a Point representing the desired wavevector.
	 * @return energies: a vector of size (numBands) of energies
	 */
	Eigen::VectorXd getEnergies(Point & point);

	/** compute the group velocity of a Bloch state.
	 * @param point: a Point representing the desired wavevector.
	 * @param bandIndex: the number identifying the band from 0 to numBands-1
	 * @return velocity: a 3d vector with the group velocity
	 */
	Eigen::Vector3d getGroupVelocity(Point & point, long & bandIndex);

	/** compute the group velocities at a given wavevector
	 * @param point: a Point representing the desired wavevector.
	 * @return velocities: a matrix of size (numBands,3) of group velocities
	 */
	Eigen::MatrixXd getGroupVelocities(Point & point);

	/** compute the band structure on a mesh of points
	 * Must provide only one of the two optional parameters
	 * @param (optional) fullPoints: k-point mesh to use for computing energies
	 * @param (optional) irreduciblePoints: k-point mesh for computing energies
	 * @return FullBandStructure: an FullBandStructure object with the values
	 * of the band structure computed over the input mesh of points.
	 *
	 */
	FullBandStructure populateBandStructure(FullPoints * fullpoints=nullptr,
			IrreduciblePoints * irreduciblePoints=nullptr);

	virtual std::tuple<Eigen::VectorXd, Eigen::Tensor<std::complex<double>,3>>
		diagonalize(Point & point);

	virtual Eigen::Tensor<std::complex<double>,3> diagonalizeVelocity(
				Point & point);
protected:
	FullPoints coarsePoints;
	FullBandStructure coarseBandStructure;

	Crystal & crystal;

	Eigen::MatrixXcd expansionCoefficients;

	long numBands;
	double cutoff;
	long numDataPoints;
	long numPositionVectors;
	double minDistance;
	std::vector<long> positionDegeneracies;
	std::vector<Eigen::Vector3d> positionVectors;
	void setPositionVectors();
	Eigen::VectorXcd getLagrangeMultipliers(Eigen::VectorXd energies);
	Eigen::VectorXcd getCoefficients(Eigen::VectorXd energies);
	std::complex<double> getStarFunction(Eigen::Vector3d & wavevector,
			long & iR);
	Eigen::Vector3cd getDerivativeStarFunction(Eigen::Vector3d & wavevector,
			long & iR);
	double getRoughnessFunction(Eigen::Vector3d position);
	const double coeff1 = 0.75; // 3/4
	const double coeff2 = 0.75;
	Eigen::Vector3d refWavevector;
	virtual std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
		diagonalizeFromCoords(Eigen::Vector3d & wavevector);
	double getEnergyFromCoords(Eigen::Vector3d & wavevector, long & bandIndex);
};

