#ifndef H0FOURIER_H
#define H0FOURIER_H

#include <string>
#include "points.h"
#include "bandstructure.h"
#include "harmonic.h"

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
	 * in the plane-wave interpolation.
	 */
	ElectronH0Fourier(Crystal & crystal_, FullPoints coarsePoints_,
			FullBandStructure<FullPoints> coarseBandStructure_, double cutoff);

	/** get the electronic energies (in Ry) at a single k-point.
	 * Energies don't have any reference value, and must be used in connection
	 * with a chemical potential.
	 * @param k: a point object with the wavevector. Must have the cartesian
	 * coordinates of the wavevector.
	 * @return tuple(energies, eigenvectors): the energies are a double vector
	 * of size (numBands). Eigenvectors of size (numBands,numBands), but are
	 * simply set to zero, since there is no diagonalization happening here.
	 */
	std::tuple<Eigen::VectorXd, Eigen::Tensor<std::complex<double>,3>>
			diagonalize(Point & point);

	/** get the electron velocities (in atomic units) at a single k-point.
	 * @param k: a Point object with the wavevector coordinates.
	 * @return velocity(numBands,numBands,3): values of the velocity operator
	 * for this state, in atomic units. Note that the off-diagonal matrix
	 * elements are set to zero, because this kind of interpolation, at the
	 * moment, doesn't have any information on the off-diagonal elements.
	 */
	Eigen::Tensor<std::complex<double>,3> diagonalizeVelocity(Point & point);

    /** Method to return that the underlying is that of an electronic Fermion.
     */
	Particle getParticle();

	/** This method constructs an electron bandstructure.
	 * @param points: the object with the list/mesh of wavevectors
	 * @param withVelocities: if true, compute the electron velocity operator.
	 * @param withEigenvectors: can only be false, as there are no eigenvectors
	 * with this kind of interpolation.
	 * @return FullBandStructure: the bandstructure object containing the
	 * complete electronic band structure.
	 */
	template<typename T>
	FullBandStructure<T> populate(T & fullPoints, bool & withVelocities,
			bool & withEigenvectors);

	/** Get the total number of bands available at ech wavevector.
	 *
	 */
	long getNumBands();
protected:
	Crystal & crystal;
	FullBandStructure<FullPoints> coarseBandStructure;
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
	Eigen::Vector3d getGroupVelocityFromCoords(Eigen::Vector3d & wavevector,
			long & bandIndex);
};

template<typename T>
FullBandStructure<T> ElectronH0Fourier::populate(T & fullPoints,
		bool & withVelocities, bool & withEigenvectors) {

	FullBandStructure<T> fullBandStructure(numBands, particle,
			withVelocities, withEigenvectors, fullPoints);

	for ( long ik=0; ik<fullBandStructure.getNumPoints(); ik++ ) {
		Point point = fullBandStructure.getPoint(ik);
		auto [ens, eigvecs] = diagonalize(point);
		fullBandStructure.setEnergies(point, ens);
		if ( withVelocities ) {
			auto vels = diagonalizeVelocity(point);
			fullBandStructure.setVelocities(point, vels);
		}
	}
	return fullBandStructure;
}

#endif
