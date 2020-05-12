#ifndef H0FOURIER_H
#define H0FOURIER_H

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
			FullBandStructure<FullPoints> coarseBandStructure_, double cutoff);

	template<typename T>
	std::tuple<Eigen::VectorXd, Eigen::Tensor<std::complex<double>,3>>
		diagonalize(Point<T> & point);

	template<typename T>
	Eigen::Tensor<std::complex<double>,3> diagonalizeVelocity(
				Point<T> & point);

	Statistics getStatistics();

	template<typename T>
	FullBandStructure<T> populate(T & fullPoints, bool & withVelocities,
			bool & withEigenvectors);

	long getNumBands();
protected:
	Crystal & crystal;
	FullBandStructure<FullPoints> coarseBandStructure;
	FullPoints coarsePoints;
	Statistics statistics;

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
std::tuple<Eigen::VectorXd, Eigen::Tensor<std::complex<double>,3>>
		ElectronH0Fourier::diagonalize(Point<T> & point) {

	Eigen::Vector3d coords = point.getCoords("cartesian");
	auto [energies,x] = diagonalizeFromCoords(coords);

	// this is to return something aligned with the phonon case
	// One should investigate how to return a null pointer
	Eigen::Tensor<std::complex<double>,3> eigvecs;
	eigvecs.setZero();

	return {energies,eigvecs};
}

template<typename T>
Eigen::Tensor<std::complex<double>,3> ElectronH0Fourier::diagonalizeVelocity(
			Point<T> & point) {
	Eigen::Tensor<std::complex<double>,3> velocity(numBands,numBands,3);
	velocity.setZero();
	Eigen::Vector3d coords = point.getCoords("cartesian");
	for ( long ib=0; ib<numBands; ib++ ) {
		Eigen::Vector3d v = getGroupVelocityFromCoords(coords,ib);
		for ( long i=0; i<3; i++ ) {
			velocity(ib,ib,i) = v(i);
		}
	}
	return velocity;
}

template<typename T>
FullBandStructure<T> ElectronH0Fourier::populate(T & fullPoints,
		bool & withVelocities, bool & withEigenvectors) {

	FullBandStructure<T> fullBandStructure(numBands, statistics,
			withVelocities, withEigenvectors, fullPoints);

	for ( long ik=0; ik<fullBandStructure.getNumPoints(); ik++ ) {
		Point<T> point = fullBandStructure.getPoint(ik);
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
