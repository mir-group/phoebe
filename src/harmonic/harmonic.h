#ifndef HARMONIC_H
#define HARMONIC_H

#include "points.h"
#include "statistics.h"
#include "bandstructure.h"

template<typename T> class FullBandStructure;

class HarmonicHamiltonian {
public:
	HarmonicHamiltonian();
	virtual Eigen::VectorXd diagonalizeEnergy(Point & point);
	std::tuple<Eigen::VectorXd,
		Eigen::Tensor<std::complex<double>,3>> diagonalize(Point & point);

	virtual Eigen::Tensor<std::complex<double>,3> diagonalizeVelocity(
			Point & point);
	const bool hasEigenvectors = true;
	virtual long getNumBands();
	virtual Statistics getStatistics();

	template<typename Arg>
	FullBandStructure<Arg> populate(Arg & fullPoints, bool & withVelocities,
			bool & withEigenvectors);
protected:
	Statistics statistics;
	Eigen::Tensor<std::complex<double>,3> internalDiagonalizeVelocity(
				Eigen::Vector3d & coords, double & delta, double & threshold);

	std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
		diagonalizeFromCoords(Eigen::Vector3d & k);
	long numBands;
};

template<typename Arg>
FullBandStructure<Arg> HarmonicHamiltonian::populate(Arg & fullPoints,
		bool & withVelocities, bool & withEigenvectors) {
	Error e("base populate not implemented");
	(void) fullPoints;
	(void) withVelocities;
	(void) withEigenvectors;
	FullBandStructure<Arg> t;
	return t;
}

#endif
