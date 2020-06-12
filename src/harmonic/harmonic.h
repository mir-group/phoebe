#ifndef HARMONIC_H
#define HARMONIC_H

#include "points.h"
#include "particle.h"
#include "bandstructure.h"

template<typename T> class FullBandStructure;

class HarmonicHamiltonian {
public:
	HarmonicHamiltonian();

	std::tuple<Eigen::VectorXd,
		Eigen::Tensor<std::complex<double>,3>> diagonalize(Point & point);

	Eigen::Tensor<std::complex<double>,3> diagonalizeVelocity(Point & point);

	const bool hasEigenvectors = true;
	virtual long getNumBands();
	virtual Particle getParticle();

	template<typename Arg>
	FullBandStructure<Arg> populate(Arg & fullPoints, bool & withVelocities,
			bool & withEigenvectors);
protected:
	Particle particle;
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
