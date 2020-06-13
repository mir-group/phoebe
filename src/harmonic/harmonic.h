#ifndef HARMONIC_H
#define HARMONIC_H

#include "points.h"
#include "particle.h"
#include "bandstructure.h"

class FullBandStructure;

class HarmonicHamiltonian {
public:
	HarmonicHamiltonian();

	virtual std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> diagonalize(
			Point & point);

	virtual Eigen::Tensor<std::complex<double>,3> diagonalizeVelocity(
			Point & point);

	const bool hasEigenvectors = true;
	virtual long getNumBands();
	virtual Particle getParticle();

	virtual FullBandStructure populate(Points & fullPoints,
			bool & withVelocities, bool & withEigenvectors);
protected:
	Particle particle;
	Eigen::Tensor<std::complex<double>,3> internalDiagonalizeVelocity(
				Eigen::Vector3d & coords, double & delta, double & threshold);

	std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
		diagonalizeFromCoords(Eigen::Vector3d & k);
	long numBands;
};

#endif
