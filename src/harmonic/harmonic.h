#ifndef HARMONIC_H
#define HARMONIC_H

#include "points.h"

class HarmonicHamiltonian {
public:
	std::tuple<Eigen::VectorXd,
		Eigen::Tensor<std::complex<double>,3>> diagonalize(Point & point);

	Eigen::Tensor<std::complex<double>,3> diagonalizeVelocity(
			Point & point);
	const bool hasEigenvectors = true;
};

#endif
