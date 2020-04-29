#ifndef HARMONIC_H
#define HARMONIC_H

#include "points.h"

class HarmonicHamiltonian {
public:
	virtual std::tuple<Eigen::VectorXd,
		Eigen::Tensor<std::complex<double>,3>> diagonalize(Point & point);

	Eigen::Tensor<std::complex<double>,3> diagonalizeVelocity(
			Point & point);
	const bool hasEigenvectors = true;
protected:
	Eigen::Tensor<std::complex<double>,3>
			internalDiagonalizeVelocity(
				Eigen::Vector3d & coords, double & delta, double & threshold);

	std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
		diagonalizeFromCoords(Eigen::Vector3d & k);
	long numBands;
};

#endif
