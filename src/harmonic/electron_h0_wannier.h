#ifndef WANNIERH0_H
#define WANNIERH0_H

#include <math.h>
#include "eigen.h"
#include "harmonic.h"
#include "points.h"

class ElectronH0Wannier : public HarmonicHamiltonian {
public:
	ElectronH0Wannier(Eigen::Matrix3d & directUnitCell_,
			Eigen::MatrixXd & crystalVectors_,
			Eigen::VectorXd & vectorsDegeneracies_,
			Eigen::Tensor<std::complex<double>,3> & h0R_);

	std::tuple<Eigen::VectorXd,
		Eigen::Tensor<std::complex<double>,3>> diagonalize(Point & point);

	Eigen::Tensor<std::complex<double>,3> diagonalizeVelocity(Point & point);
    const bool hasEigenvectors = false;
protected:
    std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
    	diagonalizeFromCoords(Eigen::Vector3d & k);

    Eigen::MatrixXd crystalVectors;
    Eigen::VectorXd vectorsDegeneracies;
    Eigen::Matrix3d directUnitCell;
    Eigen::Tensor<std::complex<double>,3> h0R;

    long numBands;
    long numVectors;
};

#endif
