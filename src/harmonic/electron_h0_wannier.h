#ifndef WANNIERH0_H
#define WANNIERH0_H

#include <math.h>
#include "eigen.h"
#include "harmonic.h"
#include "points.h"

class ElectronH0Wannier : public HarmonicHamiltonian {
public:
	ElectronH0Wannier(const Eigen::Matrix3d & directUnitCell_,
			const Eigen::MatrixXd & crystalVectors_,
			const Eigen::VectorXd & vectorsDegeneracies_,
			const Eigen::Tensor<std::complex<double>,3> & h0R_);

	std::tuple<Eigen::VectorXd,Eigen::MatrixXcd> diagonalize(
			Point & point);

	virtual Eigen::Tensor<std::complex<double>,3> diagonalizeVelocity(Point & point);
    const bool hasEigenvectors = true;
    Statistics getStatistics();
    long getNumBands();

    // copy constructor
    ElectronH0Wannier( const ElectronH0Wannier & that );
    // copy assignment
    ElectronH0Wannier & operator = ( const ElectronH0Wannier & that );
    // empty constructor
    ElectronH0Wannier();

	FullBandStructure populate(FullPoints & fullPoints,
			bool & withVelocities, bool & withEigenvectors);
protected:
    Statistics statistics;
    virtual std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
    	diagonalizeFromCoords(Eigen::Vector3d & k);

    Eigen::MatrixXd crystalVectors;
    Eigen::VectorXd vectorsDegeneracies;
    Eigen::Matrix3d directUnitCell;
    Eigen::Tensor<std::complex<double>,3> h0R;

    long numBands;
    long numVectors;
};

#endif
