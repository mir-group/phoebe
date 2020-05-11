#ifndef WANNIERH0_H
#define WANNIERH0_H

#include <math.h>
#include "eigen.h"
#include "harmonic.h"
#include "points.h"
#include "bandstructure.h"

class ElectronH0Wannier : public HarmonicHamiltonian {
public:
	ElectronH0Wannier(const Eigen::Matrix3d & directUnitCell_,
			const Eigen::MatrixXd & crystalVectors_,
			const Eigen::VectorXd & vectorsDegeneracies_,
			const Eigen::Tensor<std::complex<double>,3> & h0R_);

	std::tuple<Eigen::VectorXd,Eigen::MatrixXcd> diagonalize(
			Point & point);

	virtual Eigen::Tensor<std::complex<double>,3> diagonalizeVelocity(
			Point & point);
    const bool hasEigenvectors = true;
    Statistics getStatistics();
    long getNumBands();

    // copy constructor
    ElectronH0Wannier( const ElectronH0Wannier & that );
    // copy assignment
    ElectronH0Wannier & operator = ( const ElectronH0Wannier & that );
    // empty constructor
    ElectronH0Wannier();

    template<typename Arg>
    FullBandStructure<Arg> populate(Arg & fullPoints, bool & withVelocities,
    		bool &withEigenvectors);
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

template<typename Arg>
FullBandStructure<Arg> ElectronH0Wannier::populate(Arg & fullPoints,
		bool & withVelocities, bool & withEigenvectors) {

	FullBandStructure<Arg> fullBandStructure(numBands, statistics,
			withVelocities, withEigenvectors, fullPoints);

	for ( long ik=0; ik<fullBandStructure.getNumPoints(); ik++ ) {
		Point point = fullBandStructure.getPoint(ik);
		auto [ens, eigvecs] = diagonalize(point);
		fullBandStructure.setEnergies(point, ens);
		if ( withVelocities ) {
			auto vels = diagonalizeVelocity(point);
			fullBandStructure.setVelocities(point, vels);
		}
		//TODO: I must fix the different shape of eigenvectors w.r.t. phonons
//		if ( withEigenvectors ) {
//			fullBandStructure.setEigenvectors(point, eigvecs);
//		}
	}
	return fullBandStructure;
}

#endif
