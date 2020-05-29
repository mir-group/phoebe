#ifndef PHONONH0_H
#define PHONONH0_H

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <math.h>
#include "crystal.h"
#include "harmonic.h"
#include "points.h"
#include "statistics.h"
#include "bandstructure.h"
#include "constants.h"

class PhononH0 : public HarmonicHamiltonian {
public:
	/** Class to store the force constants and diagonalize the dynamical matrix
	 * @param crystal: the object with the information on the crystal structure
	 * @param dielectricMatrix: 3x3 matrix with the dielectric matrix
	 * @param bornCharges: real tensor of size (numAtoms,3,3) with the Born
	 * effective charges
	 * @param forceConstants: a tensor of doubles with the force constants
	 * size is (meshx, meshy, meshz, 3, 3, numAtoms, numAtoms)
	 */
	PhononH0(Crystal & crystal,
			const Eigen::MatrixXd & dielectricMatrix_,
			const Eigen::Tensor<double, 3> & bornCharges_,
			const Eigen::Tensor<double, 7> & forceConstants_);
	// copy constructor
	PhononH0( const PhononH0 & that );
	// copy assignment
	PhononH0 & operator = ( const PhononH0 & that );

	~PhononH0();

	/** get the phonon energies (in Ry) at a single q-point.
	 * @param q: a q-point in cartesian coordinates.
	 * @return tuple(energies, eigenvectors): the energies are a double vector
	 * of size (numBands=3*numAtoms). Eigenvectors are a complex tensor of
	 * size (3,numAtoms,numBands). The eigenvector is rescaled by the
	 * sqrt(masses) (masses in rydbergs)
	 */
	template<typename T>
	std::tuple<Eigen::VectorXd,
		Eigen::Tensor<std::complex<double>,3>> diagonalize(Point<T> & point);

	/** get the phonon velocities (in atomic units) at a single q-point.
	 * @param q: a Point object with the wavevector coordinates.
	 * @return velocity(numBands,numBands,3): values of the velocity operator
	 * for this stata, in atomic units.
	 */
	template<typename T>
	Eigen::Tensor<std::complex<double>,3> diagonalizeVelocity(Point<T> &point);

	/** Impose the acoustic sum rule on force constants and Born charges
	 * @param sumRule: name of the sum rule to be used
	 * Currently supported values are akin to those from Quantum ESPRESSO
	 * i.e. "simple" (for a rescaling of the diagonal elements) or "crystal"
	 * (to find the closest matrix which satisfies the sum rule)
	 */
	void setAcousticSumRule(const std::string & sumRule);

	long getNumBands();
	Statistics getStatistics();

	template<typename T>
	FullBandStructure<T> populate(T & points, bool & withVelocities,
			bool & withEigenvectors);

	// this is almost the same as diagonalize, but takes in input the
	// cartesian coordinates
	// also, we return the eigenvectors aligned with the dynamical matrix,
	// and without the mass scaling.
	virtual std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> diagonalizeFromCoords(
				Eigen::Vector3d & q);
protected:
	Statistics statistics;

	Eigen::Vector3i getCoarseGrid();
	// internal variables

	// these 3 variables might be used for extending future functionalities.
	// for the first tests, they can be left at these default values
	// in the future, we might expose them to the user input
	bool na_ifc = false;
	bool loto_2d = false;
	bool frozenPhonon = false;

	bool hasDielectric;
	long numAtoms;
	long numBands;
	Eigen::MatrixXd directUnitCell;
	Eigen::MatrixXd reciprocalUnitCell;
	double latticeParameter;
	double volumeUnitCell;
	Eigen::MatrixXi atomicSpecies;
	Eigen::VectorXd speciesMasses;
	Eigen::MatrixXd atomicPositions;
	Eigen::MatrixXd dielectricMatrix;
	Eigen::Tensor<double,3> bornCharges;
	Eigen::Vector3i qCoarseGrid;
	Eigen::Tensor<double,7> forceConstants;
	Eigen::Tensor<double, 5> wscache;
	long nr1Big, nr2Big, nr3Big;

	// private methods, used to diagonalize the Dyn matrix
	void wsinit(const Eigen::MatrixXd & unitCell);
	double wsweight(const Eigen::VectorXd & r,
			const Eigen::MatrixXd & rws);
	void longRangeTerm(Eigen::Tensor<std::complex<double>,4> & dyn,
			const Eigen::VectorXd & q,
			const long sign);
	void nonAnaliticTerm(const Eigen::VectorXd & q,
			Eigen::Tensor<std::complex<double>,4> & dyn);
	void nonAnalIFC(const Eigen::VectorXd & q,
			Eigen::Tensor<std::complex<double>, 4> & f_of_q);
	void shortRangeTerm(Eigen::Tensor<std::complex<double>, 4> & dyn,
			const Eigen::VectorXd & q,
			Eigen::Tensor<std::complex<double>, 4> & f_of_q);
	std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> dyndiag(
			Eigen::Tensor<std::complex<double>,4> & dyn);

	// methods for sum rule
	void sp_zeu(Eigen::Tensor<double,3> & zeu_u,
			Eigen::Tensor<double,3> & zeu_v, double & scal);
};

template <typename T>
FullBandStructure<T> PhononH0::populate(T & points, bool & withVelocities,
		bool & withEigenvectors) {

	FullBandStructure<T> fullBandStructure(numBands, statistics,
			withVelocities, withEigenvectors, points);

	for ( long ik=0; ik<fullBandStructure.getNumPoints(); ik++ ) {
		Point point = fullBandStructure.getPoint(ik);

		auto [ens, eigvecs] = diagonalize(point);
		fullBandStructure.setEnergies(point, ens);

		if ( withVelocities) {
			auto vels = diagonalizeVelocity(point);
			fullBandStructure.setVelocities(point, vels);
		}
		if ( withEigenvectors ) {
			fullBandStructure.setEigenvectors(point, eigvecs);
		}
	}
	return fullBandStructure;
}

template <typename T>
std::tuple<Eigen::VectorXd,
		Eigen::Tensor<std::complex<double>,3>> PhononH0::diagonalize(
				Point<T> & point) {
	Eigen::Vector3d q = point.getCoords("cartesian");
	auto [energies, eigvecTemp] = diagonalizeFromCoords(q);

	//  displacements are eigenvectors divided by sqrt(speciesMasses)

	Eigen::Tensor<std::complex<double>,3> eigenvectors(3,numAtoms,numBands);
	for ( long iband=0; iband<numBands; iband++ ) {
		for ( long iat=0; iat<numAtoms; iat++ ) {
			for ( long ipol=0; ipol<3; ipol++ ) {
				auto ind = compress2Indeces(iat,ipol,numAtoms,3);
				eigenvectors(ipol,iat,iband) = eigvecTemp(ind, iband);
			}
		}
	}
	return {energies, eigenvectors};
}

template <typename T>
Eigen::Tensor<std::complex<double>,3> PhononH0::diagonalizeVelocity(
		Point<T> & point) {
	Eigen::Tensor<std::complex<double>,3> velocity(numBands,numBands,3);
	velocity.setZero();

	// if we are working at gamma, we set all velocities to zero.
	Eigen::Vector3d coords = point.getCoords("cartesian");
	if ( coords.norm() < 1.0e-6 )  {
		return velocity;
	}

	// get the eigenvectors and the energies of the q-point
	auto [energies,eigenvectors] = diagonalizeFromCoords(coords);

	// now we compute the velocity operator, diagonalizing the expectation
	// value of the derivative of the dynamical matrix.
	// This works better than doing finite differences on the frequencies.
	double deltaQ = 1.0e-8;
	for ( long i=0; i<3; i++ ) {
		// define q+ and q- from finite differences.
		Eigen::Vector3d qPlus = coords;
		Eigen::Vector3d qMins = coords;
		qPlus(i) += deltaQ;
		qMins(i) -= deltaQ;

		// diagonalize the dynamical matrix at q+ and q-
		auto [enPlus,eigPlus] = diagonalizeFromCoords(qPlus);
		auto [enMins,eigMins] = diagonalizeFromCoords(qMins);

		// build diagonal matrices with frequencies
		Eigen::MatrixXd enPlusMat(numBands,numBands);
		Eigen::MatrixXd enMinsMat(numBands,numBands);
		enPlusMat.setZero();
		enMinsMat.setZero();
		enPlusMat.diagonal() << enPlus;
		enMinsMat.diagonal() << enMins;

		// build the dynamical matrix at the two wavevectors
		// since we diagonalized it before, A = M.U.M*
		Eigen::MatrixXcd sqrtDPlus(numBands,numBands);
		sqrtDPlus = eigPlus * enPlusMat * eigPlus.adjoint();
		Eigen::MatrixXcd sqrtDMins(numBands,numBands);
		sqrtDMins = eigMins * enMinsMat * eigMins.adjoint();

		// now we can build the velocity operator
		Eigen::MatrixXcd der(numBands,numBands);
		der = (sqrtDPlus - sqrtDMins) / ( 2. * deltaQ );

		// and to be safe, we reimpose hermiticity
		der = 0.5 * ( der + der.adjoint() );

		// now we rotate in the basis of the eigenvectors at q.
		der = eigenvectors.adjoint() * der * eigenvectors;

		for ( long ib1=0; ib1<numBands; ib1++ ) {
			for ( long ib2=0; ib2<numBands; ib2++ ) {
				velocity(ib1,ib2,i) = der(ib1,ib2);
			}
		}
	}

	// turns out that the above algorithm has problems with degenerate bands
	// so, we diagonalize the velocity operator in the degenerate subspace,

	for ( long ib=0; ib<numBands; ib++ ) {

		// first, we check if the band is degenerate, and the size of the
		// degenerate subspace
		long sizeSubspace = 1;
		for ( long ib2=ib+1; ib2<numBands; ib2++ ) {
			// I consider bands degenerate if their frequencies are the same
			// within 0.0001 cm^-1
			if ( abs(energies(ib)-energies(ib2)) > 0.0001 / ryToCmm1 ) {
				break;
			}
			sizeSubspace += 1;
		}

		if ( sizeSubspace > 1 ) {
			Eigen::MatrixXcd subMat(sizeSubspace,sizeSubspace);
			// we have to repeat for every direction
			for ( long iCart=0; iCart<3; iCart++ ) {

				// take the velocity matrix of the degenerate subspace
				for ( long i=0; i<sizeSubspace; i++ ) {
					for ( long j=0; j<sizeSubspace; j++ ) {
						subMat(i,j) = velocity(ib+i,ib+j,iCart);
					}
				}

				// reinforce hermiticity
				subMat = 0.5 * ( subMat + subMat.adjoint());

				// diagonalize the subMatrix
				Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(subMat);
				Eigen::MatrixXcd newEigvecs = eigensolver.eigenvectors();

				// rotate the original matrix in the new basis
				// that diagonalizes the subspace.
				subMat = newEigvecs.adjoint() * subMat * newEigvecs;

				// reinforce hermiticity
				subMat = 0.5 * ( subMat + subMat.adjoint());

				// substitute back
				for ( long i=0; i<sizeSubspace; i++ ) {
					for ( long j=0; j<sizeSubspace; j++ ) {
						velocity(ib+i,ib+j,iCart) = subMat(i,j);
					}
				}
			}
		}

		// we skip the bands in the subspace, since we corrected them already
		ib += sizeSubspace - 1;
	}
	return velocity;
}



#endif
