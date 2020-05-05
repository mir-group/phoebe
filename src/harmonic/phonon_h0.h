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

	/** get the phonon energies (in Ry) at a single q-point.
	 * @param q: a q-point in cartesian coordinates.
	 * @return tuple(energies, eigenvectors): the energies are a double vector
	 * of size (numBands=3*numAtoms). Eigenvectors are a complex tensor of
	 * size (3,numAtoms,numBands). The eigenvector is rescaled by the
	 * sqrt(masses) (masses in rydbergs)
	 */
	virtual std::tuple<Eigen::VectorXd,
		Eigen::Tensor<std::complex<double>,3>> diagonalize(Point & point);

	/** get the phonon velocities (in atomic units) at a single q-point.
	 * @param q: a Point object with the wavevector coordinates.
	 * @return velocity(numBands,numBands,3): values of the velocity operator
	 * for this stata, in atomic units.
	 */
	virtual Eigen::Tensor<std::complex<double>,3> diagonalizeVelocity(Point & point);

	/** Impose the acoustic sum rule on force constants and Born charges
	 * @param sumRule: name of the sum rule to be used
	 * Currently supported values are akin to those from Quantum ESPRESSO
	 * i.e. "simple" (for a rescaling of the diagonal elements) or "crystal"
	 * (to find the closest matrix which satisfies the sum rule)
	 */
	void setAcousticSumRule(const std::string sumRule);

//	long getNumBands();
	Statistics getStatistics();
protected:
	Statistics statistics;

	Eigen::Vector3i getCoarseGrid();
	// internal variables
	bool na_ifc;
	bool loto_2d;
	bool frozenPhonon;
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
	void wsinit(const Eigen::MatrixXd& unitCell);
	double wsweight(const Eigen::VectorXd& r,
			const Eigen::MatrixXd& rws);
	void longRangeTerm(Eigen::Tensor<std::complex<double>,4>& dyn,
			const Eigen::VectorXd& q,
			const long sign);
	void nonAnaliticTerm(const Eigen::VectorXd& q,
			Eigen::Tensor<std::complex<double>,4>& dyn);
	void nonAnalIFC(const Eigen::VectorXd& q,
			Eigen::Tensor<std::complex<double>, 4>& f_of_q);
	void shortRangeTerm(Eigen::Tensor<std::complex<double>, 4>& dyn,
			const Eigen::VectorXd& q,
			Eigen::Tensor<std::complex<double>, 4>& f_of_q);
	std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> dyndiag(
			Eigen::Tensor<std::complex<double>,4>& dyn);

	// this is almost the same as diagonalize, but takes in input the
	// cartesian coordinates
	// also, we return the eigenvectors aligned with the dynamical matrix,
	// and without the mass scaling.
	virtual std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> diagonalizeFromCoords(
				Eigen::Vector3d & q);

	// methods for sum rule
	void sp_zeu(Eigen::Tensor<double,3>& zeu_u,
			Eigen::Tensor<double,3>& zeu_v,
			double& scal);
};

#endif
