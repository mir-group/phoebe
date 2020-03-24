#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <math.h>
#include "crystal.h"

class PhononH0 {
public:
	PhononH0(Crystal crystal,
			const Eigen::MatrixXd& dielectricMatrix_,
			const Eigen::Tensor<double, 3>& bornCharges_,
			const Eigen::Tensor<double, 7>& forceConstants_);
	void diagonalize(const Eigen::VectorXd& q, Eigen::VectorXd& energies,
			Eigen::Tensor<std::complex<double>,3>& eigenvectors);
	void setAcousticSumRule(const std::string sumRule);

//private:

	// internal variables
	bool na_ifc;
	bool loto_2d;
	bool frozenPhonon;
	bool hasDielectric;
	int numAtoms;
	int numBands;
	Eigen::MatrixXd directUnitCell;
	Eigen::MatrixXd reciprocalUnitCell;
	double latticeParameter;
	double volumeUnitCell;
	Eigen::MatrixXi atomicSpecies;
	Eigen::VectorXd speciesMasses;
	Eigen::MatrixXd atomicPositions;
	Eigen::MatrixXd dielectricMatrix;
	Eigen::Tensor<double,3> bornCharges;
	Eigen::VectorXi qCoarseGrid;
	Eigen::Tensor<double,7> forceConstants;
	Eigen::Tensor<double, 5> wscache;
	int nr1Big, nr2Big, nr3Big;

	// private methods, used to diagonalize the Dyn matrix
	void wsinit(const Eigen::MatrixXd& unitCell);
	double wsweight(const Eigen::VectorXd& r,
			const Eigen::MatrixXd& rws);
	void longRangeTerm(Eigen::Tensor<std::complex<double>,4>& dyn,
			const Eigen::VectorXd& q,
			const int sign);
	void nonAnaliticTerm(const Eigen::VectorXd& q,
			Eigen::Tensor<std::complex<double>,4>& dyn);
	void nonAnalIFC(const Eigen::VectorXd& q,
			Eigen::Tensor<std::complex<double>, 4>& f_of_q);
	void shortRangeTerm(Eigen::Tensor<std::complex<double>, 4>& dyn,
			const Eigen::VectorXd& q,
			Eigen::Tensor<std::complex<double>, 4>& f_of_q);
	void dyndiag(Eigen::Tensor<std::complex<double>,4>& dyn,
			Eigen::VectorXd& energies,
			Eigen::Tensor<std::complex<double>,3>& z);

	// methods for sum rule
	void sp_zeu(Eigen::Tensor<double,3>& zeu_u,
			Eigen::Tensor<double,3>& zeu_v,
			double& scal);
};
