#include "points.h"

class State {
public:
	State(Point& point_,
			Eigen::VectorXd& energies_,
			Eigen::Tensor<std::complex<double>,3>& velocities_,
			Eigen::VectorXd& dnde_, Eigen::VectorXd& dndt_);
	double getEnergy(const int bandIndex);
	Eigen::VectorXd getEnergies();
	Eigen::Vector3d getGroupVelocity(const int bandIndex);
	Eigen::MatrixXd getGroupVelocities();
	Eigen::Tensor<std::complex<double>,3> getVelocities();
	Eigen::VectorXd getDndt();
	Eigen::VectorXd getDnde();
	double getWeight();
protected:
	// pointers to the bandstructure, I don't want to duplicate storage here
	Point& point; // in cryst coords
	Eigen::VectorXd energies;
	Eigen::VectorXd dndt;
	Eigen::VectorXd dnde;
	Eigen::Tensor<std::complex<double>,3> velocities;
	int numBands;
};

class ElState: public State {
public:
	ElState(Point& point_,
			Eigen::VectorXd& energies_,
			Eigen::Tensor<std::complex<double>,3>& velocities_,
			Eigen::VectorXd& dnde_, Eigen::VectorXd& dndt_);
};

class PhState: public State {
public:
	PhState(Point& point_,
			Eigen::VectorXd& energies_,
			Eigen::Tensor<std::complex<double>,3>& eigenvectors_,
			Eigen::Tensor<std::complex<double>,3>& velocities_,
			Eigen::VectorXd& dnde_,
			Eigen::VectorXd& dndt_);
	Eigen::Tensor<std::complex<double>,3> getEigenvectors();
protected:
	Eigen::Tensor<std::complex<double>,3> eigenvectors;
};
