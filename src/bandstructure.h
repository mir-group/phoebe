#include "points.h"
#include "phononH0.h"
#include "state.h"

class BaseBandStructure {
protected:
	Eigen::MatrixXd energies;
	Eigen::MatrixXcd velocities;
	int numBands = 0;
	double chemicalPotential = 0.;
	double temperature = 0.;
	Eigen::MatrixXd dndt;
	Eigen::MatrixXd dnde;
	void setOccupations();
	std::string statistics = "";
	bool useIrreducible = false;

	int getNumPoints();
	int getIndex(Eigen::Vector3d& pointCoords);
	Point getPoint(const int& pointIndex);

	FullPoints* fullPoints = nullptr;
	IrreduciblePoints* irreduciblePoints = nullptr;
public:
	BaseBandStructure(int numBands_, FullPoints* fullPoints_=nullptr,
			IrreduciblePoints* irreduciblePoints_=nullptr);

	State getStateFromPointIndex(int index);
	void setChemicalPotential(double chemPot);
	void setTemperature(double temp);
	void populate();
	void setEnergies(Eigen::Vector3d& pointCoords,
			Eigen::VectorXd& energies_);
	void setVelocities(Eigen::Vector3d& pointCoords,
			Eigen::Tensor<std::complex<double>,3>& velocities_);
};

class ElBandStructure: public BaseBandStructure {
protected:
	std::string statistics = "fermi";
public:
	ElState getStateFromPointIndex(int index);
};

class PhBandStructure: public BaseBandStructure {
public:
	PhBandStructure(int numBands_, FullPoints* fullPoints_=nullptr,
			IrreduciblePoints* irreduciblePoints_=nullptr);
	void populate(PhononH0 phononH0);
	void setEigenvectors(Eigen::Vector3d& pointCoords,
			Eigen::Tensor<std::complex<double>,3>& eigenvectors_);
	PhState getStateFromPointIndex(int index);
protected:
	std::string statistics = "bose";
	int numAtoms;
	const double chemicalPotential = 0.;
	Eigen::MatrixXcd eigenvectors;
};

