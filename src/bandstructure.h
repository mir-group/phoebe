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

	FullPoints* fullPoints = nullptr;
	IrreduciblePoints* irreduciblePoints = nullptr;
	int getIndex(Eigen::Vector3d& pointCoords);

	double homo;
	int numValenceElectrons;
public:
	BaseBandStructure(int numBands_, FullPoints* fullPoints_=nullptr,
			IrreduciblePoints* irreduciblePoints_=nullptr);
	Point getPoint(const int& pointIndex);
	int getNumPoints();
	State getStateFromPointIndex(int index);
	void setChemicalPotential(double chemPot);
	void setTemperature(double temp);
	void populate();
	void setEnergies(Eigen::Vector3d& point, Eigen::VectorXd& energies_);
	void setEnergies(Point& point, Eigen::VectorXd& energies_);
	void setVelocities(Eigen::Vector3d& pointCoords,
			Eigen::Tensor<std::complex<double>,3>& velocities_);
	int getNumBands();
	bool hasIrreduciblePoints();
	void setNumValenceElectrons(int numElectrons);
	void setHomo(double homo);
};

class ElBandStructure: public BaseBandStructure {
protected:
	std::string statistics = "fermi";
public:
	ElState getStateFromPointIndex(int index);
	Eigen::VectorXd getBandEnergies(int& bandIndex);
	ElBandStructure(int numBands_, FullPoints* fullPoints_=nullptr,
			IrreduciblePoints* irreduciblePoints_=nullptr);
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

