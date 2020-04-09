#include <string>
#include "points.h"
#include "bandstructure.h"

class ElectronH0Spline {
public:
	ElectronH0Spline(Crystal& crystal_, ElBandStructure coarseBandStructure_,
			double cutoff_);
	double getEnergy(Point& point, int& bandIndex);
	Eigen::VectorXd getEnergies(Point& point);
	Eigen::Vector3d getGroupVelocity(Point& point, int& bandIndex);
	Eigen::MatrixXd getGroupVelocities(Point& point);
	ElBandStructure populateBandStructure(FullPoints* fullpoints=nullptr,
			IrreduciblePoints* irreduciblePoints=nullptr);
private:
	ElBandStructure coarseBandStructure;
	Crystal& crystal;
	FullPoints* coarsePoints;

	Eigen::MatrixXcd expansionCoefficients;

	int numBands;
	double cutoff;
	int numDataPoints;
	int numPositionVectors;
	double minDistance;
	std::vector<Eigen::Vector3d> positionVectors;
	void setPositionVectors();
	Eigen::VectorXcd getLagrangeMultipliers(Eigen::VectorXd energies);
	Eigen::VectorXcd getCoefficients(Eigen::VectorXd energies);
	std::complex<double> getStarFunction(Eigen::Vector3d& wavevector,
			Eigen::Vector3d& position);
	Eigen::Vector3cd getDerivativeStarFunction(Eigen::Vector3d& wavevector,
			Eigen::Vector3d& position);
	double getRoughnessFunction(Eigen::Vector3d position);
	const double coeff1 = 0.75; // 3/4
	const double coeff2 = 0.75;
	Eigen::Vector3d refWavevector;
};

