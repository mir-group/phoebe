#ifndef POINTS_H
#define POINTS_H

#include "crystal.h"

class Points {
public:
	Points(Crystal& crystal_, const Eigen::Vector3i& mesh_);
	Eigen::Vector3i getMesh();

	int getNumPoints();
	int getNumIrredPoints();

	Eigen::MatrixXd getReduciblePoints(const std::string& basis);
	Eigen::MatrixXd getIrreduciblePoints(const std::string& basis)
	;
	Eigen::MatrixXd getReducibleFromIrreducible(Eigen::Vector3d point);
	Eigen::Vector3d getIrreducibleFromReducible(Eigen::Vector3d point);

	Eigen::VectorXi getIndexReducibleFromIrreducible(int indexIrr);
	int getIndexIrreducibleFromReducible(int indexRed);

	int getIndexInverted(const int& ik);
	int getIndex(const Eigen::Vector3d& point);

	Eigen::Vector3d getPoint(int index, std::string basis);
	Eigen::Vector3d getIrredPoint(int index, std::string basis);

	Eigen::Vector3d crystalToCartesian(const Eigen::Vector3d& point);
	Eigen::Vector3d cartesianToCrystal(const Eigen::Vector3d& point);
	Eigen::Vector3d crystalToWS(const Eigen::Vector3d& pointCrystal,
			const std::string& basis);
private:
	void setMesh(Eigen::Vector3i meshGrid);
	void setIrreduciblePoints();

	Crystal* crystal;
	Eigen::Vector3i mesh;
	// points are internally stored in crystal coordinates
	Eigen::MatrixXd irreduciblePoints;
	Eigen::MatrixXd reduciblePoints;
	Eigen::VectorXd irreducibleWeights;
	Eigen::VectorXi mapIrredPoints;
	Eigen::VectorXi indexIrreduciblePoints;
	int numPoints;
	int numIrredPoints;

	// for Wigner Seitz folding
	Eigen::MatrixXd gVectors;
	Eigen::MatrixXi igVectors;
};

class KPoints: public Points {
};

class QPoints: public Points {
};

#endif
