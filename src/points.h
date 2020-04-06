#ifndef POINTS_H
#define POINTS_H

#include "crystal.h"

class Points {
public:
	Points(Crystal& crystal_, const Eigen::Vector3i& mesh_,
			const Eigen::Vector3d& offset_=Eigen::Vector3d::Zero(),
			const bool useIrreducible_=false);

	std::tuple<Eigen::Vector3i, Eigen::Vector3d> getMesh();
	int getNumPoints();

	Eigen::Vector3d getPoint(int index, std::string basis="crystal");
	int getIndex(const Eigen::Vector3d& point);
	std::vector<Eigen::Vector3d> getPoints(std::string basis="crystal");

	Eigen::Vector3d crystalToCartesian(const Eigen::Vector3d& point);
	Eigen::Vector3d cartesianToCrystal(const Eigen::Vector3d& point);
	Eigen::Vector3d crystalToWS(const Eigen::Vector3d& pointCrystal,
			const std::string& basis);
	std::tuple<Eigen::Vector3i, Eigen::Vector3d> findMesh(
			Eigen::MatrixXd& points);
	double getWeight(int ik);

	const Eigen::Vector3d iteratePoints();

	// iterators over points
//	std::vector<Eigen::Vector3d>::iterator iteratePoints;
//	std::vector<Eigen::Vector3d>::iterator iteratePoints();

protected:
	void setMesh(Eigen::Vector3i mesh_, Eigen::Vector3d offset_);
	void setIrreduciblePoints();

	Eigen::Vector3d reduciblePoints(const int idx);
	Eigen::Vector3d pointsCoords(int& index);

	Crystal* crystal;
	Eigen::Vector3i mesh;
	Eigen::Vector3d offset;
	// points are internally stored in crystal coordinates
	Eigen::MatrixXd irreduciblePoints;
	Eigen::VectorXd irreducibleWeights;
	Eigen::VectorXi mapReducibleToIrreducible;
	Eigen::VectorXi mapIrreducibleToReducible;
	Eigen::VectorXi indexIrreduciblePoints;
	int numPoints;
	int numIrredPoints;
	bool useIrreducible = false;

	// for Wigner Seitz folding
	Eigen::MatrixXd gVectors;
	Eigen::MatrixXi igVectors;
};

class FullPoints: public Points {
protected:
	bool useIrreducible = false;
	Eigen::Vector3d pointsCoords(int& index);
public:
	int getIndexInverted(const int& ik);
	FullPoints(Crystal& crystal_, const Eigen::Vector3i& mesh_,
			const Eigen::Vector3d& offset_=Eigen::Vector3d::Zero());
};

class IrreduciblePoints: public Points {
protected:
	bool useIrreducible = true;
	Eigen::Vector3d pointsCoords(int& index);
public:
	IrreduciblePoints(Crystal& crystal_, const Eigen::Vector3i& mesh_,
			const Eigen::Vector3d& offset_=Eigen::Vector3d::Zero());
	int getIndex(const Eigen::Vector3d& point);
	// int getIndexInverted(const int& ik); // not sure if this makes sense
	Eigen::VectorXi getIndexReducibleFromIrreducible(int indexIrr);
	int getIndexIrreducibleFromReducible(int indexRed);
};

#endif
