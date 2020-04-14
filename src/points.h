#ifndef POINTS_H
#define POINTS_H

#include "crystal.h"

/** Class used to pass a single wavevector
 *
 */
class Point {
public:
	/** Constructor
	 * @param crystalCoords: crystal coordinates of the point
	 * @param crystalCoordsWS: crystal coordinates of the point folded in the
	 * Wigner Seitz cell
	 * @param reciprocalUnitCell: a 3x3 matrix with lattice vector in rec space
	 */
	Point(Eigen::Vector3d crystalCoords_, Eigen::Vector3d crystalCoordsWS_,
			const Eigen::Matrix3d& reciprocalUnitCell_);

	/** Get the coordinates of the k-point
	 * @param basis: either "cartesian" or "crystal"
	 * @param inWignerSeitz: default false, if true, folds point in WS cell.
	 * @return coords: a 3d vector of coordinates
	 */
	Eigen::Vector3d getCoords(std::string basis="crystal",
			bool inWignerSeitz=false);

	/** Get the weight of the k-point (used for integrations over the BZ.
	 * @return weight: a double.
	 */
	double getWeight();
private:
	Eigen::Vector3d crystalCoords;
	Eigen::Vector3d crystalCoordsWS;
	const Eigen::Matrix3d& reciprocalUnitCell;
	double weight;
};

class Points {
public:
	Points(Crystal& crystal_, const Eigen::Vector3i& mesh_,
			const Eigen::Vector3d& offset_=Eigen::Vector3d::Zero(),
			const bool useIrreducible_=false);

	std::tuple<Eigen::Vector3i, Eigen::Vector3d> getMesh();
	int getNumPoints();

	Point getPoint(const int& index);
	Eigen::Vector3d getPointCoords(int& index, std::string basis="crystal");
	int getIndex(const Eigen::Vector3d& point);
	std::vector<Eigen::Vector3d> getPointsCoords(std::string basis="crystal");

	Eigen::Vector3d crystalToCartesian(const Eigen::Vector3d& point);
	Eigen::Vector3d cartesianToCrystal(const Eigen::Vector3d& point);
	Eigen::Vector3d crystalToWS(const Eigen::Vector3d& pointCrystal,
			const std::string& basis);
	static std::tuple<Eigen::Vector3i, Eigen::Vector3d> findMesh(
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
	Eigen::Vector3d pointsCoords(const int& index);

	Crystal& crystal;
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
	Eigen::Vector3d pointsCoords(const int& index);
public:
	int getIndexInverted(const int& ik);
	FullPoints(Crystal& crystal_, const Eigen::Vector3i& mesh_,
			const Eigen::Vector3d& offset_=Eigen::Vector3d::Zero());
};

class IrreduciblePoints: public Points {
protected:
	bool useIrreducible = true;
	Eigen::Vector3d pointsCoords(const int& index);
public:
	IrreduciblePoints(Crystal& crystal_, const Eigen::Vector3i& mesh_,
			const Eigen::Vector3d& offset_=Eigen::Vector3d::Zero());
	int getIndex(const Eigen::Vector3d& point);
	// int getIndexInverted(const int& ik); // not sure if this makes sense
	Eigen::VectorXi getIndexReducibleFromIrreducible(int indexIrr);
	int getIndexIrreducibleFromReducible(int indexRed);
};

#endif
