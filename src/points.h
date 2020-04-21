#ifndef POINTS_H
#define POINTS_H

#include "crystal.h"

// Forward declarator of Points, because Point ants to store Points as member
// and Points has methods returning Point. Here we avoid recursive dependency.
class Points;

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
	Point(long index_, Eigen::Vector3d umklappVector, Points & points_);

	/** Get the coordinates of the k-point
	 * @param basis: either "cartesian" or "crystal"
	 * @param inWignerSeitz: default false, if true, folds point in WS cell.
	 * @return coords: a 3d vector of coordinates
	 */
	Eigen::Vector3d getCoords(const std::string & basis="crystal",
			const bool & inWignerSeitz=false);

	/** Get the weight of the k-point (used for integrations over the BZ.
	 * @return weight: a double.
	 */
	double getWeight();

    Point operator + (Point & b);
    Point operator - (Point & b);

    bool hasUmklapp();
private:
	Eigen::Vector3d umklappVector;
	long index;
	Points & points;
};

class Points {
public:
	Points(Crystal & crystal_, const Eigen::Vector3i & mesh_,
			const Eigen::Vector3d & offset_=Eigen::Vector3d::Zero(),
			const bool useIrreducible_=false);

	std::tuple<Eigen::Vector3i, Eigen::Vector3d> getMesh();
	long getNumPoints();

	Point getPoint(const long & index);
	Eigen::Vector3d getPointCoords(const long & index,
			const std::string & basis="crystal");
	long getIndex(const Eigen::Vector3d & point);
	std::vector<Eigen::Vector3d> getPointsCoords(const std::string & basis="crystal");
	double getWeight(const long & ik);

	Eigen::Vector3d crystalToCartesian(const Eigen::Vector3d & point);
	Eigen::Vector3d cartesianToCrystal(const Eigen::Vector3d & point);
	Eigen::Vector3d crystalToWS(const Eigen::Vector3d& pointCrystal,
			const std::string& basis);

	static std::tuple<Eigen::Vector3i, Eigen::Vector3d> findMesh(
			const Eigen::MatrixXd & points);
protected:
	void setMesh(const Eigen::Vector3i & mesh_,
			const Eigen::Vector3d & offset_);
	void setIrreduciblePoints();

	Eigen::Vector3d reduciblePoints(const long & idx);
	Eigen::Vector3d pointsCoords(const long & index);

	Crystal& crystal;
	Eigen::Vector3i mesh;
	Eigen::Vector3d offset;
	// points are internally stored in crystal coordinates
	Eigen::MatrixXd irreduciblePoints;
	Eigen::VectorXd irreducibleWeights;
	Eigen::VectorXi mapReducibleToIrreducible;
	Eigen::VectorXi mapIrreducibleToReducible;
	Eigen::VectorXi indexIrreduciblePoints;
	long numPoints = 0;
	long numIrredPoints = 0;
	bool useIrreducible = false;

	// for Wigner Seitz folding
	Eigen::MatrixXd gVectors;
	Eigen::MatrixXi igVectors;
};

class FullPoints: public Points {
protected:
	bool useIrreducible = false;
	Eigen::Vector3d pointsCoords(const long & index);
public:
	long getIndexInverted(const long & ik);
	FullPoints(Crystal & crystal_, const Eigen::Vector3i & mesh_,
			const Eigen::Vector3d & offset_=Eigen::Vector3d::Zero());
};

class IrreduciblePoints: public Points {
protected:
	bool useIrreducible = true;
	Eigen::Vector3d pointsCoords(const long & index);
public:
	IrreduciblePoints(Crystal & crystal_, const Eigen::Vector3i & mesh_,
			const Eigen::Vector3d & offset_=Eigen::Vector3d::Zero());
	long getIndex(const Eigen::Vector3d & point);
	// long getIndexInverted(const long& ik); // not sure if this makes sense
	Eigen::VectorXi getIndexReducibleFromIrreducible(const long & indexIrr);
	long getIndexIrreducibleFromReducible(const long & indexRed);
};

class ActivePoints: public FullPoints {
protected:
	FullPoints & parentPoints;
	Eigen::MatrixXd pointsList;

	Eigen::VectorXi filteredToFullIndeces;
	long fullToFilteredIndeces(const long & indexIn);

	Eigen::Vector3d pointsCoords(const long & index);
public:
	long getIndexInverted(const long & ik);
	long getIndex(const Eigen::Vector3d & coords);
	ActivePoints(Crystal & crystal_, FullPoints & parentPoints_,
			Eigen::VectorXi filter_, const Eigen::Vector3i & mesh_,
			const Eigen::Vector3d & offset_=Eigen::Vector3d::Zero());
};

#endif
