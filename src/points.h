#ifndef POINTS_H
#define POINTS_H

#include "crystal.h"
#include "eigen.h"

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

	// copy constructor
	Point( const Point & that );
	// copy assignment
	Point & operator = ( const Point & that );
	~Point();

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

    long getIndex();
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
	Points(const Points & obj); // copy constructor
	Points & operator=(const Points & obj); // assignment operator
	~Points();

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
	Crystal & getCrystal();
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
//	Eigen::Vector3d pointsCoords(const long & index);
public:
	long getIndexInverted(const long & ik);
	FullPoints(Crystal & crystal_, const Eigen::Vector3i & mesh_,
			const Eigen::Vector3d & offset_=Eigen::Vector3d::Zero());
	FullPoints(const FullPoints & obj); // copy constructor
	FullPoints & operator=(const FullPoints & obj); // assignment operator
	~FullPoints();
};

class IrreduciblePoints: public Points {
protected:
	bool useIrreducible = true;
//	Eigen::Vector3d pointsCoords(const long & index);
public:
	IrreduciblePoints(Crystal & crystal_, const Eigen::Vector3i & mesh_,
			const Eigen::Vector3d & offset_=Eigen::Vector3d::Zero());
	long getIndex(const Eigen::Vector3d & point);
	// long getIndexInverted(const long& ik); // not sure if this makes sense
	Eigen::VectorXi getIndexReducibleFromIrreducible(const long & indexIrr);
	long getIndexIrreducibleFromReducible(const long & indexRed);
};

class ActivePoints: public Points {
protected:
	FullPoints & parentPoints;
	Eigen::MatrixXd pointsList;

	VectorXl filteredToFullIndeces;
	long fullToFilteredIndeces(const long & indexIn);
	Eigen::Vector3d pointsCoords(const long & index);
public:
	long getIndexInverted(const long & ik);
	long getIndex(const Eigen::Vector3d & coords);
//	ActivePoints() = default; // default constructor, only to be used for
	// temporary object initializations.
	ActivePoints(FullPoints & parentPoints_, VectorXl filter_);
	ActivePoints(const ActivePoints & obj); // copy constructor
	ActivePoints & operator=(const ActivePoints & obj); // assignment operator
};

class PathPoints: public FullPoints {
public:
	PathPoints(Crystal & crystal_,
			const Eigen::Tensor<double,3> & pathExtrema,
			const double & delta);

	long getIndex(const Eigen::Vector3d & coords);
	long getIndexInverted(const long & ik);
protected:
	Eigen::Vector3d pointsCoords(const long & index);
	Eigen::Matrix<double,3,Eigen::Dynamic> pointsList;
};

#endif
