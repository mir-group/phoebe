#ifndef POINTS_H
#define POINTS_H

#include "crystal.h"
#include "eigen.h"
#include "exceptions.h"

const int crystalCoords_ = 0;
const int cartesianCoords_ = 1;

//// Forward declarator of Points, because Point ants to store Points as member
//// and Points has methods returning Point. Here we avoid recursive dependency.
class Points;

/** Class used to pass a single wavevector
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

	/** Get the coordinates of the k-point
	 * @param basis: either "cartesian" or "crystal"
	 * @param inWignerSeitz: default false, if true, folds point in WS cell.
	 * @return coords: a 3d vector of coordinates
	 */
	Eigen::Vector3d getCoords(const int & basis=crystalCoords_,
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
	// constructors
	Points(Crystal & crystal_, const Eigen::Vector3i & mesh_,
			const Eigen::Vector3d & offset_=Eigen::Vector3d::Zero());
	Points(const Points & obj); // copy constructor
	Points & operator=(const Points & obj); // assignment operator

	// methods that mostly stay the same for all subclasses
	std::tuple<Eigen::Vector3i, Eigen::Vector3d> getMesh();
	long getNumPoints();
	Eigen::Vector3d crystalToCartesian(const Eigen::Vector3d & point);
	Eigen::Vector3d cartesianToCrystal(const Eigen::Vector3d & point);
	Eigen::Vector3d crystalToWS(const Eigen::Vector3d& pointCrystal,
			const int & basis);
	static std::tuple<Eigen::Vector3i, Eigen::Vector3d> findMesh(
			const Eigen::Matrix<double,3,Eigen::Dynamic> & points);
	Crystal & getCrystal();

	// methods to be overwritten in subclasses
	virtual Point getPoint(const long & index);
	virtual Eigen::Vector3d getPointCoords(const long & index,
			const int & basis=crystalCoords);
	virtual long getIndex(const Eigen::Vector3d & point);
	virtual double getWeight(const long & ik);

	// note: constexpr tells the compiler that the class member is
	// available at compilation time
	static constexpr const int crystalCoords = 0;
	static constexpr const int cartesianCoords = 1;
protected:
	void setMesh(const Eigen::Vector3i &mesh_,const Eigen::Vector3d & offset_);
	Crystal & crystal;
	Eigen::Vector3i mesh;
	Eigen::Vector3d offset;
	long numPoints = 0;
	// for Wigner Seitz folding
	Eigen::MatrixXd gVectors;
	Eigen::MatrixXi igVectors;

	// methods to be overwritten
	Eigen::Vector3d reduciblePoints(const long & idx);
};

class FullPoints: public Points {
public:
	FullPoints(Crystal & crystal_, const Eigen::Vector3i & mesh_,
			const Eigen::Vector3d & offset_=Eigen::Vector3d::Zero());
	FullPoints(const FullPoints & obj); // copy constructor
	FullPoints & operator=(const FullPoints & obj); // assignment operator

	Point getPoint(const long & index);
	long getIndexInverted(const long & ik);
};

class IrreduciblePoints: public Points {
protected:
	Eigen::VectorXi mapReducibleToIrreducible;
	Eigen::VectorXi mapIrreducibleToReducible;
	// points are internally stored in crystal coordinates
	Eigen::MatrixXd irreduciblePoints;
	Eigen::VectorXd irreducibleWeights;
	Eigen::VectorXi indexIrreduciblePoints;
	long numIrredPoints = 0;

	//	Eigen::Vector3d pointsCoords(const long & index);
	void setIrreduciblePoints();
public:
	IrreduciblePoints(Crystal & crystal_, const Eigen::Vector3i & mesh_,
			const Eigen::Vector3d & offset_=Eigen::Vector3d::Zero());
	IrreduciblePoints(const IrreduciblePoints & obj); // copy constructor
	IrreduciblePoints & operator=(const IrreduciblePoints & obj); // assignment

	long getIndex(const Eigen::Vector3d & point);
	Eigen::VectorXi getIndexReducibleFromIrreducible(const long & indexIrr);
	long getIndexIrreducibleFromReducible(const long & indexRed);
	Point getPoint(const long & index);

	long getNumPoints();
	double getWeight(const long & ik);
	Eigen::Vector3d getPointCoords(const long & index,
			const int & basis=crystalCoords);
};

class ActivePoints: public Points {
protected:
	FullPoints & parentPoints;
	Eigen::MatrixXd pointsList;

	VectorXl filteredToFullIndeces;
	long fullToFilteredIndeces(const long & indexIn);
public:
	// constructors
	ActivePoints(FullPoints & parentPoints_, VectorXl filter_);
	ActivePoints(const ActivePoints & obj); // copy constructor
	ActivePoints & operator=(const ActivePoints & obj); // assignment operator

	long getIndex(const Eigen::Vector3d & coords);
	Point getPoint(const long & index);
	long getIndexInverted(const long & ik);
	Eigen::Vector3d getPointCoords(const long & index,
			const int & basis=crystalCoords);
};

class PathPoints: public FullPoints {
public:
	PathPoints(Crystal & crystal_,
			const Eigen::Tensor<double,3> & pathExtrema,
			const double & delta);
	PathPoints(const PathPoints & obj); // copy constructor
	PathPoints & operator=(const PathPoints & obj); // assignment operator

	Point getPoint(const long & index);
	long getIndex(const Eigen::Vector3d & coords);
	long getIndexInverted(const long & ik);
	Eigen::Vector3d getPointCoords(const long & index,
			const int & basis=crystalCoords);
protected:
	Eigen::Matrix<double,3,Eigen::Dynamic> pointsList;
};

#endif
