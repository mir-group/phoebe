#ifndef POINTS_H
#define POINTS_H

#include "crystal.h"
#include "eigen.h"
#include "exceptions.h"

//// Forward declarator of Points, because Point ants to store Points as member
//// and Points has methods returning Point. Here we avoid recursive dependency.
//class Points;

/** Class used to pass a single wavevector
 *
 */
template<typename T>
class Point {
public:
	/** Constructor
	 * @param crystalCoords: crystal coordinates of the point
	 * @param crystalCoordsWS: crystal coordinates of the point folded in the
	 * Wigner Seitz cell
	 * @param reciprocalUnitCell: a 3x3 matrix with lattice vector in rec space
	 */
	Point(long index_, Eigen::Vector3d umklappVector, T & points_);

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
	T & points;
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
			const std::string& basis);
	static std::tuple<Eigen::Vector3i, Eigen::Vector3d> findMesh(
			const Eigen::Matrix<double,3,Eigen::Dynamic> & points);
	Crystal & getCrystal();

	// methods to be overwritten in subclasses
	Point<Points> getPoint(const long & index);
	Eigen::Vector3d getPointCoords(const long & index,
			const std::string & basis="crystal");
	long getIndex(const Eigen::Vector3d & point);
	double getWeight(const long & ik);

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

	Point<FullPoints> getPoint(const long & index);
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
	Point<IrreduciblePoints> getPoint(const long & index);

	long getNumPoints();
	double getWeight(const long & ik);
	Eigen::Vector3d getPointCoords(const long & index,
			const std::string & basis="crystal");
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
	Point<ActivePoints> getPoint(const long & index);
	long getIndexInverted(const long & ik);
	Eigen::Vector3d getPointCoords(const long & index,
			const std::string & basis="crystal");
};

class PathPoints: public FullPoints {
public:
	PathPoints(Crystal & crystal_,
			const Eigen::Tensor<double,3> & pathExtrema,
			const double & delta);
	PathPoints(const PathPoints & obj); // copy constructor
	PathPoints & operator=(const PathPoints & obj); // assignment operator

	Point<PathPoints> getPoint(const long & index);
	long getIndex(const Eigen::Vector3d & coords);
	long getIndexInverted(const long & ik);
	Eigen::Vector3d getPointCoords(const long & index,
			const std::string & basis="crystal");
protected:
	Eigen::Matrix<double,3,Eigen::Dynamic> pointsList;
};

template<typename T>
Point<T>::Point(long index_, Eigen::Vector3d umklappVector_, T & points_)
		: points(points_) {
	umklappVector = umklappVector_;
	index = index_;
}

// copy constructor
template<typename T>
Point<T>::Point( const Point & that ) : umklappVector(that.umklappVector),
		index(that.index), points(that.points) {
}

// copy assignment
template<typename T>
Point<T> & Point<T>::operator = ( const Point & that ) {
	if ( this != &that ) {
		umklappVector = that.umklappVector;
		index = that.index;
		points = that.points;
	}
	return *this;
}

template<typename T>
Point<T>::~Point() {
}

template<typename T>
long Point<T>::getIndex() {
	return index;
}

template<typename T>
Eigen::Vector3d Point<T>::getCoords(const std::string & basis,
		const bool & inWignerSeitz) {
	if ( ( basis != "crystal" ) && ( basis != "cartesian" ) ) {
		Error e("Point getCoordinates: basis must be crystal or cartesian", 1);
	}
	Eigen::Vector3d coords;
	if ( not inWignerSeitz ) {
		Eigen::Vector3d crystalCoords = points.getPointCoords(index,"crystal");
		coords = points.crystalToWS(crystalCoords, basis);
	} else {
		coords = points.getPointCoords(index, basis);
	}
	return coords;
}

template<typename T>
double Point<T>::getWeight() {
	return points.getWeight(index);
}

template<typename T>
bool Point<T>::hasUmklapp() {
	if ( umklappVector.dot(umklappVector) < 1.0e-12 ) {
		return false;
	} else {
		return true;
	}
}

template<typename T>
Point<T> Point<T>::operator + (Point & b) {
	if ( &b.points != &points ) {
		Error e("Points sum should refer to points of the same mesh", 1);
	}
	Eigen::Vector3d coords = getCoords() + b.getCoords();
	long ik = points.getIndex(coords);
	Eigen::Vector3d umklappVector = points.getPointCoords(ik) - coords;
	Point p(ik, umklappVector, points);
    return p;
}

template<typename T>
Point<T> Point<T>::operator - (Point & b) {
	if ( &b.points != &points ) {
		Error e("Points sum should refer to points of the same mesh", 1);
	}
	Eigen::Vector3d coords = getCoords() - b.getCoords();
	long ik = points.getIndex(coords);
	Eigen::Vector3d umklappVector = points.getPointCoords(ik) - coords;
	Point p(ik, umklappVector, points);
    return p;
}

#endif
