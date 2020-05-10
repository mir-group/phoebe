#include <cmath>
#include <set>
#include <iterator>
#include "points.h"
#include "exceptions.h"
#include "eigen.h"
#include "utilities.h" // for mod()

Point::Point(long index_, Eigen::Vector3d umklappVector_, Points & points_)
		: points(points_) {
	umklappVector = umklappVector_;
	index = index_;
}

// copy constructor
Point:: Point( const Point & that ) : umklappVector(that.umklappVector),
		index(that.index), points(that.points) {
}

// copy assignment
Point & Point::operator = ( const Point & that ) {
	if ( this != &that ) {
		umklappVector = that.umklappVector;
		index = that.index;
		points = that.points;
	}
	return *this;
}

Point::~Point() {
}

long Point::getIndex() {
	return index;
}

Eigen::Vector3d Point::getCoords(const std::string & basis,
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

double Point::getWeight() {
	return points.getWeight(index);
}

bool Point::hasUmklapp() {
	if ( umklappVector.dot(umklappVector) < 1.0e-12 ) {
		return false;
	} else {
		return true;
	}
}

Point Point::operator + (Point & b) {
	if ( &b.points != &points ) {
		Error e("Points sum should refer to points of the same mesh", 1);
	}
	Eigen::Vector3d coords = getCoords() + b.getCoords();
	long ik = points.getIndex(coords);
	Eigen::Vector3d umklappVector = points.getPointCoords(ik) - coords;
	Point p(ik, umklappVector, points);
    return p;
}

Point Point::operator - (Point & b) {
	if ( &b.points != &points ) {
		Error e("Points sum should refer to points of the same mesh", 1);
	}
	Eigen::Vector3d coords = getCoords() - b.getCoords();
	long ik = points.getIndex(coords);
	Eigen::Vector3d umklappVector = points.getPointCoords(ik) - coords;
	Point p(ik, umklappVector, points);
    return p;
}

Points::Points(Crystal& crystal_, const Eigen::Vector3i & mesh_,
		const Eigen::Vector3d & offset_,
		const bool useIrreducible_) : crystal{crystal_} {
	setMesh(mesh_, offset_);

	useIrreducible = useIrreducible_;
	if ( useIrreducible ) {
		setIrreduciblePoints();
		useIrreducible = true;
	}

	// This block allows the crystalToWS functionality
	long nGx = 2;
	long nGvec = (2*nGx+1)*(2*nGx+1)*(2*nGx+1);
	Eigen::MatrixXd gVectors_(3,nGvec);
	gVectors = gVectors_;
	gVectors.setZero();
	Eigen::MatrixXi igVectors_(3,nGvec);
	igVectors = igVectors_;
	Eigen::Matrix3d reciprocalUnitCell = crystal.getReciprocalUnitCell();
	igVectors.setZero();
	Eigen::Vector3d vec;
	nGvec = 1; // we skip the first point which is G=(0,0,0)
	for ( long i1=-nGx; i1<=nGx; i1++ ) {
		for ( long i2=-nGx; i2<=nGx; i2++ ) {
			for ( long i3=-nGx; i3<=nGx; i3++ ) {
				if ( ( i1<=0 ) || ( i2<=0 ) || ( i3<=0 ) ) {
					vec(0) = (double)i1;
					vec(1) = (double)i2;
					vec(2) = (double)i3;
					gVectors.col(nGvec) = reciprocalUnitCell * vec;
					igVectors(0,nGvec) = i1;
					igVectors(1,nGvec) = i2;
					igVectors(2,nGvec) = i3;
					nGvec += 1;
				}
			}
		}
	}
}

Points::~Points() {};
FullPoints::~FullPoints() {};

IrreduciblePoints::IrreduciblePoints(Crystal & crystal_,
		const Eigen::Vector3i & mesh_, const Eigen::Vector3d & offset_) :
				Points(crystal_, mesh_, offset_, true) {
}

FullPoints::FullPoints(Crystal & crystal_, const Eigen::Vector3i & mesh_,
		const Eigen::Vector3d & offset_) : Points(crystal_, mesh_,
				offset_, false) {
}

// copy constructor
FullPoints::FullPoints(const FullPoints & that) : Points(that) {
	useIrreducible = useIrreducible;
}

// assignment operator
FullPoints & FullPoints::operator=(const FullPoints & that) {
	if ( this != &that ) {
		mesh = that.mesh;
		offset = that.offset;
		irreduciblePoints = that.irreduciblePoints;
		irreducibleWeights = that.irreducibleWeights;
		mapReducibleToIrreducible = that.mapReducibleToIrreducible;
		mapIrreducibleToReducible = that.mapIrreducibleToReducible;
		indexIrreduciblePoints = that.indexIrreduciblePoints;
		numPoints = that.numPoints;
		numIrredPoints = numIrredPoints;
		useIrreducible = useIrreducible;
		gVectors = that.gVectors;
		igVectors = that.igVectors;
	}
	return *this;
}

Eigen::Vector3d Points::pointsCoords(const long & index) {
	Eigen::Vector3d x;
	if ( useIrreducible ) {
		x = irreduciblePoints.row(index);
	} else {
		x = reduciblePoints(index);
	}
	return x;
}

Point Points::getPoint(const long & index) {
	Eigen::Vector3d p = Eigen::Vector3d::Zero();
	return Point(index, p, *this);
}

Eigen::Vector3d Points::getPointCoords(const long & index,
		const std::string & basis) {
	Eigen::Vector3d pointCrystal = pointsCoords(index);
	if ( basis == "crystal" ) {
		return pointCrystal;
	} else if ( basis == "cartesian" ) {
		Eigen::Vector3d pointCartesian = crystalToCartesian(pointCrystal);
		return pointCartesian;
	} else {
		Error e("Wrong basis for getPoint", 1);
	}
}

std::vector<Eigen::Vector3d> Points::getPointsCoords(const std::string& basis){
	if ( basis != "crystal" && basis != "cartesian" ) {
		Error e("Wrong basis for getPoint", 1);
	}

	std::vector<Eigen::Vector3d> allPoints;

	if (  useIrreducible ) {
		for ( long ik=0; ik<numPoints; ik++ ) {
			if ( basis == "crystal" ) {
				allPoints.push_back(irreduciblePoints.row(ik));
			} else {
				allPoints.push_back(crystalToCartesian(
						irreduciblePoints.row(ik)));
			}
		}
	} else {
		for ( long ik=0; ik<numPoints; ik++ ) {
			if ( basis == "crystal" ) {
				allPoints.push_back(pointsCoords(ik));
			} else {
				allPoints.push_back(crystalToCartesian(pointsCoords(ik)));
			}
		}
	}
	return allPoints;
}

Eigen::Vector3d Points::crystalToWS(const Eigen::Vector3d & pointCrystal,
		const std::string & basis) {

	Eigen::Vector3d pointCart = crystalToCartesian(pointCrystal);

	double norm2 = pointCart.transpose() * pointCart;
	double thisNorm2;
	long iws = 0;
	for ( long iG=1; iG<gVectors.cols(); iG++ ) {
		thisNorm2 = (pointCart+gVectors.col(iG)).transpose()
							  * (pointCart+gVectors.col(iG));
		if ( thisNorm2 < norm2 - 1.0e-12 ) {
			norm2 = thisNorm2;
			iws = iG;
		}
	}

	if ( basis != "crystal" && basis != "cartesian" ) {
		Error e("Wrong input to Wigner Seitz folding", 1);
	}

	if ( basis == "crystal" ) {
		Eigen::Vector3i igVec = igVectors.col(iws);
		Eigen::Vector3d pointCrystalWS = pointCrystal;
		for ( long i=0; i<3; i++) {
			pointCrystalWS(i) += igVec(i);
		}
		return pointCrystalWS;
	} else {
		Eigen::Vector3d pointCartesianWS = pointCart + gVectors.col(iws);
		return pointCartesianWS;
	}

	return pointCart;
}

long Points::getNumPoints() {
	if ( useIrreducible ) {
		return numIrredPoints;
	} else {
		return numPoints;
	}
}

double Points::getWeight(const long & ik) {
	if ( useIrreducible ) {
		return irreducibleWeights(ik);
	} else {
		return 1. / (double)numPoints;
	}
}

long Points::getIndex(const Eigen::Vector3d & point) {
	// given a point coordinate, finds its index in the points list
	// input point must be in crystal coordinates!
	Eigen::Vector3d p;
	// multiply by grid, so that p now contains integers
	p(0) = ( point(0) - offset(0) ) * mesh(0);
	p(1) = ( point(1) - offset(1) ) * mesh(1);
	p(2) = ( point(2) - offset(2) ) * mesh(2);
	// fold in BZ
	long i = mod( long(round(p(0))) , mesh(0));
	long j = mod(long(round(p(1))) , mesh(1));
	long k = mod(long(round(p(2))) , mesh(2));
	long ik = i * mesh(2) * mesh(1) + j * mesh(2) + k;
	return ik;
}

long IrreduciblePoints::getIndex(const Eigen::Vector3d & point) {
	long ik = getIndex(point);
	long ikIrr = mapReducibleToIrreducible(ik);
	return ikIrr;
}

long FullPoints::getIndexInverted(const long & ik) {
	// given the index of point k, return the index of point -k
	Eigen::Vector3d point = reduciblePoints(ik);
	long ikx = (long)round(( point(0) - offset(0) ) * mesh(0));
	long iky = (long)round(( point(1) - offset(1) ) * mesh(1));
	long ikz = (long)round(( point(2) - offset(2) ) * mesh(2));

	long ikxm = - mod(ikx , mesh(0));
	long ikym = - mod(iky , mesh(1));
	long ikzm = - mod(ikz , mesh(2));

	long ikm = ikxm * mesh(2) * mesh(1) + ikym * mesh(2) + ikzm;
	return ikm;
}

void Points::setMesh(const Eigen::Vector3i & mesh_,
		const Eigen::Vector3d & offset_) {

	// validate the mesh and then store it
	if ( mesh_(0) <= 0 ) {
		Error e("meshGrid(0) <= 0, should be positive", 1);
	}
	if ( mesh_(1) <= 0 ) {
		Error e("meshGrid(1) <= 0, should be positive", 1);
	}
	if ( mesh_(2) <= 0 ) {
		Error e("meshGrid(2) <= 0, should be positive", 1);
	}
	mesh = mesh_;
	numPoints = mesh(0) * mesh(1) * mesh(2);

	if ( offset_(0) < 0. && offset_(0) >= 1. ) {
		Error e("offset(0) should be 0 <= offset < 1", 1);
	}
	if ( offset_(1) < 0. && offset_(1) >= 1. ) {
		Error e("offset(1) should be 0 <= offset < 1", 1);
	}
	if ( offset_(2) < 0. && offset_(2) >= 1. ) {
		Error e("offset(2) should be 0 <= offset < 1", 1);
	}
	offset = offset_;

	// I won't set the list of reducible BZ points
	// we generate it on the fly
}

Eigen::Vector3d Points::reduciblePoints(const long & idx) {
	// find 3 indexes from the single index
	long ikz = idx / (mesh(0) * mesh(1));
    long idx_ = idx - (ikz * mesh(0) * mesh(1));
    long iky = idx_ / mesh(0);
    long ikx = mod(idx_ , mesh(0));
	Eigen::Vector3d p;
	p(0) = ikx / (double)mesh(0) + offset(0);
	p(1) = iky / (double)mesh(1) + offset(1);
	p(2) = ikz / (double)mesh(2) + offset(2);
	return p;
}

std::tuple<Eigen::Vector3i, Eigen::Vector3d> Points::getMesh() {
	return {mesh, offset};
}

Eigen::Vector3d Points::crystalToCartesian(const Eigen::Vector3d & point) {
	return crystal.getReciprocalUnitCell() * point;
}

Eigen::Vector3d Points::cartesianToCrystal(const Eigen::Vector3d & point) {
	Eigen::Vector3d p = crystal.getReciprocalUnitCell().inverse() * point;
	return p;
}

void Points::setIrreduciblePoints() {
	// generate a list of points compatible with the symmetries of the system

	const double eps = 1.0e-5;

	Eigen::VectorXd tmpWeight(numPoints);
	Eigen::VectorXi equiv(numPoints);

	// equiv(nk) =nk : k-point nk is not equivalent to any previous k-point
	// equiv(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)

	for ( long nk=0; nk<numPoints; nk++ ) {
		equiv(nk) = nk;
	}

	std::vector<Eigen::Matrix3d> symms = crystal.getSymmetryMatrices();

	Eigen::Vector3d rotatedPoint;
	Eigen::Vector3d thisPoint;
	bool inTheList;
	long ix, iy, iz, n;
	double xx, yy, zz;
	Eigen::Matrix3d s;

	for ( long ik=0; ik<numPoints; ik++ ) {
		// check if this k-point has already been found equivalent to another
		if ( equiv(ik) == ik ) {
			tmpWeight(ik) = 1.;
			// check if there are equivalent k-point to this in the list
			// (excepted those previously found to be equivalent to another)
			// check both k and -k
			for ( auto s : symms ) {
				thisPoint = reduciblePoints(ik);
				rotatedPoint = s * thisPoint;

				for ( long i=0; i<3; i++ ) {
					rotatedPoint(i) -= round( rotatedPoint(i) );
				}
				xx = rotatedPoint(0) * mesh(0) - offset(0);
				yy = rotatedPoint(1) * mesh(1) - offset(1);
				zz = rotatedPoint(2) * mesh(2) - offset(2);
				inTheList = ( abs(xx-round(xx)) <=eps &&
						abs(yy-round(yy)) <=eps &&
						abs(zz-round(zz)) <=eps );
				if ( inTheList ) {
					ix = mod(long(round(rotatedPoint(0)*mesh(0) - offset(0) + 2*mesh(0))), mesh(0) );
					iy = mod(long(round(rotatedPoint(1)*mesh(1) - offset(1) + 2*mesh(1))), mesh(1) );
					iz = mod(long(round(rotatedPoint(2)*mesh(2) - offset(2) + 2*mesh(2))), mesh(2) );

					n = iz + iy*mesh(2) + ix*mesh(1)*mesh(2);
					if ( n>ik && equiv(n)==n ) {
						equiv(n) = ik;
						tmpWeight(ik) += 1.;
					} else {
						if ( equiv(n)!=ik || n<ik ) {
							Error e("Error in finding irred kpoints",1);
						}
					}
				}
			}
		}
	}

	// count irreducible points and order them

	long numIrredPoints_ = 0;
	for ( long ik=0; ik<numPoints; ik++ ) {
		if ( equiv(ik)==ik ) {
			numIrredPoints_ += 1;
		}
	}
	numIrredPoints = numIrredPoints_;

	Eigen::MatrixXd xk(numIrredPoints,3);
	Eigen::VectorXd wk(numIrredPoints);
	Eigen::VectorXi indexIrreduciblePoints_(numIrredPoints);

	long i = 0;
	for ( long j=0; j<numPoints; j++ ) {
		if ( equiv(j)==j ) {
			wk(i) = tmpWeight(j);
			xk.row(i) = reduciblePoints(j);
			indexIrreduciblePoints_(i) = j;
			i += 1;
		}
	}

	// normalize weights to one
	wk /= wk.sum();

	Eigen::VectorXi invEquiv;
	invEquiv.setZero();
	i = 0;
	for ( long ik=0; ik<numPoints; ik++ ) {
		if ( equiv(ik) == ik ) {
			invEquiv(i) = ik;
			i += 1;
		}
	}

	// save as class properties
	irreduciblePoints = xk;
	irreducibleWeights = wk;
	mapReducibleToIrreducible = equiv;
	indexIrreduciblePoints = indexIrreduciblePoints_;
	mapIrreducibleToReducible = invEquiv;
}

Eigen::VectorXi IrreduciblePoints::getIndexReducibleFromIrreducible(
		const long & indexIrr) {
	std::vector<int> indexVec;
	long sizeStar = 0;
	for ( long ik=0; ik<numPoints; ik++ ) {
		if ( mapReducibleToIrreducible(ik) == indexIrr ) {
			indexVec.push_back(ik);
			sizeStar += 1;
		}
	}
	Eigen::VectorXi star(sizeStar);
	for ( long ik=0; ik<sizeStar; ik++ ) {
		star(ik) = indexVec[ik];
	}
	return star;
}

long IrreduciblePoints::getIndexIrreducibleFromReducible(const long& indexRed){
	long ik = mapReducibleToIrreducible(indexRed);
	return ik;
}

std::tuple<Eigen::Vector3i, Eigen::Vector3d> Points::findMesh(
		const Eigen::MatrixXd & testPoints) {
	// given a list of kpoints, figures out the mesh and offset
	// input points must be in crystal coordinates
	Eigen::Vector3i mesh_(3);
	mesh_.setZero();
	Eigen::Vector3d offset_(3);
	offset_.setZero();

	if ( 3 != testPoints.cols() ) {
		Error e("Wrong format specified for points to findMesh",1);
	}
	long numTestPoints = testPoints.rows();

	double value;
	for ( long iCart=0; iCart<3; iCart++ ) {
		std::set<double> s; // note that sets are ordered
		for ( long i = 0; i < numTestPoints; i++ ) {
			// round double to 6 decimal digits
			value = std::ceil(testPoints(i,iCart) * pow(10.,8)) / pow(10.,8);
			// fold number in [0,1[
			while ( value < 0. ) {
				value += 1.;
			}
			while ( value >= 1. ) {
				value -= 1.;
			}
			s.insert( value );
		}

		if ( s.size() > 0 ) {
			mesh_(iCart) = s.size();
			auto first = s.begin(); // get iterator to 1st element of set s
			offset_(iCart) = *first;
		} else {
			Error e("findMesh error: something is wrong", 1);
		}
	}

	if ( numTestPoints != mesh_(0)*mesh_(1)*mesh_(2) ) {
		Error e("Mesh of points seems incomplete", 1);
	}
	return {mesh_, offset_};
}

Crystal & Points::getCrystal() {
	return crystal;
}

ActivePoints::ActivePoints(FullPoints & parentPoints_,
		VectorXl filter) : Points(parentPoints_.getCrystal(),
				std::get<0>(parentPoints_.getMesh()),
				std::get<1>(parentPoints_.getMesh())),
				parentPoints(parentPoints_) {
	// this contain the list of indices of the points in the FullPoints class
	// which we want to include in the ActivePoints class
	filteredToFullIndeces = filter;

	numPoints = filteredToFullIndeces.size();

	// we then construct the list of points
	Eigen::MatrixXd pointsList_(numPoints,3);

	Eigen::Vector3d x;

	long ik;
	for ( long ikNew=0; ikNew<numPoints; ikNew++ ) {
		ik = filteredToFullIndeces(ikNew);
		x = parentPoints.getPointCoords(ik);
		pointsList_.row(ikNew) = x;
	}
	pointsList = pointsList_;
}

ActivePoints::ActivePoints(const ActivePoints & that) : Points(that),
		parentPoints(that.parentPoints) {
	// copy constructor
	pointsList = that.pointsList;
	filteredToFullIndeces = that.filteredToFullIndeces;
}

// copy assignment operator
ActivePoints & ActivePoints::operator=(const ActivePoints & that) {
	if ( this != &that ) {
		parentPoints = that.parentPoints;
		pointsList = that.pointsList;
		filteredToFullIndeces = that.filteredToFullIndeces;
	}
	return *this;
}

Points::Points(const Points & that) : crystal(that.crystal) {
	mesh = that.mesh;
	offset = that.offset;
	irreduciblePoints = that.irreduciblePoints;
	irreducibleWeights = that.irreducibleWeights;
	mapReducibleToIrreducible = that.mapReducibleToIrreducible;
	mapIrreducibleToReducible = that.mapIrreducibleToReducible;
	indexIrreduciblePoints = that.indexIrreduciblePoints;
	numPoints = that.numPoints;
	numIrredPoints = numIrredPoints;
	useIrreducible = useIrreducible;
	gVectors = that.gVectors;
	igVectors = that.igVectors;
}

Points & Points::operator=(const Points & that) { // assignment operator
	if ( this != &that ) {
		mesh = that.mesh;
		offset = that.offset;
		irreduciblePoints = that.irreduciblePoints;
		irreducibleWeights = that.irreducibleWeights;
		mapReducibleToIrreducible = that.mapReducibleToIrreducible;
		mapIrreducibleToReducible = that.mapIrreducibleToReducible;
		indexIrreduciblePoints = that.indexIrreduciblePoints;
		numPoints = that.numPoints;
		numIrredPoints = numIrredPoints;
		useIrreducible = useIrreducible;
		gVectors = that.gVectors;
		igVectors = that.igVectors;
	}
}

Eigen::Vector3d ActivePoints::pointsCoords(const long & index) {
	return pointsList.row(index);
}

long ActivePoints::getIndex(const Eigen::Vector3d & coords) {
	long indexFull = parentPoints.getIndex(coords);
	long ik = fullToFilteredIndeces(indexFull);
	return ik;
}

long ActivePoints::getIndexInverted(const long & ik) {
	long indexFull = filteredToFullIndeces(ik);
	long indexFullInverted = parentPoints.getIndexInverted(indexFull);
	long ikInv = fullToFilteredIndeces(indexFullInverted);
	return ikInv;
}

long ActivePoints::fullToFilteredIndeces(const long & indexIn) {
	// Note: this function could obviously be made much faster if you could
	// save in memory the map of every point in the Full list into the
	// ActivePoints list (and 0 if there's no mapping.
	// But the list of kpoints might be too large!
	for ( long ik=0; ik<numPoints; ik++ ) {
		if ( indexIn == filteredToFullIndeces(ik) ) {
			return ik;
		}
	}
	Error e("Couldn't find the desired kpoint", 1);
}
