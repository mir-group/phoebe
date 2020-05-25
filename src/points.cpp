#include <cmath>
#include <set>
#include <iterator>
#include "points.h"
#include "eigen.h"
#include "constants.h"
#include "utilities.h" // for mod()

// default constructors

Points::Points(Crystal& crystal_, const Eigen::Vector3i & mesh_,
		const Eigen::Vector3d & offset_) : crystal{crystal_} {

	setMesh(mesh_, offset_);

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

FullPoints::FullPoints(Crystal & crystal_, const Eigen::Vector3i & mesh_,
		const Eigen::Vector3d & offset_) : Points(crystal_, mesh_, offset_) {
}

IrreduciblePoints::IrreduciblePoints(Crystal & crystal_,
		const Eigen::Vector3i & mesh_, const Eigen::Vector3d & offset_) :
				Points(crystal_, mesh_, offset_) {
	setIrreduciblePoints();
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
	Eigen::MatrixXd pointsList_(3,numPoints);

	Eigen::Vector3d x;

	long ik;
	for ( long ikNew=0; ikNew<numPoints; ikNew++ ) {
		ik = filteredToFullIndeces(ikNew);
		x = parentPoints.getPointCoords(ik);
		pointsList_.col(ikNew) = x;
	}
	pointsList = pointsList_;
}

PathPoints::PathPoints(Crystal & crystal_,
		const Eigen::Tensor<double,3> & pathExtrema, const double & delta) :
		FullPoints(crystal_, Eigen::Vector3i::Constant(1),
				Eigen::Vector3d::Zero()){

	// build the list of points
	std::vector<Eigen::Vector3d> points;
	Eigen::Vector3d p0 , p1;//, thisDelta;

	// initialize the path
	p0(0) = pathExtrema(0,0,0);
	p0(1) = pathExtrema(0,0,1);
	p0(2) = pathExtrema(0,0,2);
	points.push_back(p0);

	// we loop over the segments provided in user input
	for ( long i=0; i<pathExtrema.dimension(0); i++ ) {
		// load coords of the extrema of the segment
		p0(0) = pathExtrema(i,0,0);
		p0(1) = pathExtrema(i,0,1);
		p0(2) = pathExtrema(i,0,2);
		p1(0) = pathExtrema(i,1,0);
		p1(1) = pathExtrema(i,1,1);
		p1(2) = pathExtrema(i,1,2);

		// delta may not divide the interval exactly
		// so, we find the closest one
		long nk = abs(long( (p1-p0).norm() / delta ));

		// now we build the points of the segment
		std::vector<Eigen::Vector3d> segmentPoints;
		for ( long j=0; j<=nk; j++ ) {
			Eigen::Vector3d thisP;
			thisP(0) = (p1(0) - p0(0)) / nk * j + p0(0);
			thisP(1) = (p1(1) - p0(1)) / nk * j + p0(1);
			thisP(2) = (p1(2) - p0(2)) / nk * j + p0(2);
			segmentPoints.push_back(thisP);
		}

		// now we update the original list of points
		// but we have to check to not add the same extrema twice
		// work with the first element
		Eigen::Vector3d lastP = points[points.size()-1];
		Eigen::Vector3d thisP = segmentPoints[0];
		bool addIt = (thisP - lastP).norm() > epsilon8 ;
		if ( addIt ) {
			points.push_back(thisP);
		}

		// select first element of the segment
		auto p = std::begin(segmentPoints);
		++p;
		// and then add all the rest of the segment
		for ( auto end=std::end(segmentPoints); p!=end; ++p ) {
			// iterate over the rest of the container
			points.push_back(*p);
		}
	}

	numPoints = points.size();
	pointsList = Eigen::MatrixXd::Zero(3,numPoints);
	long i = 0;
	for ( auto p : points ) {
		pointsList.col(i) = p;
		i += 1;
	}
}

// copy constructors

// copy constructor
Points::Points(const Points & that) : crystal(that.crystal), mesh(that.mesh),
	offset(that.offset), numPoints(that.numPoints), gVectors(that.gVectors),
	igVectors(that.igVectors) {
}

// copy constructor
FullPoints::FullPoints(const FullPoints & that) : Points(that) {
}

// copy constructor
IrreduciblePoints::IrreduciblePoints(const IrreduciblePoints & that) :
		Points(that),
		mapReducibleToIrreducible(that.mapReducibleToIrreducible),
		mapIrreducibleToReducible(that.mapIrreducibleToReducible),
		irreduciblePoints(that.irreduciblePoints),
		irreducibleWeights(that.irreducibleWeights),
		indexIrreduciblePoints(that.indexIrreduciblePoints),
		numIrredPoints(that.numIrredPoints) {
}

// copy constructor
ActivePoints::ActivePoints(const ActivePoints & that) : Points(that),
		parentPoints(that.parentPoints), pointsList(that.pointsList),
		filteredToFullIndeces(that.filteredToFullIndeces) {
}

// copy constructor
PathPoints::PathPoints(const PathPoints & that) : FullPoints(that),
		pointsList(that.pointsList) {
}

// copy assignment operator

Points & Points::operator=(const Points & that) { // assignment operator
	if ( this != &that ) {
		crystal = that.crystal;
		mesh = that.mesh;
		offset = that.offset;
		numPoints = that.numPoints;
		gVectors = that.gVectors;
		igVectors = that.igVectors;
	}
	return *this;
}

// assignment operator
FullPoints & FullPoints::operator=(const FullPoints & that) {
	if ( this != &that ) {
		crystal = that.crystal;
		mesh = that.mesh;
		offset = that.offset;
		numPoints = that.numPoints;
		gVectors = that.gVectors;
		igVectors = that.igVectors;
	}
	return *this;
}

// assignment operator
IrreduciblePoints & IrreduciblePoints::operator=(
		const IrreduciblePoints & that) {
	if ( this != &that ) {
		crystal = that.crystal;
		mesh = that.mesh;
		offset = that.offset;
		numPoints = that.numPoints;
		gVectors = that.gVectors;
		igVectors = that.igVectors;
		mapReducibleToIrreducible = that.mapReducibleToIrreducible;
		mapIrreducibleToReducible = that.mapIrreducibleToReducible;
		irreduciblePoints = that.irreduciblePoints;
		irreducibleWeights = that.irreducibleWeights;
		indexIrreduciblePoints = that.indexIrreduciblePoints;
		numIrredPoints = that.numIrredPoints;
	}
	return *this;
}

// copy assignment operator
ActivePoints & ActivePoints::operator=(const ActivePoints & that) {
	if ( this != &that ) {
		crystal = that.crystal;
		mesh = that.mesh;
		offset = that.offset;
		numPoints = that.numPoints;
		gVectors = that.gVectors;
		igVectors = that.igVectors;
		parentPoints = that.parentPoints;
		filteredToFullIndeces = that.filteredToFullIndeces;
		pointsList = that.pointsList;
	}
	return *this;
}

// assignment operator
PathPoints & PathPoints::operator=(const PathPoints & that) {
	if ( this != &that ) {
		crystal = that.crystal;
		mesh = that.mesh;
		offset = that.offset;
		numPoints = that.numPoints;
		gVectors = that.gVectors;
		igVectors = that.igVectors;
		pointsList = that.pointsList;
	}
	return *this;
}

// getPoint methods

//Point<Points> Points::getPoint(const long & index) {
//	Eigen::Vector3d p = getPointCoords(index);
//	return Point<Points>(index, p, *this);
//}

Point<FullPoints> FullPoints::getPoint(const long & index) {
	Eigen::Vector3d p = getPointCoords(index);
	return Point<FullPoints>(index, p, *this);
}

Point<IrreduciblePoints> IrreduciblePoints::getPoint(const long & index) {
	Eigen::Vector3d p = getPointCoords(index);
	return Point<IrreduciblePoints>(index, p, *this);
}

Point<ActivePoints> ActivePoints::getPoint(const long & index) {
	Eigen::Vector3d p = getPointCoords(index);
	return Point<ActivePoints>(index, p, *this);
}

Point<PathPoints> PathPoints::getPoint(const long & index) {
	Eigen::Vector3d p = getPointCoords(index);
	return Point<PathPoints>(index, p, *this);
}

// getPointCoords is the tool to find the coordinates of a point

Eigen::Vector3d Points::getPointCoords(const long & index,
		const std::string & basis) {
	if ( basis != "crystal" && basis != "cartesian" ) {
		Error e("Wrong basis for getPoint", 1);
	}
	Eigen::Vector3d pointCrystal;
	pointCrystal = reduciblePoints(index);
	if ( basis == "crystal" ) {
		return pointCrystal;
	} else {
		Eigen::Vector3d pointCartesian = crystalToCartesian(pointCrystal);
		return pointCartesian;
	}
}

//FullPoints just like Points

Eigen::Vector3d IrreduciblePoints::getPointCoords(const long & index,
		const std::string & basis) {
	if ( basis != "crystal" && basis != "cartesian" ) {
		Error e("Wrong basis for getPoint", 1);
	}
	Eigen::Vector3d pointCrystal = irreduciblePoints.col(index);
	if ( basis == "crystal" ) {
		return pointCrystal;
	} else {
		Eigen::Vector3d pointCartesian = crystalToCartesian(pointCrystal);
		return pointCartesian;
	}
}
Eigen::Vector3d ActivePoints::getPointCoords(const long & index,
		const std::string & basis) {
	if ( basis != "crystal" && basis != "cartesian" ) {
		Error e("Wrong basis for getPoint", 1);
	}
	Eigen::Vector3d pointCrystal = pointsList.col(index);
	if ( basis == "crystal" ) {
		return pointCrystal;
	} else {
		Eigen::Vector3d pointCartesian = crystalToCartesian(pointCrystal);
		return pointCartesian;
	}
}
Eigen::Vector3d PathPoints::getPointCoords(const long & index,
		const std::string & basis) {
	if ( basis != "crystal" && basis != "cartesian" ) {
		Error e("Wrong basis for getPoint", 1);
	}
	Eigen::Vector3d pointCrystal = pointsList.col(index);
	if ( basis == "crystal" ) {
		return pointCrystal;
	} else {
		Eigen::Vector3d pointCartesian = crystalToCartesian(pointCrystal);
		return pointCartesian;
	}
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

// getIndex methods

long Points::getIndex(const Eigen::Vector3d & point) {
	// given a point coordinate, finds its index in the points list
	// input point must be in crystal coordinates!
	Eigen::Vector3d p;
	// multiply by grid, so that p now contains integers
	p(0) = ( point(0) - offset(0) ) * mesh(0);
	p(1) = ( point(1) - offset(1) ) * mesh(1);
	p(2) = ( point(2) - offset(2) ) * mesh(2);
	// fold in BZ
	long i = mod(long(round(p(0))) , mesh(0));
	long j = mod(long(round(p(1))) , mesh(1));
	long k = mod(long(round(p(2))) , mesh(2));
	long ik = i * mesh(2) * mesh(1) + j * mesh(2) + k;
	return ik;
}

// FullPoints::getIndex is just like getIndex

long IrreduciblePoints::getIndex(const Eigen::Vector3d & point) {
	// untested
	long ik = getIndex(point);
	long ikIrr = mapReducibleToIrreducible(ik);
	return ikIrr;
}

long ActivePoints::getIndex(const Eigen::Vector3d & coords) {
	// we take advantage of the fact that the parent points have an order
	long indexFull = parentPoints.getIndex(coords);
	long ik = fullToFilteredIndeces(indexFull);
	return ik;
}

long PathPoints::getIndex(const Eigen::Vector3d & coords) {
	// in this case there is no order, so we just search through a loop
	long counter = 0;
	for ( counter=0; counter<numPoints; counter++ ) {
		if ( (pointsList.col(counter)-coords).norm() < 1.0e-8 ) {
			break;
		}
	}
	return counter;
}

// methods to find the index of the point -k, given k

long FullPoints::getIndexInverted(const long & ik) {
	// given the index of point k, return the index of point -k
	Eigen::Vector3d point = reduciblePoints(ik);
	long ikx = (long)round(( point(0) - offset(0) ) * mesh(0));
	long iky = (long)round(( point(1) - offset(1) ) * mesh(1));
	long ikz = (long)round(( point(2) - offset(2) ) * mesh(2));

	long ikxm = mod(-ikx , mesh(0));
	long ikym = mod(-iky , mesh(1));
	long ikzm = mod(-ikz , mesh(2));

	long ikm = ikxm * mesh(2) * mesh(1) + ikym * mesh(2) + ikzm;
	return ikm;
}

long ActivePoints::getIndexInverted(const long & ik) {
	long indexFull = filteredToFullIndeces(ik);
	long indexFullInverted = parentPoints.getIndexInverted(indexFull);
	long ikInv = fullToFilteredIndeces(indexFullInverted);
	return ikInv;
}

// change of basis methods

Eigen::Vector3d Points::crystalToCartesian(const Eigen::Vector3d & point) {
	return crystal.getReciprocalUnitCell() * point;
}

Eigen::Vector3d Points::cartesianToCrystal(const Eigen::Vector3d & point) {
	Eigen::Vector3d p = crystal.getReciprocalUnitCell().inverse() * point;
	return p;
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

Crystal & Points::getCrystal() {
	return crystal;
}

long Points::getNumPoints() {
	return numPoints;
}

long IrreduciblePoints::getNumPoints() {
	return numIrredPoints;
}

double Points::getWeight(const long & ik) {
	(void) ik;
	return 1. / (double)numPoints;
}

double IrreduciblePoints::getWeight(const long & ik) {
	return irreducibleWeights(ik);
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

std::tuple<Eigen::Vector3i, Eigen::Vector3d> Points::getMesh() {
	return {mesh, offset};
}

std::tuple<Eigen::Vector3i, Eigen::Vector3d> Points::findMesh(
		const Eigen::Matrix<double,3,Eigen::Dynamic> & testPoints) {
	// given a list of kpoints, figures out the mesh and offset
	// input points must be in crystal coordinates
	Eigen::Vector3i mesh_(3);
	mesh_.setZero();
	Eigen::Vector3d offset_(3);
	offset_.setZero();

	long numTestPoints = testPoints.cols();
	for ( long iCart=0; iCart<3; iCart++ ) {
		std::set<double> s; // note that sets are ordered
		for ( long i = 0; i < numTestPoints; i++ ) {
			double value = testPoints(iCart,i);
			// fold number in [0,1[
			while ( value < 0. ) {
				value += 1.;
			}
			while ( value >= 1. ) {
				value -= 1.;
			}
			s.insert( value );
		}

		// a few things to remember for comparison.
		// * there might be cases where deltaK is a periodic number (e.g. 1./6)
		// * we are working with floats, so there can be rounding errors
		//   and thus s might contain duplicates

		// determine if there is a single point
		bool isSingle = true;
		for ( auto it=s.begin(); it!= s.end(); it++ ) {
			if ( it != s.begin() ) { // if not in the first case
				// the first condition is the trivial check on having a
				// different point. The second condition is to exclude cases
				// such that s = [0., 0.99999]
				if ( (*it - *s.begin() > 1.0e-6) &&
						(1.-*it + *s.begin() > 1.0e-6) ) {
					isSingle = false;
					break;
				}
			}
		}

		if ( isSingle ) { // no mesh, just a single point
			mesh_(iCart) = 1;
			offset_(iCart) = *s.begin();
		} else {
			double delta = 0.;
			for ( auto it=s.begin(); it!= s.end(); it++ ) {
				if ( *it - *s.begin() > 1.0e-6 ) { // avoid duplicates
					delta = *it - *s.begin();
					break;
				}
			}
			mesh_(iCart) = long(round(1./delta+0.1));
			offset_(iCart) = *s.begin();
		}
	}

	if ( numTestPoints != mesh_(0)*mesh_(1)*mesh_(2) ) {
		Error e("Mesh of points seems incomplete", 1);
	}
	return {mesh_, offset_};
}

void IrreduciblePoints::setIrreduciblePoints() {
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

	long numIrredPoints = 0;
	for ( long ik=0; ik<numPoints; ik++ ) {
		if ( equiv(ik)==ik ) {
			numIrredPoints += 1;
		}
	}

	Eigen::MatrixXd xk(3,numIrredPoints);
	Eigen::VectorXd wk(numIrredPoints);
	Eigen::VectorXi indexIrreduciblePoints_(numIrredPoints);

	long i = 0;
	for ( long j=0; j<numPoints; j++ ) {
		if ( equiv(j)==j ) {
			wk(i) = tmpWeight(j);
			xk.col(i) = reduciblePoints(j);
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

// maps between two internal lists of points

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

long ActivePoints::fullToFilteredIndeces(const long & indexIn) {
	// Note: this function could obviously be made much faster if you could
	// save in memory the map of every point in the Full list into the
	// ActivePoints list (and 0 if there's no mapping.
	// But the list of kpoints might be too large!
	long target = -1;
	for ( long ik=0; ik<numPoints; ik++ ) {
		if ( indexIn == filteredToFullIndeces(ik) ) {
			target = ik;
			break;
		}
	}
	if ( target == -1 ) {
		Error e("Couldn't find the desired kpoint");
	}
	return target;
}
