#include <cmath>
#include "points.h"
#include "exceptions.h"
#include <set>

Eigen::Vector3d crystalToCartesian(const Eigen::Vector3d& point,
		const Eigen::Matrix3d& reciprocalUnitCell) {
	return reciprocalUnitCell * point;
}

Point::Point(Eigen::Vector3d crystalCoords_, Eigen::Vector3d crystalCoordsWS_,
		const Eigen::Matrix3d& reciprocalUnitCell_)
		: reciprocalUnitCell(reciprocalUnitCell_) {
	crystalCoords = crystalCoords_;
	crystalCoordsWS = crystalCoordsWS_;
}

Eigen::Vector3d Point::getCoords(std::string basis, bool inWignerSeitz) {
	if ( (basis != "crystal") && (basis != "cartesian") ) {
		Error e("Point getCoordinates: basis must be crystal or cartesian", 1);
	}

	if ( not inWignerSeitz ) {
		if ( basis == "crystal" ) {
			return crystalCoords;
		} else {
			return crystalToCartesian(crystalCoords, reciprocalUnitCell);
		}
	} else {
		if ( basis == "crystal" ) {
			return crystalCoordsWS;
		} else {
			return crystalToCartesian(crystalCoordsWS, reciprocalUnitCell);
		}
	}
}

double Point::getWeight() {
	return weight;
}

Points::Points(Crystal& crystal_, const Eigen::Vector3i& mesh_,
		const Eigen::Vector3d& offset_,
		const bool useIrreducible_) : crystal{crystal_} {
	setMesh(mesh_, offset_);

	useIrreducible = useIrreducible_;
	if ( useIrreducible ) {
		setIrreduciblePoints();
		useIrreducible = true;
	}

	// This block allows the crystalToWS functionality
	int nGx = 2;
	int nGvec = (2*nGx+1)*(2*nGx+1)*(2*nGx+1);
	Eigen::MatrixXd gVectors_(3,nGvec);
	gVectors = gVectors_;
	gVectors.setZero();
	Eigen::MatrixXi igVectors_(3,nGvec);
	igVectors = igVectors_;
	Eigen::Matrix3d reciprocalUnitCell = crystal.getReciprocalUnitCell();
	igVectors.setZero();
	Eigen::Vector3d vec;
	nGvec = 1; // we skip the first point which is G=(0,0,0)
	for ( int i1=-nGx; i1<=nGx; i1++ ) {
		for ( int i2=-nGx; i2<=nGx; i2++ ) {
			for ( int i3=-nGx; i3<=nGx; i3++ ) {
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

IrreduciblePoints::IrreduciblePoints(Crystal& crystal_,
		const Eigen::Vector3i& mesh_,
		const Eigen::Vector3d& offset_) : Points(crystal_, mesh_,
				offset_, true) {
}

FullPoints::FullPoints(Crystal& crystal_, const Eigen::Vector3i& mesh_,
		const Eigen::Vector3d& offset_) : Points(crystal_, mesh_,
				offset_, false) {
}

Eigen::Vector3d Points::pointsCoords(const int& index){
	Error e("Base Points class doesn't have pointsCoords implemented", 1);
	Eigen::Vector3d x;
	return x.setZero();
}

Eigen::Vector3d IrreduciblePoints::pointsCoords(const int& index) {
	Eigen::Vector3d x = irreduciblePoints.row(index);
	return x;
}

Eigen::Vector3d FullPoints::pointsCoords(const int& index) {
	Eigen::Vector3d x = reduciblePoints(index);
	return x;
}

Point Points::getPoint(const int& index) {
	Eigen::Vector3d p = pointsCoords(index);
	return Point(p,crystalToWS(p,"crystal"),crystal.getReciprocalUnitCell());
}

Eigen::Vector3d Points::getPointCoords(int& index, std::string basis) {
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

std::vector<Eigen::Vector3d> Points::getPointsCoords(std::string basis) {
	if ( basis != "crystal" && basis != "cartesian" ) {
		Error e("Wrong basis for getPoint", 1);
	}

	std::vector<Eigen::Vector3d> allPoints;

	if (  useIrreducible ) {
		for ( int ik=0; ik<numPoints; ik++ ) {
			if ( basis == "crystal" ) {
				allPoints.push_back(irreduciblePoints.row(ik));
			} else {
				allPoints.push_back(crystalToCartesian(irreduciblePoints.row(ik)));
			}
		}
	} else {
		for ( int ik=0; ik<numPoints; ik++ ) {
			if ( basis == "crystal" ) {
				allPoints.push_back(pointsCoords(ik));
			} else {
				allPoints.push_back(crystalToCartesian(pointsCoords(ik)));
			}
		}
	}
	return allPoints;
}

Eigen::Vector3d Points::crystalToWS(const Eigen::Vector3d& pointCrystal,
		const std::string& basis) {

	Eigen::Vector3d pointCart = crystalToCartesian(pointCrystal);

	double norm2 = pointCart.transpose() * pointCart;
	double thisNorm2;
	int iws = 0;
	for ( int iG=1; iG<gVectors.cols(); iG++ ) {
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
		for ( int i=0; i<3; i++) {
			pointCrystalWS(i) += igVec(i);
		}
		return pointCrystalWS;
	} else {
		Eigen::Vector3d pointCartesianWS = pointCart + gVectors.col(iws);
		return pointCartesianWS;
	}

	return pointCart;
}

int Points::getNumPoints() {
	if ( useIrreducible ) {
		return numIrredPoints;
	} else {
		return numPoints;
	}
}

double Points::getWeight(int ik) {
	if ( useIrreducible ) {
		return irreducibleWeights(ik);
	} else {
		return 1. / (double)numPoints;
	}
}

int Points::getIndex(const Eigen::Vector3d& point) {
	// given a point coordinate, finds its index in the points list
	// input point must be in crystal coordinates!
	Eigen::Vector3d p;
	// multiply by grid
	p(0) = ( point(0) - offset(0) ) * mesh(0);
	p(1) = ( point(1) - offset(1) ) * mesh(1);
	p(2) = ( point(2) - offset(2) ) * mesh(2);
	// fold in BZ
	p(0) -= (double)mesh(0) * floor(round(p(0))/(double)mesh(0));
	p(1) -= (double)mesh(1) * floor(round(p(1))/(double)mesh(1));
	p(2) -= (double)mesh(2) * floor(round(p(2))/(double)mesh(2));
	int i = (int)round(p(0));
	int j = (int)round(p(1));
	int k = (int)round(p(2));
	int ik = i * mesh(2) * mesh(1) + j * mesh(2) + k;
	return ik;
}

int IrreduciblePoints::getIndex(const Eigen::Vector3d& point) {
	int ik = getIndex(point);
	int ikIrr = mapReducibleToIrreducible(ik);
	return ikIrr;
}

int FullPoints::getIndexInverted(const int& ik) {
	// given the index of point k, return the index of point -k
	Eigen::Vector3d point = reduciblePoints(ik);
	int ikx = (int)round(( point(0) - offset(0) ) * mesh(0));
	int iky = (int)round(( point(1) - offset(1) ) * mesh(1));
	int ikz = (int)round(( point(2) - offset(2) ) * mesh(2));

	int ikxm = - ikx;
	int ikym = - iky;
	int ikzm = - ikz;

	ikxm -= (double)mesh(0) * floor((double)ikxm/(double)mesh(0));
	ikym -= (double)mesh(1) * floor((double)ikym/(double)mesh(1));
	ikzm -= (double)mesh(2) * floor((double)ikzm/(double)mesh(2));

	int ikm = ikxm * mesh(2) * mesh(1) + ikym * mesh(2) + ikzm;
	return ikm;
}

void Points::setMesh(Eigen::Vector3i mesh_, Eigen::Vector3d offset_) {

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

Eigen::Vector3d Points::reduciblePoints(const int idx) {
	// find 3 indexes from the single index
	int ikz = idx / (mesh(0) * mesh(1));
    int idx_ = idx - (ikz * mesh(0) * mesh(1));
    int iky = idx_ / mesh(0);
    int ikx = idx_ % mesh(0);
	Eigen::Vector3d p;
	p(0) = ikx / (double)mesh(0) + offset(0);
	p(1) = iky / (double)mesh(1) + offset(1);
	p(2) = ikz / (double)mesh(2) + offset(2);
	return p;
}

std::tuple<Eigen::Vector3i, Eigen::Vector3d> Points::getMesh() {
	return {mesh, offset};
}

Eigen::Vector3d Points::crystalToCartesian(const Eigen::Vector3d& point) {
	return crystal.getReciprocalUnitCell() * point;
}

Eigen::Vector3d Points::cartesianToCrystal(const Eigen::Vector3d& point) {
	Eigen::Vector3d p = point.transpose() * crystal.getReciprocalUnitCell();
	return p;
}

void Points::setIrreduciblePoints() {
	// generate a list of points compatible with the symmetries of the system

	const double eps = 1.0e-5;

	Eigen::VectorXd tmpWeight(numPoints);
	Eigen::VectorXi equiv(numPoints);

	// equiv(nk) =nk : k-point nk is not equivalent to any previous k-point
	// equiv(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)

	for ( int nk=0; nk<numPoints; nk++ ) {
		equiv(nk) = nk;
	}

	std::vector<Eigen::Matrix3d> symms = crystal.getSymmetryMatrices();

	Eigen::Vector3d rotatedPoint;
	Eigen::Vector3d thisPoint;
	bool inTheList;
	int ix, iy, iz, n;
	double xx, yy, zz;
	Eigen::Matrix3d s;

	for ( int ik=0; ik<numPoints; ik++ ) {
		// check if this k-point has already been found equivalent to another
		if ( equiv(ik) == ik ) {
			tmpWeight(ik) = 1.;
			// check if there are equivalent k-point to this in the list
			// (excepted those previously found to be equivalent to another)
			// check both k and -k
			for ( auto s : symms ) {
				thisPoint = reduciblePoints(ik);
				rotatedPoint = s * thisPoint;

				for ( int i=0; i<3; i++ ) {
					rotatedPoint(i) -= round( rotatedPoint(i) );
				}
				xx = rotatedPoint(0) * mesh(0);
				yy = rotatedPoint(1) * mesh(1);
				zz = rotatedPoint(2) * mesh(2);
				inTheList = ( abs(xx-round(xx)) <=eps &&
						abs(yy-round(yy)) <=eps &&
						abs(zz-round(zz)) <=eps );
				if ( inTheList ) {
					ix = (int)round( ( rotatedPoint(0)+2 )*mesh(0) ) % mesh(0);
					iy = (int)round( ( rotatedPoint(1)+2 )*mesh(1) ) % mesh(1);
					iz = (int)round( ( rotatedPoint(2)+2 )*mesh(2) ) % mesh(2);
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

	int numIrredPoints_ = 0;
	for ( int ik=0; ik<numPoints; ik++ ) {
		if ( equiv(ik)==ik ) {
			numIrredPoints_ += 1;
		}
	}
	numIrredPoints = numIrredPoints_;

	Eigen::MatrixXd xk(numIrredPoints,3);
	Eigen::VectorXd wk(numIrredPoints);
	Eigen::VectorXi indexIrreduciblePoints_(numIrredPoints);

	int i = 0;
	for ( int j=0; j<numPoints; j++ ) {
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
	for ( int ik=0; ik<numPoints; ik++ ) {
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

Eigen::VectorXi IrreduciblePoints::getIndexReducibleFromIrreducible(int indexIrr) {
	std::vector<int> indexVec;
	int sizeStar = 0;
	for ( int ik=0; ik<numPoints; ik++ ) {
		if ( mapReducibleToIrreducible(ik) == indexIrr ) {
			indexVec.push_back(ik);
			sizeStar += 1;
		}
	}
	Eigen::VectorXi star(sizeStar);
	for ( int ik=0; ik<sizeStar; ik++ ) {
		star(ik) = indexVec[ik];
	}
	return star;
}

int IrreduciblePoints::getIndexIrreducibleFromReducible(int indexRed) {
	int ik = mapReducibleToIrreducible(indexRed);
	return ik;
}

std::tuple<Eigen::Vector3i, Eigen::Vector3d> Points::findMesh(
		Eigen::MatrixXd& testPoints) {
	// given a list of kpoints, figures out the mesh and offset
	// input points must be in crystal coordinates
	Eigen::Vector3i mesh_(3);
	mesh_.setZero();
	Eigen::Vector3d offset_(3);
	offset_.setZero();

	if ( 3 != testPoints.cols() ) {
		Error e("Wrong format specified for points to findMesh",1);
	}
	int numTestPoints = testPoints.rows();

	double value;
	for ( int iCart=0; iCart<3; iCart++ ) {
		std::set<double> s; // note that sets are ordered
		for ( int i = 0; i < numTestPoints; i++ ) {
			// round double to 6 decimal digits
			value = std::ceil(testPoints(i,iCart) * 100000000.0) / 100000000.0;
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
