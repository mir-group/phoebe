#include <cmath>
#include "points.h"
#include "exceptions.h"

Points::Points(Crystal* crystal_, const Eigen::Vector3i& mesh_) {
	crystal = crystal_;
	setMesh(mesh_);
	setIrreduciblePoints();

	// This block allows the crystalToWS functionality
	int nGx = 2;
	int nGvec = (2*nGx+1)*(2*nGx+1)*(2*nGx+1);
	Eigen::MatrixXd gVectors_(3,nGvec);
	gVectors = gVectors_;
	gVectors.setZero();
	Eigen::MatrixXi igVectors_(3,nGvec);
	igVectors = igVectors_;
	Eigen::Matrix3d reciprocalUnitCell = crystal->getReciprocalUnitCell();
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

Eigen::Vector3d Points::getPoint(int index, std::string basis) {
	Eigen::Vector3d pointCrystal = reduciblePoints.row(index);
	if ( basis == "crystal" ) {
		return pointCrystal;
	} else if ( basis == "cartesian" ) {
		Eigen::Vector3d pointCartesian = crystalToCartesian(pointCrystal);
		return pointCartesian;
	} else {
		Error e("Wrong basis for getPoint", 1);
	}
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
	return numPoints;
}

int Points::getNumIrredPoints() {
	return numIrredPoints;
}

Eigen::MatrixXd Points::getReduciblePoints(const std::string& basis) {
	if ( basis == "crystal" ) {
		return reduciblePoints;
	} else if ( basis == "cartesian" ) {
		Eigen::MatrixXd xk(numPoints,3);
		for ( int ik=0; ik<numPoints; ik++ ) {
			xk.row(ik) = crystalToCartesian(reduciblePoints.row(ik));
		}
		return xk;
	} else {
		Error e("Wrong basis specified in getReduciblePoints", 1);
	}
}

Eigen::MatrixXd Points::getIrreduciblePoints(const std::string& basis) {
	if ( basis == "crystal" ) {
		return irreduciblePoints;
	} else if ( basis == "cartesian" ) {
		Eigen::MatrixXd xk(numIrredPoints,3);
		for ( int ik=0; ik<numIrredPoints; ik++ ) {
			xk.row(ik) = crystalToCartesian(irreduciblePoints.row(ik));
		}
		return xk;

	} else {
		Error e("Wrong basis specified in getIrreduciblePoints", 1);
	}
}

int Points::getIndex(const Eigen::Vector3d& point) {
	// given a point coordinate, finds its index in the points list
	// input point must be in crystal coordinates!

	Eigen::Vector3d p;
	// multiply by grid
	p(0) = point(0) * mesh(0);
	p(1) = point(1) * mesh(1);
	p(2) = point(2) * mesh(2);
	// fold in BZ
	p(0) -= (double)mesh(0) * floor(round(p(0))/(double)mesh(0));
	p(1) -= (double)mesh(1) * floor(round(p(1))/(double)mesh(1));
	p(2) -= (double)mesh(2) * floor(round(p(2))/(double)mesh(2));
	int i = (int)round(p(0));
	int j = (int)round(p(1));
	int k = (int)round(p(2));
	int ik = i * mesh(2) * mesh(1) + j * mesh(2) + k;
	return  ik;
}

int Points::getIndexInverted(int ik) {
	// given the index of point k, return the index of point -k
	Eigen::Vector3d point = reduciblePoints.row(ik);
	int ikx = (int)round(point(0) * mesh(0));
	int iky = (int)round(point(1) * mesh(1));
	int ikz = (int)round(point(2) * mesh(2));

	int ikxm = - ikx;
	int ikym = - iky;
	int ikzm = - ikz;

	ikxm -= (double)mesh(0) * floor((double)ikxm/(double)mesh(0));
	ikym -= (double)mesh(1) * floor((double)ikym/(double)mesh(1));
	ikzm -= (double)mesh(2) * floor((double)ikzm/(double)mesh(2));

	int ikm = ikxm * mesh(2) * mesh(1) + ikym * mesh(2) + ikzm;
	return ikm;
}

void Points::setMesh(Eigen::Vector3i mesh_) {

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

	// set the list of reducible BZ points
	Eigen::MatrixXd reduciblePoints_(numPoints,3);
	int ik = 0;
	for ( int ikx=0; ikx<mesh(0); ikx++ ) {
		for ( int iky=0; iky<mesh(1); iky++ ) {
			for ( int ikz=0; ikz<mesh(2); ikz++ ) {
				reduciblePoints_(ik,0) = ikx / (double)mesh(0);
				reduciblePoints_(ik,1) = iky / (double)mesh(1);
				reduciblePoints_(ik,2) = ikz / (double)mesh(2);
				ik += 1;
			}
		}
	}
	reduciblePoints = reduciblePoints_;
}

Eigen::Vector3i Points::getMesh() {
	return mesh;
}

Eigen::Vector3d Points::crystalToCartesian(const Eigen::Vector3d& point) {
	Eigen::Vector3d p = crystal->getReciprocalUnitCell() * point;
	return p;
}

Eigen::Vector3d Points::cartesianToCrystal(const Eigen::Vector3d& point) {
	Eigen::Vector3d p = point.transpose() * crystal->getReciprocalUnitCell();
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

	std::vector<Eigen::Matrix3d> symms = crystal->getSymmetryMatrices();

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
			for ( int is=0; is<symms.size(); is++ ) {
				n=0;
				s = symms[is];
				thisPoint = reduciblePoints.row(ik);
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

	int i = 0;
	for ( int j=0; j<numPoints; j++ ) {
		if ( equiv(j)==j ) {
			wk(i) = tmpWeight(j);
			xk.row(i) = reduciblePoints.row(j);
			i += 1;
		}
	}

	// normalize weights to one
	wk /= wk.sum();

	// save as class properties
	reduciblePoints = xk;
	reducibleWeights = wk;
}
