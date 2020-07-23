#include <cmath>
#include <set>
#include <math.h>
#include "points.h"
#include "eigen.h"
#include "utilities.h" // for mod()

// default constructors

Points::Points(Crystal &crystal_, const Eigen::Vector3i &mesh_,
               const Eigen::Vector3d &offset_) :
    crystal{crystal_} {

  setMesh(mesh_, offset_);

  // This block allows the folding to the wigner seitz zone
  long nGx = 2;
  long nGvec = (2 * nGx + 1) * (2 * nGx + 1) * (2 * nGx + 1);
  gVectors = Eigen::MatrixXd::Zero(3, nGvec);
  Eigen::Matrix3d reciprocalUnitCell = crystal.getReciprocalUnitCell();
  Eigen::Vector3d vec;
  nGvec = 1; // we skip the first point which is G=(0,0,0)
  for (long i1 = -nGx; i1 <= nGx; i1++) {
    for (long i2 = -nGx; i2 <= nGx; i2++) {
      for (long i3 = -nGx; i3 <= nGx; i3++) {
        if ((i1 <= 0) || (i2 <= 0) || (i3 <= 0)) {
          Eigen::Vector3d vec;
          vec << i1, i2, i3;
          gVectors.col(nGvec) = reciprocalUnitCell * vec;
          nGvec += 1;
        }
      }
    }
  }
}

// copy constructor
Points::Points(const Points &that) :
    crystal(that.crystal), mesh(that.mesh), offset(that.offset),
    numPoints(that.numPoints), gVectors(that.gVectors) {
}

// copy assignment operator
Points &Points::operator=(const Points &that) { // assignment operator
  if (this != &that) {
    crystal = that.crystal;
    mesh = that.mesh;
    offset = that.offset;
    numPoints = that.numPoints;
    gVectors = that.gVectors;
  }
  return *this;
}

Point Points::getPoint(const long &index) {
  return Point(*this, index);
}

// getPointCoords is the tool to find the coordinates of a point

Eigen::Vector3d Points::getPointCoords(const long &index, const int &basis) {
  if (basis != crystalCoords && basis != cartesianCoords) {
    Error e("Wrong basis for getPoint");
  }
  Eigen::Vector3d pointCrystal;
  pointCrystal = reduciblePoints(index);
  if (basis == crystalCoords) {
    return pointCrystal;
  } else {
    Eigen::Vector3d pointCartesian = crystalToCartesian(pointCrystal);
    return pointCartesian;
  }
}

Eigen::Vector3d Points::reduciblePoints(const long &idx) {
  // find 3 indexes from the single index
  long ikz = idx / (mesh(0) * mesh(1));
  long idx_ = idx - (ikz * mesh(0) * mesh(1));
  long iky = idx_ / mesh(0);
  long ikx = mod(idx_, mesh(0));
  Eigen::Vector3d p;
  p(0) = ikx / (double) mesh(0) + offset(0);
  p(1) = iky / (double) mesh(1) + offset(1);
  p(2) = ikz / (double) mesh(2) + offset(2);
  return p;
}

long Points::getIndex(const Eigen::Vector3d &point) {
  // given a point coordinate, finds its index in the points list
  // input point must be in crystal coordinates!
  Eigen::Vector3d p;
  // multiply by grid, so that p now contains integers
  p(0) = (point(0) - offset(0)) * mesh(0);
  p(1) = (point(1) - offset(1)) * mesh(1);
  p(2) = (point(2) - offset(2)) * mesh(2);
  // fold in BZ
  long i = mod(long(round(p(0))), mesh(0));
  long j = mod(long(round(p(1))), mesh(1));
  long k = mod(long(round(p(2))), mesh(2));
  long ik = k * mesh(0) * mesh(1) + j * mesh(0) + i;
  return ik;
}

// change of basis methods

Eigen::Vector3d Points::crystalToCartesian(const Eigen::Vector3d &point) {
  return crystal.getReciprocalUnitCell() * point;
}

Eigen::Vector3d Points::cartesianToCrystal(const Eigen::Vector3d &point) {
  Eigen::Vector3d
      p = crystal.getReciprocalUnitCell().inverse() * point;
  return p;
}

Eigen::Vector3d Points::bzToWs(const Eigen::Vector3d &point,
        const int &basis) {

  Eigen::Vector3d pointCrystal = point;
  if (basis == cartesianCoords) {
    pointCrystal = cartesianToCrystal(pointCrystal);
  }

  // fold to BZ first
  pointCrystal = foldToBz(pointCrystal, crystalCoords);

  Eigen::Vector3d pointCart = crystalToCartesian(pointCrystal);

  double norm2 = pointCart.squaredNorm();
  long iws = 0;
  for (long iG = 1; iG < gVectors.cols(); iG++) {
    double thisNorm2 = (pointCart + gVectors.col(iG)).squaredNorm();
    if (thisNorm2 < norm2 - 1.0e-12) {
      norm2 = thisNorm2;
      iws = iG;
    }
  }

  Eigen::Vector3d point2 = pointCart + gVectors.col(iws);
  if (basis == crystalCoords) {
    point2 = cartesianToCrystal(point2);
  }
  return point2;
}

Eigen::Vector3d Points::foldToBz(const Eigen::Vector3d &point,
        const int &basis) {
  if (basis != crystalCoords && basis != cartesianCoords) {
    Error e("Wrong input to Wigner Seitz folding");
  }

  Eigen::Vector3d pointCrystal = point;
  if ( basis == cartesianCoords) {
    pointCrystal = cartesianToCrystal(pointCrystal);
  }

  for (int i : {0,1,2}) {
    // note: the small double added takes into account for rounding errors
    // which may arise e.g. when the points come from rotations, or machine
    // precision errors in the conversion cartesian to crystal
    // is only a problem in the (unlikely) case of huge k-point meshes (10^8)
    pointCrystal(i) -= std::floor(pointCrystal(i) + 1.0e-8);
  }

  Eigen::Vector3d point2 = pointCrystal;
  if ( basis == cartesianCoords) {
    point2 = crystalToCartesian(point2);
  }
  return point2;
}

Crystal &Points::getCrystal() {
  return crystal;
}

long Points::getNumPoints() {
  return numPoints;
}

double Points::getWeight(const long &ik) {
  (void) ik;
  return 1. / (double) numPoints;
}

void Points::setMesh(const Eigen::Vector3i &mesh_,
                     const Eigen::Vector3d &offset_) {

  // validate the mesh and then store it
  if (mesh_(0) <= 0) {
    Error e("meshGrid(0) <= 0, should be positive", 1);
  }
  if (mesh_(1) <= 0) {
    Error e("meshGrid(1) <= 0, should be positive", 1);
  }
  if (mesh_(2) <= 0) {
    Error e("meshGrid(2) <= 0, should be positive", 1);
  }
  mesh = mesh_;
  numPoints = mesh(0) * mesh(1) * mesh(2);

  if (offset_(0) < 0. && offset_(0) >= 1.) {
    Error e("offset(0) should be 0 <= offset < 1", 1);
  }
  if (offset_(1) < 0. && offset_(1) >= 1.) {
    Error e("offset(1) should be 0 <= offset < 1", 1);
  }
  if (offset_(2) < 0. && offset_(2) >= 1.) {
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
    const Eigen::Matrix<double, 3, Eigen::Dynamic> &testPoints) {
  // given a list of kpoints, figures out the mesh and offset
  // input points must be in crystal coordinates
  Eigen::Vector3i mesh_(3);
  mesh_.setZero();
  Eigen::Vector3d offset_(3);
  offset_.setZero();

  long numTestPoints = testPoints.cols();
  for (long iCart = 0; iCart < 3; iCart++) {
    std::set<double> s; // note that sets are ordered
    for (long i = 0; i < numTestPoints; i++) {
      double value = testPoints(iCart, i);
      // fold number in [0,1[
      while (value < 0.) {
        value += 1.;
      }
      while (value >= 1.) {
        value -= 1.;
      }
      s.insert(value);
    }

    // a few things to remember for comparison.
    // * there might be cases where deltaK is a periodic number (e.g. 1./6)
    // * we are working with floats, so there can be rounding errors
    //   and thus s might contain duplicates

    // determine if there is a single point
    bool isSingle = true;
    for (auto it = s.begin(); it != s.end(); it++) {
      if (it != s.begin()) { // if not in the first case
        // the first condition is the trivial check on having a
        // different point. The second condition is to exclude cases
        // such that s = [0., 0.99999]
        if ((*it - *s.begin() > 1.0e-6)
            && (1. - *it + *s.begin() > 1.0e-6)) {
          isSingle = false;
          break;
        }
      }
    }

    if (isSingle) { // no mesh, just a single point
      mesh_(iCart) = 1;
      offset_(iCart) = *s.begin();
    } else {
      double delta = 0.;
      for (auto it = s.begin(); it != s.end(); it++) {
        if (*it - *s.begin() > 1.0e-6) { // avoid duplicates
          delta = *it - *s.begin();
          break;
        }
      }
      mesh_(iCart) = long(round(1. / delta + 0.1));
      offset_(iCart) = *s.begin();
    }
  }

  if (numTestPoints != mesh_(0) * mesh_(1) * mesh_(2)) {
    Error e("Mesh of points seems incomplete");
  }
  return {mesh_, offset_};
}

std::tuple<int,Eigen::Matrix3d> Points::getRotationToIrreducible(
    const Eigen::Vector3d &x, const int & basis) {
  int ik = getIndex(x);
  Eigen::Matrix3d identity;
  identity.setIdentity();
  (void) basis;
  return {ik, identity};
}

std::vector<Eigen::Matrix3d> Points::getRotationsStar(const int & ik) {
  (void) ik;
  std::vector<Eigen::Matrix3d> rotations;
  Eigen::Matrix3d identity;
  identity.setIdentity();
  rotations.push_back(identity);
  return rotations;
}

int Points::irreducibleToReducible(const int &ikIrr) {
  return ikIrr;
}

////////////////////////////////////////////////////////////////

Point::Point(Points &points_, long index_, Eigen::Vector3d umklappVector_) :
    points(points_) {
  umklappVector = umklappVector_;
  index = index_;
}

// copy constructor
Point::Point(const Point &that) :
    umklappVector(that.umklappVector), index(that.index), points(
    that.points) {
}

// copy assignment
Point &Point::operator=(const Point &that) {
  if (this != &that) {
    umklappVector = that.umklappVector;
    index = that.index;
    points = that.points;
  }
  return *this;
}

long Point::getIndex() {
  return index;
}

Eigen::Vector3d Point::getCoords(const int &basis, const bool &inWignerSeitz) {
  if ((basis != crystalCoords_) && (basis != cartesianCoords_)) {
    Error e("Point getCoordinates: basis must be crystal or cartesian");
  }

  Eigen::Vector3d coords = points.getPointCoords(index, basis);
  if (inWignerSeitz) {
    coords = points.bzToWs(coords, basis);
  }
  return coords;

//  Eigen::Vector3d coords;
//  if (not inWignerSeitz) {
//    Eigen::Vector3d crysCoords = points.getPointCoords(index,
//                                                       crystalCoords_);
//    coords = points.crystalToWS(crysCoords, basis);
//  } else {
//    coords = points.getPointCoords(index, basis);
//  }
//  return coords;
}

double Point::getWeight() {
  return points.getWeight(index);
}

bool Point::hasUmklapp() {
  if (umklappVector.norm() < 1.0e-8) {
    return false;
  } else {
    return true;
  }
}

Point Point::operator+(Point &b) {
  if (&b.points != &points) {
    Error e("Points sum should refer to points of the same mesh");
  }
  Eigen::Vector3d coords = getCoords() + b.getCoords();
  long ik = points.getIndex(coords);
  Eigen::Vector3d umklappVector = points.getPointCoords(ik) - coords;
  Point p(points, ik, umklappVector);
  return p;
}

Point Point::operator-(Point &b) {
  if (&b.points != &points) {
    Error e("Points sum should refer to points of the same mesh");
  }
  Eigen::Vector3d coords = getCoords() - b.getCoords();
  long ik = points.getIndex(coords);
  Eigen::Vector3d umklappVector = points.getPointCoords(ik) - coords;
  Point p(points, ik, umklappVector);
  return p;
}
