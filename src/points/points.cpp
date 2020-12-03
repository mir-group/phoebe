#include "points.h"
#include "eigen.h"
#include "mpiHelper.h"
#include "utilities.h" // for mod()
#include <cmath>
#include <math.h>
#include <set>

const int Points::crystalCoords = 0;
const int Points::cartesianCoords = 1;

// default constructors

Points::Points(Crystal &crystal_, const Eigen::Vector3i &mesh_,
               const Eigen::Vector3d &offset_)
    : crystal{crystal_} {

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
Points::Points(const Points &that)
    : crystal(that.crystal), mesh(that.mesh), offset(that.offset),
      numPoints(that.numPoints), gVectors(that.gVectors),
      rotationMatricesCrystal(that.rotationMatricesCrystal),
      rotationMatricesCartesian(that.rotationMatricesCartesian),
      mapEquivalenceRotationIndex(that.mapEquivalenceRotationIndex),
      mapIrreducibleToReducibleList(that.mapIrreducibleToReducibleList),
      mapReducibleToIrreducibleList(that.mapReducibleToIrreducibleList),
      numIrrPoints(that.numIrrPoints), irreducibleStars(that.irreducibleStars),
      equiv(that.equiv) {}

// copy assignment operator
Points &Points::operator=(const Points &that) { // assignment operator
  if (this != &that) {
    crystal = that.crystal;
    mesh = that.mesh;
    offset = that.offset;
    numPoints = that.numPoints;
    gVectors = that.gVectors;
    rotationMatricesCrystal = that.rotationMatricesCrystal;
    rotationMatricesCartesian = that.rotationMatricesCartesian;
    mapEquivalenceRotationIndex = that.mapEquivalenceRotationIndex;
    mapIrreducibleToReducibleList = that.mapIrreducibleToReducibleList;
    mapReducibleToIrreducibleList = that.mapReducibleToIrreducibleList;
    numIrrPoints = that.numIrrPoints;
    irreducibleStars = that.irreducibleStars;
    equiv = that.equiv;
  }
  return *this;
}

Point Points::getPoint(const long &index) { return Point(*this, index); }

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
  int ikz = idx / (mesh(0) * mesh(1));
  int idx_ = idx - (ikz * mesh(0) * mesh(1));
  int iky = idx_ / mesh(0);
  int ikx = mod(idx_, mesh(0));
  Eigen::Vector3d p;
  p(0) = ikx / (double)mesh(0) + offset(0);
  p(1) = iky / (double)mesh(1) + offset(1);
  p(2) = ikz / (double)mesh(2) + offset(2);
  return p;
}

long Points::getIndex(const Eigen::Vector3d &point) {
  long ik = isPointStored(point);
  if (ik == -1) {
    Error e("getIndex found a point not falling on the wavevector mesh");
  }
  return ik;
}

long Points::isPointStored(const Eigen::Vector3d &crystalCoords) {

  if (numPoints == mesh.prod()) { // full list is faster
    Eigen::Vector3i p;
    // multiply by grid, so that p now contains integers
    double diff = 0.;
    for (int i : {0, 1, 2}) {
      // bring the point to integer coordinates
      double x = (crystalCoords(i) - offset(i)) * mesh(i);
      // check that p is indeed a point commensurate to the mesh.
      diff += round(x) - x;
      // fold in Brillouin zone in range [0,mesh-1]
      p(i) = mod(int(round(x)), mesh(i));
    }
    if (diff >= 1.0e-6) {
      long ik = -1;
      return ik;
    }
    long ik = p(2) * mesh(0) * mesh(1) + p(1) * mesh(0) + p(0);
    return ik;
  } else {
    for (long ikTest = 0; ikTest < numPoints; ikTest++) {
      Eigen::Vector3d kTest = getPointCoords(ikTest, Points::crystalCoords);
      Eigen::Vector3d diff = kTest - crystalCoords;
      for (int i : {0, 1, 2}) {
        diff(i) -= std::floor(diff(i) + 1.0e-8);
      }
      if (diff.squaredNorm() < 1.0e-5) {
        return ikTest;
      }
    }
    long ik = -1;
    return ik;
  }
}

// change of basis methods

Eigen::Vector3d Points::crystalToCartesian(const Eigen::Vector3d &point) {
  return crystal.getReciprocalUnitCell() * point;
}

Eigen::Vector3d Points::cartesianToCrystal(const Eigen::Vector3d &point) {
  Eigen::Vector3d p = crystal.getReciprocalUnitCell().inverse() * point;
  return p;
}

Eigen::Vector3d Points::bzToWs(const Eigen::Vector3d &point, const int &basis) {

  Eigen::Vector3d pointCrystal = point;
  if (basis == cartesianCoords) {
    pointCrystal = cartesianToCrystal(pointCrystal);
  }

  // fold to BZ first
  pointCrystal = foldToBz(pointCrystal, crystalCoords);

  Eigen::Vector3d pointCart = crystalToCartesian(pointCrystal);

  double norm2 = pointCart.squaredNorm();
  long iws = 0;
  Eigen::Vector3d constShift;
  constShift.setConstant(1.0e-8);
  for (long iG = 1; iG < gVectors.cols(); iG++) {
    Eigen::Vector3d thisG = gVectors.col(iG);
    double thisNorm2 = (pointCart + thisG + constShift).squaredNorm();
    if (thisNorm2 < norm2) {
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
  if (basis == cartesianCoords) {
    pointCrystal = cartesianToCrystal(pointCrystal);
  }

  for (int i : {0, 1, 2}) {
    // note: the small double added takes into account for rounding errors
    // which may arise e.g. when the points come from rotations, or machine
    // precision errors in the conversion cartesian to crystal
    // is only a problem in the (unlikely) case of huge k-point meshes (10^8)
    pointCrystal(i) -= std::floor(pointCrystal(i) + 1.0e-8);
  }

  Eigen::Vector3d point2 = pointCrystal;
  if (basis == cartesianCoords) {
    point2 = crystalToCartesian(point2);
  }
  return point2;
}

Crystal &Points::getCrystal() { return crystal; }

long Points::getNumPoints() { return numPoints; }

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

std::tuple<Eigen::Vector3i, Eigen::Vector3d>
Points::findMesh(const Eigen::Matrix<double, 3, Eigen::Dynamic> &testPoints) {
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
        if ((*it - *s.begin() > 1.0e-6) && (1. - *it + *s.begin() > 1.0e-6)) {
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

////////////////////////////////////////////////////////////////

Point::Point(Points &points_, long index_, Eigen::Vector3d umklappVector_)
    : points(points_) {
  umklappVector = umklappVector_;
  index = index_;
}

// copy constructor
Point::Point(const Point &that)
    : umklappVector(that.umklappVector), index(that.index),
      points(that.points) {}

// copy assignment
Point &Point::operator=(const Point &that) {
  if (this != &that) {
    umklappVector = that.umklappVector;
    index = that.index;
    points = that.points;
  }
  return *this;
}

long Point::getIndex() { return index; }

Eigen::Vector3d Point::getCoords(const int &basis, const bool &inWignerSeitz) {
  if ((basis != crystalCoords_) && (basis != cartesianCoords_)) {
    Error e("Point getCoordinates: basis must be crystal or cartesian");
  }

  Eigen::Vector3d coords = points.getPointCoords(index, basis);
  if (inWignerSeitz) {
    coords = points.bzToWs(coords, basis);
  }
  return coords;
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

void Points::setIrreduciblePoints(std::vector<Eigen::MatrixXd> *groupVelocities) {
  // go through the list of wavevectors and find the irreducible points

  const double epsilon = 1.0e-5;

  equiv.resize(numPoints);
  for (long i = 0; i < numPoints; i++) {
    equiv(i) = i;
  }
  // equiv(i) = i : k-point i is not equivalent to any other (irreducible)
  // equiv(i)!=i : k-point i is equivalent to k-point equiv(nk)

  std::vector<SymmetryOperation> symms = crystal.getSymmetryOperations();
  {
    rotationMatricesCrystal.resize(0);
    rotationMatricesCartesian.resize(0);
    Eigen::Matrix3d bg = crystal.getReciprocalUnitCell();
    for (auto symm : symms) {
      Eigen::Matrix3d rotation = symm.rotation;
      rotationMatricesCrystal.push_back(rotation);
      rotation = bg * rotation * bg.inverse();
      rotationMatricesCartesian.push_back(rotation);
    }
  }

  // the identity is assumed to be the first symmetry. Let's check it
  {
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    double x = (rotationMatricesCrystal[0] - I).norm();
    if (x > epsilon) {
      Error e("Development bug: Identity should be the first symmetry");
    }
  }

  for (long ik = 0; ik < numPoints; ik++) {
    // check if this k-point has already been found equivalent to another
    if (equiv(ik) == ik) {
      // check if there are equivalent k-point to this in the list
      // (excepted those previously found to be equivalent to another)
      // check both k and -k
      std::vector<int> possibleRotations;
      for (auto symm : symms) {
        Eigen::Matrix3d rot = symm.rotation;
        Eigen::Vector3d rotatedPoint =
            rot * getPointCoords(ik, Points::crystalCoords);
        long ikRot = isPointStored(rotatedPoint);
        if (ikRot >= 0 && equiv(ikRot) == ikRot) {
          if (ikRot >= ik) {
            equiv(ikRot) = ik;
          } else {
            if (equiv(ikRot) != ik || ikRot < ik) {
              Error e("Error in finding irred kpoints");
            }
          }
        }
      }
    }
  }

  // here we find the rotations mapping the reducible point ot the irreducible
  // Note: I want to see how velocities rotate.
  mapEquivalenceRotationIndex = Eigen::VectorXi::Zero(numPoints);
  mapEquivalenceRotationIndex.setConstant(-1);
  if ( groupVelocities == nullptr ) {
    for (long ikRed = 0; ikRed < numPoints; ikRed++) {
      if (equiv(ikRed) == ikRed) { // identity
        mapEquivalenceRotationIndex(ikRed) = 0; // checked above it's identity
      } else {
        long ikIrr = equiv(ikRed);

//        Eigen::Vector3d kRed = getPointCoords(ikRed, Points::crystalCoords);
        Eigen::Vector3d kIrr = getPointCoords(ikIrr, Points::crystalCoords);

        for (unsigned int is = 0; is < symms.size(); is++) {
          Eigen::Matrix3d rot = rotationMatricesCrystal[is];
          Eigen::Vector3d rotatedPoint = rot * kIrr;
          long ikRot = isPointStored(rotatedPoint);
          if (ikRot == ikRed) {
            mapEquivalenceRotationIndex(ikRed) = is;
            break;
          }
        }
      }
    }
  } else {
    if ( long((*groupVelocities).size()) != numPoints ) {
      Error e("setIrreducible: velocities not aligned with the full grid");
    }

    for (long ikRed = 0; ikRed < numPoints; ikRed++) {
      if (equiv(ikRed) == ikRed) { // identity
        mapEquivalenceRotationIndex(ikRed) = 0; // checked above it's identity
      } else {
        long ikIrr = equiv(ikRed);

//        Eigen::Vector3d kRed = getPointCoords(ikRed, Points::crystalCoords);
        Eigen::Vector3d kIrr = getPointCoords(ikIrr, Points::crystalCoords);

        Eigen::MatrixXd irrVels = (*groupVelocities)[ikIrr];
        Eigen::MatrixXd redVels = (*groupVelocities)[ikRed];

        if (irrVels.rows() != redVels.rows()) {
          Error e("Different number of bands at two equivalent points");
        }

        std::vector<int> isSelects;
        std::vector<double> diffs;

        for (unsigned int is = 0; is < symms.size(); is++) {
          Eigen::Vector3d rotatedPoint = rotationMatricesCrystal[is] * kIrr;
          long ikRot = isPointStored(rotatedPoint);
          if (ikRot != ikRed) {
            continue;
          }

          int numBands = irrVels.rows();
          double diff = 0.;
          for (int ib=0; ib<numBands; ib++) {
            Eigen::Vector3d thisIrrVel = irrVels.row(ib);
            Eigen::Vector3d thisRedVel = redVels.row(ib);
            Eigen::Vector3d rotVel = rotationMatricesCartesian[is] * thisIrrVel;
            diff += (rotVel - thisRedVel).squaredNorm();
          }
          isSelects.push_back(is);
          diffs.push_back(diff);
        }

        int minElementIndex = std::min_element(diffs.begin(),diffs.end())
                              - diffs.begin();

        mapEquivalenceRotationIndex(ikRed) = isSelects[minElementIndex];

      }
    }
  }

  // count number of irreducible points
  numIrrPoints = 0;
  for (long ik = 0; ik < numPoints; ik++) {
    if (equiv(ik) == ik) {
      numIrrPoints += 1;
    }
  }

  // this allows us to map (ikIrr in the irrPoints) -> (ikRed in fullpoints)
  mapIrreducibleToReducibleList = Eigen::VectorXi::Zero(numIrrPoints);
  long i = 0;
  for (long j = 0; j < numPoints; j++) {
    if (equiv(j) == j) { // if is irreducible
      mapIrreducibleToReducibleList(i) = j;
      i += 1;
    }
  }

  // this allows us to map (ikRed in fullpoints) -> (ikIrr in the irrPoints)
  // basically the inverse of mapIrrToRedList
  mapReducibleToIrreducibleList = Eigen::VectorXi(numPoints);
  for (long ik = 0; ik < numPoints; ik++) {
    int ikIrr = equiv(ik); // map to the irreducible in fullPoints
    long ikIrr2 = -1;
    for (long i = 0; i < numIrrPoints; i++) {
      if (mapIrreducibleToReducibleList(i) == ikIrr) {
        ikIrr2 = i;
        break;
      }
    }
    if (ikIrr2 == -1) {
      Error e("Failed building irreducible points mapRedToIrrList");
    }
    mapReducibleToIrreducibleList(ik) = ikIrr2;
  }

  for (long ikIrr = 0; ikIrr < numIrrPoints; ikIrr++) {
    std::vector<long> thisStar;
    long ikIrrFull = mapIrreducibleToReducibleList(ikIrr);
    for (long ik = 0; ik < numPoints; ik++) {
      if (equiv(ik) == ikIrrFull) {
        thisStar.push_back(ik);
      }
    }
    irreducibleStars.push_back(thisStar);
  }
}

std::vector<long> Points::irrPointsIterator() {
  std::vector<long> iter;
  if (numIrrPoints == 0) {
    for (long ik = 0; ik < numPoints; ik++) {
      iter.push_back(ik);
    }
  } else {
    for (long ik = 0; ik < numPoints; ik++) {
      if (equiv(ik) == ik) {
        iter.push_back(ik);
      }
    }
  }
  return iter;
}

std::vector<long> Points::parallelIrrPointsIterator() {
  auto v = irrPointsIterator();
  //
  auto divs = mpi->divideWork(v.size());
  long start = divs[0];
  long stop = divs[1];
  //
  std::vector<long> iter(v.begin() + start, v.begin() + stop);
  return iter;
}

long Points::asIrreducibleIndex(const long &ik) {
  if (numIrrPoints > 0) { // we set the irreducible points
    assert(ik < numPoints && ik > 0);
    return mapReducibleToIrreducibleList(ik);
  } else { // no symmetries set, only the identity symmetry exists
    return ik;
  }
}

long Points::asReducibleIndex(const long &ik) {
  if (numIrrPoints > 0) { // we set the irreducible points
    assert(ik < numIrrPoints && ik > 0);
    return mapIrreducibleToReducibleList(ik);
  } else { // no symmetries set, only the identity symmetry exists
    return ik;
  }
}

std::tuple<long, Eigen::Matrix3d>
Points::getRotationToIrreducible(const Eigen::Vector3d &x, const int &basis) {
  Eigen::Vector3d xCryst;
  if (basis == cartesianCoords) {
    xCryst = cartesianToCrystal(x);
  } else {
    xCryst = x;
  }
  // find the index in the list
  long ik = getIndex(xCryst);

  if (numIrrPoints > 0) { // if irreducible points have been set
    // find rotation such that rotation * qFull = qRed
    Eigen::Matrix3d rot;
    if (basis == crystalCoords) {
      rot = rotationMatricesCrystal[mapEquivalenceRotationIndex(ik)].inverse();
    } else {
      rot = rotationMatricesCartesian[mapEquivalenceRotationIndex(ik)].inverse();
    }
    // also, we add the index of the irreducible point to which x is mapped
    return {equiv(ik), rot};
  } else {
    Eigen::Matrix3d identity;
    identity.setIdentity();
    return {ik, identity};
  }
}

std::vector<Eigen::Matrix3d> Points::getRotationsStar(const long &ik) {
  if (numIrrPoints > 0) {
    // first, map ik to [0,N_irr]
    long ikIrr = mapReducibleToIrreducibleList(ik);
    std::vector<long> starIndices = irreducibleStars[ikIrr];
    std::vector<Eigen::Matrix3d> rotations;
    for (long ikFull : starIndices) {
      Eigen::Matrix3d rot;
      rot = rotationMatricesCartesian[mapEquivalenceRotationIndex(ikFull)];
      rotations.push_back(rot);
    }
    // list of rotations such that qStar = rotation * qIrr
    return rotations;
  } else {
    std::vector<Eigen::Matrix3d> rotations;
    Eigen::Matrix3d identity;
    identity.setIdentity();
    rotations.push_back(identity);
    return rotations;
  }
}

std::vector<long> Points::getReduciblesFromIrreducible(const long &ik) {
  long ikIrr = mapReducibleToIrreducibleList(ik);
  return irreducibleStars[ikIrr];
}
