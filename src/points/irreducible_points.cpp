#include "irreducible_points.h"
#include "eigen.h"
#include "utilities.h"
#include <iomanip>

IrreduciblePoints::IrreduciblePoints(Points &parentPoints_) :
    Points(parentPoints_.getCrystal(), std::get<0>(parentPoints_.getMesh()),
           std::get<1>(parentPoints_.getMesh())), parentPoints(parentPoints_) {

  // generate a list of points compatible with the symmetries of the system

  const double epsilon = 1.0e-5;

  int numParentPoints = parentPoints.getNumPoints();

  Eigen::VectorXd tmpWeight(numParentPoints);
  Eigen::VectorXi equiv(numParentPoints);
  for (long i = 0; i < numParentPoints; i++) {
    equiv(i) = i;
  }
  // equiv(nk) =nk : k-point nk is not equivalent to any previous k-point
  // equiv(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)

  std::vector<SymmetryOperation> symms = crystal.getSymmetryOperations();
  {
    Eigen::Matrix3d bg = crystal.getReciprocalUnitCell();
    for (auto symm : symms) {
      Eigen::Matrix3d rotation = symm.rotation;
      rotationMatricesCrystal.push_back(rotation);
      rotation = bg * rotation * bg.inverse();
      rotationMatricesCartesian.push_back(rotation);
    }
  }
  mapEquivalenceRotationIndex = Eigen::VectorXi::Zero(numParentPoints);
  mapEquivalenceRotationIndex.setConstant(0);

  // the identity is assumed to be the first symmetry. Let's check it
  {
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    double x = (rotationMatricesCrystal[0] - I).norm();
    if (x > epsilon) {
      Error e("Development bug: Identity should be the first symmetry");
    }
  }

  for (long ik = 0; ik < numParentPoints; ik++) {
    // check if this k-point has already been found equivalent to another
    if (equiv(ik) == ik) {
      tmpWeight(ik) = 1.;
      // check if there are equivalent k-point to this in the list
      // (excepted those previously found to be equivalent to another)
      // check both k and -k
      int is = -1;
      for (auto symm : symms) {
        is += 1;

        Eigen::Matrix3d s = symm.rotation;

        Eigen::Vector3d rotatedPoint =
            s * parentPoints.getPointCoords(ik, Points::crystalCoords);

        double xx = rotatedPoint(0) * mesh(0) - offset(0);
        double yy = rotatedPoint(1) * mesh(1) - offset(1);
        double zz = rotatedPoint(2) * mesh(2) - offset(2);
        bool inTheList = (abs(xx - round(xx)) <= epsilon
            && abs(yy - round(yy)) <= epsilon
            && abs(zz - round(zz)) <= epsilon);
        if (inTheList) { // i.e. can be rotated to another point in the grid
          int n = parentPoints.getIndex(rotatedPoint);
          if (n > ik && equiv(n) == n) {
            equiv(n) = ik;
            tmpWeight(ik) += 1.;
            mapEquivalenceRotationIndex(n) = is;
          } else {
            if (equiv(n) != ik || n < ik) {
              Error e("Error in finding irred kpoints");
            }
          }
        }
      }
    }
  }

  // count irreducible points and order them

  numPoints = 0;
  for (long ik = 0; ik < numParentPoints; ik++) {
    if (equiv(ik) == ik) {
      numPoints += 1;
    }
  }

  pointsList = Eigen::MatrixXd::Zero(3, numPoints);
  weights = Eigen::VectorXd::Zero(numPoints);

  // this allows us to map (ikIrr in the irrPoints) -> (ikRed in fullpoints)
  mapIrreducibleToReducible = Eigen::VectorXi(numPoints);
  long i = 0;
  for (long j = 0; j < numParentPoints; j++) {
    if (equiv(j) == j) { // if is irreducible
      weights(i) = tmpWeight(j);
      pointsList.col(i) = parentPoints.getPointCoords(j, Points::crystalCoords);
      mapIrreducibleToReducible(i) = j;
      i += 1;
    }
  }

  // this allows us to map (ikRed in fullpoints) -> (ikIrr in the irrPoints)
  mapReducibleToIrreducible = Eigen::VectorXi(numParentPoints);
  for (int ik=0; ik<numParentPoints; ik++) {
    int ikIrr = equiv(ik); // map to the irreducible in fullPoints
    Eigen::Vector3d kFull = parentPoints.getPointCoords(ikIrr, Points::crystalCoords);
    i = getIndex(kFull);
    mapReducibleToIrreducible(ik) = i;
  }

  // normalize weights to one
  weights /= weights.sum();

  for (int ikIrr = 0; ikIrr < numPoints; ikIrr++) {
    std::vector<int> thisStar;
    int ikIrrFull = mapIrreducibleToReducible(ikIrr);
    for (int ik = 0; ik < numParentPoints; ik++) {
      if (equiv(ik) == ikIrrFull) {
        thisStar.push_back(ik);
      }
    }
    irreducibleStars.push_back(thisStar);
  }

  // Test 1:
  // rotation from qIrr to qRed in crystal coordinates
  for ( int ikRed=0; ikRed<numParentPoints; ikRed++ ) {
    auto qRed = parentPoints.getPointCoords(ikRed, crystalCoords);
    qRed = foldToBz(qRed, crystalCoords);
    auto qIrr = getPointCoords(mapReducibleToIrreducible(ikRed), crystalCoords);

    // Equivalently:
    auto qIrr2 = parentPoints.getPointCoords(equiv(ikRed), crystalCoords);
    if ( (qIrr2 - qIrr).norm() > 1.0e-8 ) {
      Error e("Unexpected behavior");
    }

    qIrr = foldToBz(qIrr, crystalCoords);

    auto rotation = rotationMatricesCrystal[mapEquivalenceRotationIndex(ikRed)];
    Eigen::Vector3d qRot = rotation * qIrr;
    qRot = foldToBz(qRot, crystalCoords);

    if ( (qRed - qRot).norm() > 1.0e-8 ) {
      Error e("Unexpected behavior");
    }
  }

  // Test 2
  // same as above, in cartesian coordinates
  for ( int ikRed=0; ikRed<numParentPoints; ikRed++ ) {
    auto qRed = parentPoints.getPointCoords(ikRed, cartesianCoords);
//    qRed = foldToBz(qRed, cartesianCoords);
    auto qIrr = getPointCoords(mapReducibleToIrreducible(ikRed), cartesianCoords);

    // Equivalently:
    auto qIrr2 = parentPoints.getPointCoords(equiv(ikRed), cartesianCoords);
    if ( (qIrr2 - qIrr).norm() > 1.0e-8 ) {
      Error e("Unexpected behavior 1");
    }

    qIrr = foldToBz(qIrr, cartesianCoords);

    auto rotation = rotationMatricesCartesian[mapEquivalenceRotationIndex(ikRed)];
    Eigen::Vector3d qRot = rotation * qIrr;
//    qRot = foldToBz(qRot, cartesianCoords);

    qRed = bzToWs(qRed, cartesianCoords);
    qRot = bzToWs(qRot, cartesianCoords);

    if ( (qRed - qRot).norm() > 1.0e-8 ) {
      Error e("Unexpected behavior 2");
    }
  }
}

// copy constructor
IrreduciblePoints::IrreduciblePoints(const IrreduciblePoints &that) :
    Points(that),
    parentPoints(that.parentPoints),
    mapReducibleToIrreducible(that.mapReducibleToIrreducible),
    mapIrreducibleToReducible(that.mapIrreducibleToReducible),
    irreducibleStars(that.irreducibleStars),
    mapEquivalenceRotationIndex(that.mapEquivalenceRotationIndex),
    rotationMatricesCrystal(that.rotationMatricesCrystal),
    rotationMatricesCartesian(that.rotationMatricesCartesian),
    pointsList(that.pointsList),
    weights(that.weights) {
}

// assignment operator
IrreduciblePoints &IrreduciblePoints::operator=(const IrreduciblePoints &that) {
  Points::operator=(that);
  if (this != &that) {
    parentPoints = that.parentPoints;
    mapReducibleToIrreducible = that.mapReducibleToIrreducible;
    mapIrreducibleToReducible = that.mapIrreducibleToReducible;
    irreducibleStars = that.irreducibleStars;
    mapEquivalenceRotationIndex = that.mapEquivalenceRotationIndex;
    rotationMatricesCartesian = that.rotationMatricesCartesian;
    rotationMatricesCrystal = that.rotationMatricesCrystal;
    pointsList = that.pointsList;
    weights = that.weights;
  }
  return *this;
}

Point IrreduciblePoints::getPoint(const long &index) {
  return Point(*this, index);
}

Eigen::Vector3d IrreduciblePoints::getPointCoords(const long &index,
                                                  const int &basis) {
  if (basis != Points::crystalCoords && basis != Points::cartesianCoords) {
    Error e("Wrong basis for getPoint");
  }
  Eigen::Vector3d pointCrystal = pointsList.col(index);
  if (basis == Points::crystalCoords) {
    return pointCrystal;
  } else {
    Eigen::Vector3d pointCartesian = crystalToCartesian(pointCrystal);
    return pointCartesian;
  }
}

// I don't assume to be working with a uniform mesh
long IrreduciblePoints::getIndex(const Eigen::Vector3d &coords) {
  // in this case there is no order, so we just search through a loop
  bool found = false;
  long counter;
  for (counter = 0; counter < numPoints; counter++) {
    if ((pointsList.col(counter) - coords).norm() < 1.0e-8) {
      found = true;
      break;
    }
  }
  if (!found) {
    Error e("Point not found in Irreducible Points");
  }
  return counter;
}

double IrreduciblePoints::getWeight(const long &ik) {
  return weights(ik);
}

std::tuple<int, Eigen::Matrix3d> IrreduciblePoints::getRotationToIrreducible(
    const Eigen::Vector3d &x, const int & basis) {

  // x in input is the coordinates of a wavevector in the full grid

  // first we bring it in cartesian coordinates
  Eigen::Vector3d pointCrystal = x;
  if (basis == cartesianCoords) {
    pointCrystal = cartesianToCrystal(pointCrystal);
  }

  // find index of x in the full grid
  int ikFull = parentPoints.getIndex(pointCrystal);

  // find index of irreducible point in the irreducible grid
  int ikIrr = mapReducibleToIrreducible(ikFull);

  // find rotation such that rotation * qFull = qRed
  int iRot = mapEquivalenceRotationIndex(ikFull);
  Eigen::Matrix3d rotation;
  if (basis == crystalCoords) {
    rotation = rotationMatricesCrystal[iRot].inverse();
  } else {
    rotation = rotationMatricesCartesian[iRot].inverse();
  }
  return {ikIrr, rotation};
}

std::vector<Eigen::Matrix3d> IrreduciblePoints::getRotationsStar(const int &ik) {
  // note: ik is the index of the irreducible point in the irreducible list

  std::vector<int> fullIndeces = irreducibleStars[ik];

  std::vector<Eigen::Matrix3d> rotations;
  for (int ikFull : fullIndeces) {
    Eigen::Matrix3d rot;
    rot = rotationMatricesCartesian[mapEquivalenceRotationIndex(ikFull)];
    rotations.push_back(rot);
  }
  // list of rotations such that qStar = rotation * qIrr
  return rotations;
}

Points IrreduciblePoints::getParentPoints() {
  return parentPoints;
}

int IrreduciblePoints::irreducibleToReducible(const int &ikIrr) {
  return mapIrreducibleToReducible(ikIrr);
}
