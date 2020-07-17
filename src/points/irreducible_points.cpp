#include "irreducible_points.h"
#include "eigen.h"
#include "utilities.h"

IrreduciblePoints::IrreduciblePoints(Points &parentPoints_) :
    parentPoints(parentPoints_),  Points(parentPoints_.getCrystal(), std::get<0>(parentPoints_.getMesh()),
        std::get<1>(parentPoints_.getMesh())) {

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

  for ( auto symm : symms ) {
    Eigen::Matrix3d rotation = symm.rotation;
    Eigen::Matrix3d bg = crystal.getReciprocalUnitCell();
    rotation = bg * rotation * bg.inverse();
    rotationMatrices.push_back(rotation);
  }

  mapEquivalenceRotationIndex = Eigen::VectorXi::Zero(numParentPoints);
  mapEquivalenceRotationIndex.setConstant(1);

  // the identity is assumed to be the first symmetry. Let's check it
  {
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    double x = (rotationMatrices[0] - I).norm();
    if ( x > epsilon ) {
      Error e("Development bug: rotation matrix should be the first symmetry");
    }
  }

  for (long ik = 0; ik < numParentPoints; ik++) {
    // check if this k-point has already been found equivalent to another
    if (equiv(ik) == ik) {
      tmpWeight(ik) = 1.;
      // check if there are equivalent k-point to this in the list
      // (excepted those previously found to be equivalent to another)
      // check both k and -k
      int is=-1;
      for (auto symm : symms) {
        is += 1;

        Eigen::Matrix3d s = symm.rotation;

        Eigen::Vector3d rotatedPoint = s * parentPoints.getPointCoords(ik, Points::crystalCoords);

        for (long i = 0; i < 3; i++) {
          rotatedPoint(i) -= round(rotatedPoint(i));
        }
        double xx = rotatedPoint(0) * mesh(0) - offset(0);
        double yy = rotatedPoint(1) * mesh(1) - offset(1);
        double zz = rotatedPoint(2) * mesh(2) - offset(2);
        bool inTheList = (abs(xx - round(xx)) <= epsilon
                       && abs(yy - round(yy)) <= epsilon
                       && abs(zz - round(zz)) <= epsilon);
        if (inTheList) {
          int ix = mod(long(round(rotatedPoint(0) * mesh(0)
                                  - offset(0) + 2 * mesh(0))), mesh(0));
          int iy = mod(long(round(rotatedPoint(1) * mesh(1)
                                  - offset(1) + 2 * mesh(1))), mesh(1));
          int iz = mod(long(round(rotatedPoint(2) * mesh(2)
                                  - offset(2) + 2 * mesh(2))), mesh(2));
          int n = iz + iy * mesh(2) + ix * mesh(1) * mesh(2);
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

  long numPoints = 0;
  for ( long ik = 0; ik < numParentPoints; ik++) {
    if (equiv(ik) == ik) {
      numPoints += 1;
    }
  }

  pointsList = Eigen::MatrixXd::Zero(3, numPoints);
  weights = Eigen::VectorXd::Zero(numPoints);

  mapIrreducibleToReducible = Eigen::VectorXi(numPoints);
  mapReducibleToIrreducible = Eigen::VectorXi(numParentPoints);

  long i = 0;
  for (long j = 0; j < numParentPoints; j++) {
    if (equiv(j) == j) {
      weights(i) = tmpWeight(j);
      pointsList.col(i) = parentPoints.getPointCoords(j,Points::crystalCoords);

      mapIrreducibleToReducible(i) = j;
      mapReducibleToIrreducible(j) = i;

      i += 1;
    }
  }

  // normalize weights to one
  weights /= weights.sum();

  for ( int ikIrr=0; ikIrr<numPoints; ikIrr++ ) {
    std::vector<int> thisStar;
    int ikIrrFull = mapIrreducibleToReducible(ikIrr);
    for ( int ik=0; ik<numParentPoints; ik++ ) {
      if ( equiv(ik) == ikIrrFull) {
        thisStar.push_back(ik);
      }
    }
    irreducibleStars.push_back(thisStar);
  }

}

// copy constructor
IrreduciblePoints::IrreduciblePoints(const IrreduciblePoints &that) :
        Points(that), parentPoints(that.parentPoints),
        mapReducibleToIrreducible(that.mapReducibleToIrreducible),
        mapIrreducibleToReducible(that.mapIrreducibleToReducible) {
}

// assignment operator
IrreduciblePoints& IrreduciblePoints::operator=(const IrreduciblePoints &that) {
    if (this != &that) {
        crystal = that.crystal;
        mesh = that.mesh;
        offset = that.offset;
        numPoints = that.numPoints;
        gVectors = that.gVectors;
        igVectors = that.igVectors;
        mapReducibleToIrreducible = that.mapReducibleToIrreducible;
        mapIrreducibleToReducible = that.mapIrreducibleToReducible;
    }
    return *this;
}

Point IrreduciblePoints::getPoint(const long &index) {
    return Point(*this, index);
}

Eigen::Vector3d IrreduciblePoints::getPointCoords(const long &index,
        const int &basis) {
    if (basis != crystalCoords && basis != cartesianCoords) {
        Error e("Wrong basis for getPoint");
    }
    Eigen::Vector3d pointCrystal = pointsList.col(index);
    if (basis == crystalCoords) {
        return pointCrystal;
    } else {
        Eigen::Vector3d pointCartesian = crystalToCartesian(pointCrystal);
        return pointCartesian;
    }
}

long IrreduciblePoints::getIndex(const Eigen::Vector3d &point) {
    // untested
    long ik = getIndex(point);
    long ikIrr = mapReducibleToIrreducible(ik);
    return ikIrr;
}

double IrreduciblePoints::getWeight(const long &ik) {
    return weights(ik);
}

Eigen::VectorXi IrreduciblePoints::getIndexReducibleFromIrreducible(
        const long &indexIrr) {
    std::vector<int> indexVec;
    long sizeStar = 0;
    for (long ik = 0; ik < numPoints; ik++) {
        if (mapReducibleToIrreducible(ik) == indexIrr) {
            indexVec.push_back(ik);
            sizeStar += 1;
        }
    }
    Eigen::VectorXi star(sizeStar);
    for (long ik = 0; ik < sizeStar; ik++) {
        star(ik) = indexVec[ik];
    }
    return star;
}

long IrreduciblePoints::getIndexIrreducibleFromReducible(const long &indexRed) {
    long ik = mapReducibleToIrreducible(indexRed);
    return ik;
}
