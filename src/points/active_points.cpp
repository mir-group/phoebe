#include "active_points.h"
#include "eigen.h"

ActivePoints::ActivePoints(Points &parentPoints_, Eigen::VectorXi filter) :
        Points(parentPoints_.getCrystal(),
                std::get < 0 > (parentPoints_.getMesh()),
                std::get < 1 > (parentPoints_.getMesh())), parentPoints(
                parentPoints_) {
  // this contain the list of indices of the points in the FullPoints class
  // which we want to include in the ActivePoints class
  filteredToFullIndices = filter;

  numPoints = filteredToFullIndices.size();

  // we then construct the list of points
  Eigen::MatrixXd pointsList_(3, numPoints);

  Eigen::Vector3d x;

  long ik;
  for (long ikNew = 0; ikNew < numPoints; ikNew++) {
    ik = filteredToFullIndices(ikNew);
    x = parentPoints.getPointCoords(ik);
    pointsList_.col(ikNew) = x;
  }
  pointsList = pointsList_;
}

// copy constructor
ActivePoints::ActivePoints(const ActivePoints &that) :
        Points(that), parentPoints(that.parentPoints), pointsList(
                that.pointsList), filteredToFullIndices(
                that.filteredToFullIndices) {
}

// copy assignment operator
ActivePoints& ActivePoints::operator=(const ActivePoints &that) {
  if (this != &that) {
    crystal = that.crystal;
    mesh = that.mesh;
    offset = that.offset;
    numPoints = that.numPoints;
    gVectors = that.gVectors;
    parentPoints = that.parentPoints;
    filteredToFullIndices = that.filteredToFullIndices;
    pointsList = that.pointsList;
  }
  return *this;
}

Point ActivePoints::getPoint(const long &index) {
  return Point(*this, index);
}

Points ActivePoints::getParentPoints() {
  return parentPoints;
}

Eigen::Vector3d ActivePoints::getPointCoords(const long &index,
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

long ActivePoints::getIndex(const Eigen::Vector3d &coords) {
  // we take advantage of the fact that the parent points have an order
  long indexFull = parentPoints.getIndex(coords);
  long ik = fullToFilteredIndices(indexFull);
  return ik;
}

long ActivePoints::fullToFilteredIndices(const long &indexIn) {
  // Note: this function could obviously be made much faster if you could
  // save in memory the map of every point in the Full list into the
  // ActivePoints list (and 0 if there's no mapping.
  // But the list of kpoints might be too large!
  long target = -1;
  for (long ik = 0; ik < numPoints; ik++) {
    if (indexIn == filteredToFullIndices(ik)) {
      target = ik;
      break;
    }
  }
  if (target == -1) {
    Error e("Couldn't find the desired kpoint");
  }
  return target;
}
