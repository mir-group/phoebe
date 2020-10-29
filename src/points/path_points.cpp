#include "path_points.h"
#include "constants.h"
#include "eigen.h"

PathPoints::PathPoints(Crystal &crystal_,
                       const Eigen::Tensor<double, 3> &pathExtrema,
                       const double &delta)
    : FullPoints(crystal_, Eigen::Vector3i::Constant(1),
                 Eigen::Vector3d::Zero()) {

  // build the list of points
  std::vector<Eigen::Vector3d> points;
  Eigen::Vector3d p0, p1; //, thisDelta;

  // initialize the path
  p0(0) = pathExtrema(0, 0, 0);
  p0(1) = pathExtrema(0, 0, 1);
  p0(2) = pathExtrema(0, 0, 2);
  points.push_back(p0);

  // we loop over the segments provided in user input
  for (long i = 0; i < pathExtrema.dimension(0); i++) {
    // load coords of the extrema of the segment
    p0(0) = pathExtrema(i, 0, 0);
    p0(1) = pathExtrema(i, 0, 1);
    p0(2) = pathExtrema(i, 0, 2);
    p1(0) = pathExtrema(i, 1, 0);
    p1(1) = pathExtrema(i, 1, 1);
    p1(2) = pathExtrema(i, 1, 2);

    // delta may not divide the interval exactly
    // so, we find the closest one
    long nk = abs(long((p1 - p0).norm() / delta));

    // now we build the points of the segment
    std::vector<Eigen::Vector3d> segmentPoints;
    for (long j = 0; j <= nk; j++) {
      Eigen::Vector3d thisP;
      thisP(0) = (p1(0) - p0(0)) / nk * j + p0(0);
      thisP(1) = (p1(1) - p0(1)) / nk * j + p0(1);
      thisP(2) = (p1(2) - p0(2)) / nk * j + p0(2);
      segmentPoints.push_back(thisP);
    }

    // now we update the original list of points
    // but we have to check to not add the same extrema twice
    // work with the first element
    Eigen::Vector3d lastP = points[points.size() - 1];
    Eigen::Vector3d thisP = segmentPoints[0];
    bool addIt = (thisP - lastP).norm() > epsilon8;
    if (addIt) {
      points.push_back(thisP);
    }

    // select first element of the segment
    auto p = std::begin(segmentPoints);
    ++p;
    // and then add all the rest of the segment
    for (auto end = std::end(segmentPoints); p != end; ++p) {
      // iterate over the rest of the container
      points.push_back(*p);
    }
  }

  numPoints = points.size();
  pointsList = Eigen::MatrixXd::Zero(3, numPoints);
  long i = 0;
  for (auto p : points) {
    pointsList.col(i) = p;
    i += 1;
  }
}

// copy constructor
PathPoints::PathPoints(const PathPoints &that)
    : FullPoints(that), pointsList(that.pointsList) {}

// assignment operator
PathPoints &PathPoints::operator=(const PathPoints &that) {
  if (this != &that) {
    crystal = that.crystal;
    mesh = that.mesh;
    offset = that.offset;
    numPoints = that.numPoints;
    gVectors = that.gVectors;
    pointsList = that.pointsList;
  }
  return *this;
}

Point PathPoints::getPoint(const long &index) { return Point(*this, index); }

Eigen::Vector3d PathPoints::getPointCoords(const long &index,
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

long PathPoints::getIndex(const Eigen::Vector3d &coords) {
  // in this case there is no order, so we just search through a loop
  long counter = 0;
  for (counter = 0; counter < numPoints; counter++) {
    if ((pointsList.col(counter) - coords).norm() < 1.0e-8) {
      break;
    }
  }
  return counter;
}
