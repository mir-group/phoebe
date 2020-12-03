#include "full_points.h"
#include "eigen.h"
#include "utilities.h"

FullPoints::FullPoints(Crystal &crystal_, const Eigen::Vector3i &mesh_,
                       const Eigen::Vector3d &offset_)
    : Points(crystal_, mesh_, offset_) {}

// copy constructor
FullPoints::FullPoints(const FullPoints &that) : Points(that) {}

// assignment operator
FullPoints &FullPoints::operator=(const FullPoints &that) {
  Points::operator=(that);
  return *this;
}

Point FullPoints::getPoint(const long &index) { return Point(*this, index); }

long FullPoints::isPointStored(const Eigen::Vector3d &crystalCoords) {
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
}
