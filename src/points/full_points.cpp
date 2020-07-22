#include "full_points.h"
#include "eigen.h"
#include "utilities.h"

FullPoints::FullPoints(Crystal &crystal_, const Eigen::Vector3i &mesh_,
        const Eigen::Vector3d &offset_) :
        Points(crystal_, mesh_, offset_) {
}

// copy constructor
FullPoints::FullPoints(const FullPoints &that) :
        Points(that) {
}

// assignment operator
FullPoints& FullPoints::operator=(const FullPoints &that) {
    if (this != &that) {
        crystal = that.crystal;
        mesh = that.mesh;
        offset = that.offset;
        numPoints = that.numPoints;
        gVectors = that.gVectors;
        igVectors = that.igVectors;
    }
    return *this;
}

Point FullPoints::getPoint(const long &index) {
    return Point(*this, index);
}
