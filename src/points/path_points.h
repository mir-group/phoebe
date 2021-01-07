#ifndef PATH_POINTS_H
#define PATH_POINTS_H

#include "points.h"
#include "full_points.h"

/** Class for storing wavevectors along a path of the Brillouin zone.
 */
class PathPoints: public FullPoints {
public:
    /** Default constructor.
     * @param crystal: the crystal object that defines the Brillouin zone.
     * @param pathExtrema: extrema of segments defining a path of points.
     * PathExtrema has size (numSegments, 2, 3), where the first index spans
     * the segments of the path, 2 is the beginning and the last element of the
     * segment, and 3 is the cartesian coordinates of the extrema.
     * Point coordinates must be in crystal coordinates.
     * Example: if we want a path G->X->L, we have: 2 segments: G->X and X->L
     * pathExtrema(0,0,:) = G. coordinates
     * pathExtrema(0,1,:) = X. coordinates
     * pathExtrema(1,0,:) = X. coordinates
     * pathExtrema(1,1,:) = L. coordinates
     * @param delta: segments are populated with points, with a distance
     * "delta" between different points in crystal coordinates.
     */
    PathPoints(Crystal &crystal_, const Eigen::Tensor<double, 3> &pathExtrema,
            const double &delta);

    /** Copy constructor
     */
    PathPoints(const PathPoints &obj);

    /** Copy assignment operator
     */
    PathPoints& operator=(const PathPoints &obj);

    /** Returns a Point object given its integer index.
     * The Point object is used to move around the code the coordinates of the
     * wavevector.
     * @param index: the integer wavevector index ranging in [0,numPoints[
     * @return point: a point object.
     */
    Point getPoint(const int &index);

    /** Get the wavevector index given the crystal coordinates of a wavevector.
     * @param point: the wavevector in crystal coordinates.
     * @return index: the index of the wavevector in the range [0,numPoints[
     */
    int getIndex(const Eigen::Vector3d &coordinates);

    /** Get the coordinates of a wavevector from its index.
     * @param index: the index of the desired wavevector.
     * @param basis: specify the basis to be used for the output coordinates.
     * Either Points::crystalCoords or Points::cartesianCoords.
     * @return wavevector: the coordinates of the desired wavevector.
     */
    Eigen::Vector3d getPointCoords(const int &index, const int &basis =
            Points::crystalCoords);
protected:
    Eigen::Matrix<double, 3, Eigen::Dynamic> pointsList;
};

#endif
