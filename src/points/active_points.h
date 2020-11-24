#ifndef ACTIVEPOINTS_H
#define ACTIVEPOINTS_H

#include "points.h"

/** Class for storing an "active" list of wavevectors, i.e. a selection of
 * points taken from a monkhorst-pack grid of points.
 */
class ActivePoints : public Points {
 protected:
  Points &parentPoints;
  Eigen::MatrixXd pointsList;

  Eigen::VectorXi filteredToFullIndeces;
  long fullToFilteredIndeces(const long &indexIn);
 public:
  /** Default constructor
   * @param parentPoints: the "parent" Points object, from which we take a
   * selection of points
   * @param filter: a vector of integers of "filtered" points, i.e. the
   * indices of the wavevectors in parentPoints that we want to keep.
   */
  ActivePoints(Points &parentPoints_, Eigen::VectorXi filter_);

  /** Copy constructor
   */
  ActivePoints(const ActivePoints &obj);

  /** Copy assignment operator
   */
  ActivePoints &operator=(const ActivePoints &obj);

  /** Get the wavevector index given the crystal coordinates of a wavevector.
   * @param point: the wavevector in crystal coordinates.
   * @return index: the index of the wavevector in the range [0,numPoints[
   */
  long getIndex(const Eigen::Vector3d &coords);

  /** Returns a Point object given its integer index.
   * The Point object is used to move around the code the coordinates of the
   * wavevector.
   * @param index: the integer wavevector index ranging in [0,numPoints[
   * @return point: a point object.
   */
  Point getPoint(const long &index);

  /** Returns the Points object from which the ActivePoints have been built.
   */
  Points getParentPoints();

  /** Get the coordinates of a wavevector from its index.
   * @param index: the index of the desired wavevector.
   * @param basis: specify the basis to be used for the output coordinates.
   * Either Points::crystalCoords or Points::cartesianCoords.
   * @return wavevector: the coordinates of the desired wavevector.
   */
  Eigen::Vector3d getPointCoords(const long &index, const int &basis =
  crystalCoords);
};

#endif
