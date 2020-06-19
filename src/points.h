#ifndef POINTS_H
#define POINTS_H

#include "crystal.h"
#include "eigen.h"
#include "exceptions.h"

const int crystalCoords_ = 0;
const int cartesianCoords_ = 1;

//// Forward declarator of Points, because Point ants to store Points as member
//// and Points has methods returning Point. Here we avoid recursive dependency.
class Points;

/** Class used to pass a single wavevector
 */
class Point {
public:
  /** Constructor
   * @param index: integer index of the wavevector in the Points object.
   * @param umklappVector: the crystal coordinates of a possible Umklapp
   * vector, used for example in sum or differences between wavevectors.
   * @param points: the points object that this Point object belongs to.
   */
  Point(Points &points_, long index_,
        Eigen::Vector3d umklappVector = Eigen::Vector3d::Zero());

  /** copy constructor
   */
  Point(const Point &that);

  /** copy assignment operator
   */
  Point &operator=(const Point &that);

  /** Get the coordinates of the k-point
   * @param basis: either "cartesian" or "crystal"
   * @param inWignerSeitz: default false, if true, folds point in WS cell.
   * @return coords: a 3d vector of coordinates
   */
  Eigen::Vector3d getCoords(const int &basis = crystalCoords_,
                            const bool &inWignerSeitz = false);

  /** Get the weight of the k-point (used for integrations over the BZ.
   * @return weight: a double.
   */
  double getWeight();

  /** Sum of two wavevectors (this + b)
   * The vector is folded in the Wigner Seitz zone with an Umklapp vector.
   */
  Point operator+(Point &b);

  /** Difference of two wavevectors (this-b).
   * The vector is folded in the Wigner Seitz zone with an Umklapp vector.
   */
  Point operator-(Point &b);

  /** checks whether the Point has been constructed with an Umklapp vector
   */
  bool hasUmklapp();

  /** Get the index of the wavevector in the referenced Points object.
   */
  long getIndex();

private:
  Eigen::Vector3d umklappVector;
  long index;
  Points &points;
};

/** Base class for storing wavevectors in the Brillouin zone.
 * The base class mostly refers to a uniform grid sampling the Brillouin zone.
 * Different subclasses will specialize to different kinds of wavevector sets.
 */
class Points {
public:
  /** Default constructor.
   * @param crystal: the crystal object that defines the Brillouin zone.
   * @param mesh: grid size of the Monkhorst-pack.
   * @param offset: the offset of the grid w.r.t. the Gamma point. Offset
   * ranges in values from 0 to 1. 0.5 means half displacement of the grid.
   */
  Points(Crystal &crystal_, const Eigen::Vector3i &mesh_,
         const Eigen::Vector3d &offset_ = Eigen::Vector3d::Zero());

  /** Copy constructor
   */
  Points(const Points &obj);

  /** Copy assignment operator
   */
  Points &operator=(const Points &obj);

  /** Returns the details of the Monkhorst-pack mesh of wavevectors.
   * @return <grid, offset>: grid is the 3D integration mesh used for the
   * Brillouin zone sampling, and the 3D offset of the grid w.r.t. the
   * Gamma point (values of offset from 0 to 1).
   */
  std::tuple<Eigen::Vector3i, Eigen::Vector3d> getMesh();

  /** Returns the number of wavevectors stored in the points class.
   */
  long getNumPoints();

  /** Converts a wavevector from crystal to cartesian coordinates.
   * @param point: the input wavevector in crystal coordinates.
   * @return wavevector: the output wavevector in cartesian coordinates.
   */
  Eigen::Vector3d crystalToCartesian(const Eigen::Vector3d &point);

  /** Converts a wavevector from cartesian to crystal coordinates.
   * @param point: the input wavevector in cartesian coordinates.
   * @return wavevector: the output wavevector in crystal coordinates.
   */
  Eigen::Vector3d cartesianToCrystal(const Eigen::Vector3d &point);

  /** Folds a wavevector in crystal coordinates to the Wigner Seitz zone.
   * @param pointCrystal: the crystal coordinates of a wavevecto
   * @param basis: basis (Points::cartesianCoordinates or
   * Points::crystalCoordinates) in which to return the folded wavevector.
   * @return wavevector: the wavevector coordinates folded in the WS zone.
   */
  Eigen::Vector3d crystalToWS(const Eigen::Vector3d &pointCrystal,
                              const int &basis);

  /** Given a list of points, finds the monkhorst-pack mesh.
   * @param points: a matrix (3,numPoints) of wavevectors in crystal
   * coordinates.
   * @return mesh,offset: the monkhorst-pack grid size and the offset w.r.t.
   * the gamma point, that generates the input points list.
   * Calls an error if the mesh doesn't correspond to a complete list.
   */
  static std::tuple<Eigen::Vector3i, Eigen::Vector3d>
  findMesh(const Eigen::Matrix<double, 3, Eigen::Dynamic> &points);

  /** Returns a reference to the crystal object on which the Points are
   * defined.
   */
  Crystal &getCrystal();

  /** Returns a Point object given its integer index.
   * The Point object is used to move around the code the coordinates of the
   * wavevector.
   * @param index: the integer wavevector index ranging in [0,numPoints[
   * @return point: a point object.
   */
  virtual Point getPoint(const long &index);

  /** Get the coordinates of a wavevector from its index.
   * @param index: the index of the desired wavevector.
   * @param basis: specify the basis to be used for the output coordinates.
   * Either Points::crystalCoords or Points::cartesianCoords.
   * @return wavevector: the coordinates of the desired wavevector.
   */
  virtual Eigen::Vector3d getPointCoords(const long &index,
                                         const int &basis = crystalCoords);

  /** Get the wavevector index given the crystal coordinates of a wavevector.
   * @param point: the wavevector in crystal coordinates.
   * @return index: the index of the wavevector in the range [0,numPoints[
   */
  virtual long getIndex(const Eigen::Vector3d &point);

  /** Returns the value of the k-point weight, i.e. the weight that needs to
   * be used in integrations of the Brillouin zone. Simply a constant if we
   * are not using symmetries.
   * @param ik: the index of the wavevector.
   * @return weight: a weight, normalized to one over the full BZ integration
   */
  virtual double getWeight(const long &ik);

  // note: constexpr tells the compiler that the class member is
  // available at compilation time
  static const int crystalCoords;
  static const int cartesianCoords;

protected:
  void setMesh(const Eigen::Vector3i &mesh_, const Eigen::Vector3d &offset_);
  Crystal &crystal;
  Eigen::Vector3i mesh;
  Eigen::Vector3d offset;
  long numPoints = 0;
  // for Wigner Seitz folding
  Eigen::MatrixXd gVectors;
  Eigen::MatrixXi igVectors;

  // methods to be overwritten
  Eigen::Vector3d reduciblePoints(const long &idx);
};

class FullPoints : public Points {
public:
  /** Default constructor.
   * @param crystal: the crystal object that defines the Brillouin zone.
   * @param mesh: grid size of the Monkhorst-pack.
   * @param offset: the offset of the grid w.r.t. the Gamma point. Offset
   * ranges in values from 0 to 1. 0.5 means half displacement of the grid.
   */
  FullPoints(Crystal &crystal_, const Eigen::Vector3i &mesh_,
             const Eigen::Vector3d &offset_ = Eigen::Vector3d::Zero());

  /** Copy constructor
   */
  FullPoints(const FullPoints &obj); // copy constructor

  /** Copy assignment operator
   */
  FullPoints &operator=(const FullPoints &obj); // assignment operator

  /** Returns a Point object given its integer index.
   * The Point object is used to move around the code the coordinates of the
   * wavevector.
   * @param index: the integer wavevector index ranging in [0,numPoints[
   * @return point: a point object.
   */
  Point getPoint(const long &index);
};

class IrreduciblePoints : public Points {
protected:
  Eigen::VectorXi mapReducibleToIrreducible;
  Eigen::VectorXi mapIrreducibleToReducible;
  // points are internally stored in crystal coordinates
  Eigen::MatrixXd irreduciblePoints;
  Eigen::VectorXd irreducibleWeights;
  Eigen::VectorXi indexIrreduciblePoints;
  long numIrredPoints = 0;

  //	Eigen::Vector3d pointsCoords(const long & index);
  void setIrreduciblePoints();

public:
  IrreduciblePoints(Crystal &crystal_, const Eigen::Vector3i &mesh_,
                    const Eigen::Vector3d &offset_ = Eigen::Vector3d::Zero());
  IrreduciblePoints(const IrreduciblePoints &obj); // copy constructor
  IrreduciblePoints &operator=(const IrreduciblePoints &obj); // assignment

  long getIndex(const Eigen::Vector3d &point);
  Eigen::VectorXi getIndexReducibleFromIrreducible(const long &indexIrr);
  long getIndexIrreducibleFromReducible(const long &indexRed);
  Point getPoint(const long &index);

  long getNumPoints();
  double getWeight(const long &ik);
  Eigen::Vector3d getPointCoords(const long &index,
                                 const int &basis = crystalCoords);
};

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
  Eigen::Vector3d getPointCoords(const long &index,
                                 const int &basis = crystalCoords);
};

/** Class for storing wavevectors along a path of the Brillouin zone.
 */
class PathPoints : public FullPoints {
public:
  /** Default constructor.
   * @param crystal: the crystal object that defines the Brillouin zone.
   * @param pathExtrema: extrema of segments defining a path of points.
   * PathExtrema has size (numSegments, 2, 3), where the first index spans
   * the segments of the path, 2 is the beginning and the last element of the
   * segment, and 3 is the cartesian coordinates of the extrema.
   * Point coordinates must be in crystal coordinates.
   * Example: if we want a path G->X->L, we have: 2 segments: G->X and X->L
   * pathExtrema(0,0,:) = G.coords
   * pathExtrema(0,1,:) = X.coords
   * pathExtrema(1,0,:) = X.coords
   * pathExtrema(1,1,:) = L.coords
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
  PathPoints &operator=(const PathPoints &obj);

  /** Returns a Point object given its integer index.
   * The Point object is used to move around the code the coordinates of the
   * wavevector.
   * @param index: the integer wavevector index ranging in [0,numPoints[
   * @return point: a point object.
   */
  Point getPoint(const long &index);

  /** Get the wavevector index given the crystal coordinates of a wavevector.
   * @param point: the wavevector in crystal coordinates.
   * @return index: the index of the wavevector in the range [0,numPoints[
   */
  long getIndex(const Eigen::Vector3d &coords);

  /** Get the coordinates of a wavevector from its index.
   * @param index: the index of the desired wavevector.
   * @param basis: specify the basis to be used for the output coordinates.
   * Either Points::crystalCoords or Points::cartesianCoords.
   * @return wavevector: the coordinates of the desired wavevector.
   */
  Eigen::Vector3d getPointCoords(const long &index,
                                 const int &basis = crystalCoords);

protected:
  Eigen::Matrix<double, 3, Eigen::Dynamic> pointsList;
};

#endif
