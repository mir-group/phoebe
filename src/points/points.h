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
  //    Eigen::Vector3d crystalToWS(const Eigen::Vector3d &pointCrystal,
  //            const int &basis);

  Eigen::Vector3d bzToWs(const Eigen::Vector3d &point, const int &basis);

  Eigen::Vector3d foldToBz(const Eigen::Vector3d &pointCrystal,
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

  // like getIndex, but returns -1 if point not found
  virtual long isPointStored(const Eigen::Vector3d &crystalCoords);

  // note: constexpr tells the compiler that the class member is
  // available at compilation time
  static const int crystalCoords;
  static const int cartesianCoords;

  // given a wavevector of the reducible list in crystal coordinates,
  // finds the integer index ikIrr of the irreducible point in the irreducible
  // list. Provides also the rotation matrix, in cartesian coordinates, such
  // that rotation * kIrr = kRed
  virtual void setIrreduciblePoints(
      std::vector<Eigen::MatrixXd> *groupVelocities = nullptr);
  std::vector<long> irrPointsIterator();
  std::vector<long> parallelIrrPointsIterator();
  long asIrreducibleIndex(const long &ik);
  long asReducibleIndex(const long &ik);
  virtual std::tuple<long, Eigen::Matrix3d>
  getRotationToIrreducible(const Eigen::Vector3d &x,
                           const int &basis = cartesianCoords);
  virtual std::vector<Eigen::Matrix3d> getRotationsStar(const long &ik);
  virtual std::vector<long> getReduciblesFromIrreducible(const long &ik);

protected:
  void setMesh(const Eigen::Vector3i &mesh_, const Eigen::Vector3d &offset_);
  Crystal &crystal;
  Eigen::Vector3i mesh;
  Eigen::Vector3d offset;
  long numPoints = 0;
  // for Wigner Seitz folding
  Eigen::MatrixXd gVectors;

  // methods to be overwritten
  Eigen::Vector3d reduciblePoints(const long &idx);

  //------------------------------------------
  // Symmetries

  std::vector<Eigen::Matrix3d> rotationMatricesCrystal;
  std::vector<Eigen::Matrix3d> rotationMatricesCartesian;

  // vector of size 0 to numPoints: given index i, mERI(i) is the index of
  // the rotation R such that R*kIrr = k(i)
  Eigen::VectorXi mapEquivalenceRotationIndex;

  // vector of size 0 to numIrrPoints: given an irreducible wavevector i,
  // where i runs from 0 to numIrrPoints,
  // mITRL(i) is the index of the irreducible point in the reducible list
  Eigen::VectorXi mapIrreducibleToReducibleList;

  // vector of size 0 to numPoints. Given a reducible wavevector i,
  // mRTI
  Eigen::VectorXi mapReducibleToIrreducibleList;
  long numIrrPoints = 0;

  // vector of size numIrrPoints. Given an irreducible point i
  // (i from 0,numIrrPoints), gets the list of reducible point indices
  // that are symmetry equivalent to point i.
  std::vector<std::vector<long>> irreducibleStars;

  // if equiv(i) == i, point is irreducible, otherwise, equiv(i) gives the
  // index of the irreducible equivalent index.
  Eigen::VectorXi equiv;
  //------------------------------------------
};

#endif