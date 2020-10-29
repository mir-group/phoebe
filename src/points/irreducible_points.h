#ifndef IRREDPOINTS_H
#define IRREDPOINTS_H

#include "points.h"

class IrreduciblePoints: public Points {
protected:
  Points &parentPoints;

  // these two maps allow to jump from an index on the full list to the
  // corresponding irreducible point in the irreducible list.
  Eigen::VectorXi mapReducibleToIrreducible;
  Eigen::VectorXi mapIrreducibleToReducible;

  // note: irreducibleStars[ikIrr], where ikIrr is an index in
  // [0,numIrredPoints], returns a vector with the integers indexes of ik,
  // indicating which vectors in the full grid belong to this set of equivalent
  // points
  std::vector<std::vector<int>> irreducibleStars;
  Eigen::VectorXi mapEquivalenceRotationIndex;

  // these allow the rotation of points, also in cartesian coordinates
  std::vector<Eigen::Matrix3d> rotationMatricesCrystal;
  std::vector<Eigen::Matrix3d> rotationMatricesCartesian;

  Eigen::MatrixXd pointsList;
  Eigen::VectorXd weights;
public:

  IrreduciblePoints(Points &parentPoints_);
  IrreduciblePoints(const IrreduciblePoints &obj); // copy constructor
  IrreduciblePoints& operator=(const IrreduciblePoints &obj); // assignment


  long getIndex(const Eigen::Vector3d &point);

//  Eigen::VectorXi getIndexReducibleFromIrreducible(const long &indexIrr);
//  long getIndexIrreducibleFromReducible(const long &indexRed);

  Point getPoint(const long &index);

  double getWeight(const long &ik);
  Eigen::Vector3d getPointCoords(const long &index,
      const int &basis=crystalCoords);

  Points getParentPoints();

  // given a wavevector of the reducible list in crystal coordinates,
  // finds the integer index ikIrr of the irreducible point in the irreducible
  // list. Provides also the rotation matrix, in cartesian coordinates, such
  // that rotation * kIrr = kRed
  std::tuple<int,Eigen::Matrix3d> getRotationToIrreducible(
      const Eigen::Vector3d &x, const int & basis=crystalCoords);
  virtual std::vector<Eigen::Matrix3d> getRotationsStar(const int & ik);
  int irreducibleToReducible(const int &ikIrr);
//  int reducibleToIrreducible(const int &ikRed);
};

#endif
