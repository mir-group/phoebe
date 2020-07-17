#ifndef IRREDPOINTS_H
#define IRREDPOINTS_H

#include "points.h"

class IrreduciblePoints: public Points {
protected:
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

  std::vector<Eigen::Matrix3d> rotationMatrices; //in cartesian reciprocal space

  Eigen::MatrixXd pointsList;
  Points &parentPoints;
  Eigen::VectorXd weights;
public:

  IrreduciblePoints(Points &parentPoints_);
  IrreduciblePoints(const IrreduciblePoints &obj); // copy constructor
  IrreduciblePoints& operator=(const IrreduciblePoints &obj); // assignment

  long getIndex(const Eigen::Vector3d &point);

  Eigen::VectorXi getIndexReducibleFromIrreducible(const long &indexIrr);
  long getIndexIrreducibleFromReducible(const long &indexRed);

  Point getPoint(const long &index);

  long getNumPoints();
  double getWeight(const long &ik);
  Eigen::Vector3d getPointCoords(const long &index,
      const int &basis=crystalCoords);

  Points getParentPoints();
};

#endif
