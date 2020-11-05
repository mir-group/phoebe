#ifndef INTERACTION_EPA_H
#define INTERACTION_EPA_H

#include "context.h"
#include "eigen.h"

class InteractionEpa {

private:
  Eigen::VectorXd elEnergies;
  Eigen::VectorXd phEnergies;
  Eigen::Tensor<double, 3> elPhMatAverage;

public:
  // default constructor
  InteractionEpa(Eigen::VectorXd &elEnergies_,
                 Eigen::VectorXd &phEnergies_,
                 Eigen::Tensor<double, 3> &elPhMatAverage_);

  // copy constructor
  InteractionEpa(const InteractionEpa &that);

  // assignment operator overload
  InteractionEpa &operator=(const InteractionEpa &that);

  Eigen::VectorXd getPhEnergies();
  Eigen::VectorXd getElEnergies();
  Eigen::Tensor<double, 3> getElPhMatAverage();

  static InteractionEpa parseEpaCoupling(Context &context);
};

#endif
