#ifndef INTERACTION_EPA_H
#define INTERACTION_EPA_H

#include "context.h"
#include "eigen.h"

class InteractionEpa {

private:
  int numBandGroups;
  Eigen::VectorXd bandExtrema;
  Eigen::VectorXd binSize;
  Eigen::VectorXi numBins;
  Eigen::VectorXd phFreqAverage;
  Eigen::Tensor<double, 4> elPhMatAverage;

public:
  // default constructor
  InteractionEpa(int &numBandGroups_, Eigen::VectorXd &bandExtrema_,
                 Eigen::VectorXd &binSize_, Eigen::VectorXi &numBins_,
                 Eigen::VectorXd &phFreqAverage_,
                 Eigen::Tensor<double, 4> &elPhMatAverage_);

  // copy constructor
  InteractionEpa(const InteractionEpa &that);

  // assignment operator overload
  InteractionEpa &operator=(const InteractionEpa &that);

  int getNumBandGroups();
  Eigen::VectorXd getBandExtrema();
  Eigen::VectorXd getBinSize();
  Eigen::VectorXi getNumBins();
  Eigen::VectorXd getPhFreqAverage();
  Eigen::Tensor<double, 4> getElPhMatAverage();

  static InteractionEpa parseEpaCoupling(Context &context);
};

#endif
