#ifndef EPA_TRANSPORT_APP_H
#define EPA_TRANSPORT_APP_H

#include "app.h"
#include "basevector_bte.h"

class TransportEpaApp : public App {
public:
  void run(Context &context);

private:
  Eigen::Tensor<double, 3>
  calcEnergyProjVelocity(Context &context, BaseBandStructure &bandStructure,
                         const Eigen::VectorXd &energies);
  BaseVectorBTE getScatteringRates(
      Context &context, StatisticsSweep &statisticsSweep,
      FullBandStructure &fullBandStructure, Eigen::VectorXd &energies);
};

#endif