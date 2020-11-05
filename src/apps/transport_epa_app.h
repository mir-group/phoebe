#ifndef EPA_TRANSPORT_APP_H
#define EPA_TRANSPORT_APP_H

#include "app.h"
#include "basevector_bte.h"
#include "delta_function.h"

class TransportEpaApp : public App {
public:
  void run(Context &context);
  void checkRequirements(Context &context);

private:
  Eigen::Tensor<double, 3>
  calcEnergyProjVelocity(Context &context, BaseBandStructure &bandStructure,
                         const Eigen::VectorXd &energies,
                         TetrahedronDeltaFunction &tetrahedrons);
  BaseVectorBTE getScatteringRates(
      Context &context, StatisticsSweep &statisticsSweep,
      FullBandStructure &fullBandStructure, Eigen::VectorXd &energies,
      TetrahedronDeltaFunction &tetrahedrons);
};

#endif