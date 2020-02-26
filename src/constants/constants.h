#include <complex>
#include <math.h>
using namespace std;

namespace constants {
//
// mathematical constants
//
const double pi = 3.14159265358979323846;
const double twoPi = 2. * pi;
const double one = 1.;
const double zero = 0.;
const complex <double> complexZero (0.,0.);
const complex <double> complexOne (1.,0.);
const complex <double> complexI (0.,1.);
//
// small numbers
//
const double epsilon6 = 1.e-6;
const double epsilon8 = 1.e-8;
const double epsilon16 = 1.e-16;
//
// Physical constants
//
const double speedLightSi   = 299792458;
const double electronSi     = 1.602176487e-19;
const double hPlanckSi      = 6.62607015e-34;
const double hBarSi         = hPlanckSi / twoPi;
const double electronVoltSi = 1.6021766208e-19;
const double hartreeSi      = 4.359744650e-18;
const double bohrRadiusSi   = 0.52917721067e-10;
const double rydbergSi      = hartreeSi / 2.;
const double timeRy         = hBarSi / rydbergSi;
const double kBoltzmannSi   = 1.38064852e-23;
const double kBoltzmannRy   = kBoltzmannSi / rydbergSi;

const double temperatureAuToSi = 1. / kBoltzmannRy;
const double energyAuToEv   = rydbergSi / electronVoltSi;
const double distanceAngToSi = 1.0e-10;
const double distanceAuToSi = bohrRadiusSi;
const double timeAuToFs = timeRy * 1e15 * twoPi;
const double mevToPs = hPlanckSi * 1e15 / electronVoltSi;
const double distanceAuToAng = bohrRadiusSi / distanceAngToSi;
const double ryToCmm1 = 0.01 / ( hPlanckSi / electronVoltSi * speedLightSi );
//
// units for transport coefficients
//
const double thConductivityAuToSi = kBoltzmannSi * rydbergSi / bohrRadiusSi / hBarSi;
const double viscosityAuToSi = rydbergSi / hPlanckSi * twoPi * pow(bohrRadiusSi,2);
const double elConductivityAuToSi = pow(electronSi,2) / hBarSi / bohrRadiusSi;
const double thermopowerAuToSi = - kBoltzmannRy / electronSi * rydbergSi;
const double peltierAuToSi = - rydbergSi / electronSi;

}

