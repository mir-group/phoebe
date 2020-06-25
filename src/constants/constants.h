#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <complex>
#include <math.h>
//
// mathematical constants
//
const double pi = 3.14159265358979323846;
const double sqrtPi = sqrt(3.14159265358979323846);
const double twoPi = 2. * pi;
const double fourPi = 2. * pi;
const double one = 1.;
const double zero = 0.;
const std::complex<double> complexZero = { 0., 0. };
const std::complex<double> complexOne = { 1., 0. };
const std::complex<double> complexI = { 0., 1. };
//
// small numbers
//
const double epsilon6 = 1.e-6;
const double epsilon8 = 1.e-8;
const double epsilon16 = 1.e-16;
//
// Physical constants
//
const double speedLightSi = 299792458;
const double electronSi = 1.602176487e-19;
const double hPlanckSi = 6.62607015e-34;
const double hBarSi = hPlanckSi / twoPi;
const double electronVoltSi = 1.6021766208e-19;
const double hartreeSi = 4.359744650e-18;
const double bohrRadiusSi = 0.52917721067e-10;
const double rydbergSi = hartreeSi / 2.;
const double timeRy = hBarSi / rydbergSi;
const double kBoltzmannSi = 1.38064852e-23;
const double kBoltzmannRy = kBoltzmannSi / rydbergSi;
const double electronMassSi = 9.10938215e-31; // in Kg
const double amuSi = 1.660538782e-27; // in Kg

const double e2 = 2.;

//
// conversion
//
const double temperatureAuToSi = 1. / kBoltzmannRy;
const double energyRyToEv = rydbergSi / electronVoltSi;
const double distanceAngToSi = 1.0e-10;
const double distanceRyToSi = bohrRadiusSi;
const double distanceBohrToCm = bohrRadiusSi * 100.;
const double distanceBohrToMum = bohrRadiusSi * 1000000.;
const double timeRyToFs = timeRy * 1e15 * twoPi;
const double mevToPs = hPlanckSi * 1e15 / electronVoltSi;
const double distanceBohrToAng = bohrRadiusSi / distanceAngToSi;
const double massAmuToRy = amuSi / electronMassSi / 2.;
const double massRyToAmu = 1. / massAmuToRy;

const double TeraHertzAu = hPlanckSi / twoPi / hartreeSi * 1.0e+12;

const double ryToTHz = 1. / TeraHertzAu / 2. / twoPi;
const double ryToCmm1 = 1.0e10 * ryToTHz / speedLightSi;

// velocity tested for phonons
const double velocityRyToSi = distanceRyToSi * rydbergSi / hBarSi;

//
// units for transport coefficients
//
const double thConductivityAuToSi = kBoltzmannSi * rydbergSi / bohrRadiusSi
        / hBarSi;
const double viscosityAuToSi = rydbergSi / hPlanckSi * twoPi
        * pow(bohrRadiusSi, 2);
const double elConductivityAuToSi = pow(electronSi, 2) / hBarSi / bohrRadiusSi;
const double thermopowerAuToSi = -kBoltzmannRy / electronSi * rydbergSi;
const double peltierAuToSi = -rydbergSi / electronSi;

#endif
