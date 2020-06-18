#include <iomanip>
#include "phonon_viscosity.h"
#include "constants.h"

PhononViscosity::PhononViscosity(StatisticsSweep &statisticsSweep_,
        Crystal &crystal_, BaseBandStructure &bandStructure_) :
        Observable(statisticsSweep_, crystal_), bandStructure(bandStructure_) {

    tensordxdxdxd = Eigen::Tensor<double, 5>(numCalcs, dimensionality,
            dimensionality, dimensionality, dimensionality);
    tensordxdxdxd.setZero();
}
;

// copy constructor
PhononViscosity::PhononViscosity(const PhononViscosity &that) :
        Observable(that), bandStructure(that.bandStructure) {
}

// copy assigmnent
PhononViscosity& PhononViscosity::operator =(const PhononViscosity &that) {
    Observable::operator=(that);
    if (this != &that) {
        bandStructure = that.bandStructure;
    }
    return *this;
}

void PhononViscosity::calcRTA(VectorBTE &tau) {
    double norm = 1. / bandStructure.getNumPoints(true)
            / crystal.getVolumeUnitCell(dimensionality);

    auto particle = bandStructure.getParticle();
    tensordxdxdxd.setZero();

    auto excludeIndeces = tau.excludeIndeces;

    for (long is = 0; is < bandStructure.getNumStates(); is++) {
        auto en = bandStructure.getEnergy(is);

        auto vel = bandStructure.getGroupVelocity(is);
        auto q = bandStructure.getWavevector(is);

        // skip the acoustic phonons
//		if ( q.norm()==0. && en<0.1 / ryToCmm1 ) continue;
        if (std::find(excludeIndeces.begin(), excludeIndeces.end(), is)
                != excludeIndeces.end())
            continue;

        for (long iCalc = 0; iCalc < numCalcs; iCalc++) {

            auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
            double temperature = calcStat.temperature;
            double chemPot = calcStat.chemicalPotential;
            double bosep1 = particle.getPopPopPm1(en, temperature, chemPot);

            for (long i = 0; i < dimensionality; i++) {
                for (long j = 0; j < dimensionality; j++) {
                    for (long k = 0; k < dimensionality; k++) {
                        for (long l = 0; l < dimensionality; l++) {
                            tensordxdxdxd(iCalc, i, j, k, l) += q(i) * vel(j)
                                    * q(k) * vel(l) * bosep1
                                    * tau.data(iCalc, is) / temperature * norm;
                        }
                    }
                }
            }
        }
    }
}

void PhononViscosity::calcFromRelaxons(Vector0 &vector0, VectorBTE &relTimes,
        PhScatteringMatrix &sMatrix, Eigen::MatrixXd &eigenvectors) {

    if (numCalcs > 1) {
        Error e("Viscosity for relaxons only for 1 temperature");
    }

    // we decide to skip relaxon states
    // 1) there is a relaxon with zero (or epsilon) eigenvalue -> infinite tau
    // 2) if we include (3) acoustic modes at gamma, we have 3 zero eigenvalues
    //    because we set some matrix rows/cols to zero
    long firstState = 1;
    firstState += relTimes.excludeIndeces.size();

    double volume = crystal.getVolumeUnitCell(dimensionality);
    double numPoints = double(bandStructure.getNumPoints(true));
    long numStates = bandStructure.getNumStates();
    auto particle = bandStructure.getParticle();

    // to simplify, here I do everything considering there is a single
    // temperature (due to memory constraints)

    Eigen::VectorXd A(dimensionality);
    A.setZero();

    long iCalc = 0;
    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStat.temperature;
    double chemPot = calcStat.chemicalPotential;

    for (long is = firstState; is < numStates; is++) {
        auto en = bandStructure.getEnergy(is);
        double bosep1 = particle.getPopPopPm1(en, temp, chemPot); // = n(n+1)
        auto q = bandStructure.getWavevector(is);
        for (long idim = 0; idim < dimensionality; idim++) {
            A(idim) += bosep1 * q(idim) * q(idim);
        }
    }
    A /= temp * numPoints * volume;

    Eigen::MatrixXd driftEigenvector(3, numStates);
    driftEigenvector.setZero();
    for (long is = firstState; is < numStates; is++) {
        auto en = bandStructure.getEnergy(is);
        double bosep1 = particle.getPopPopPm1(en, temp, chemPot); // = n(n+1)
        auto q = bandStructure.getWavevector(is);
        for (auto idim : { 0, 1, 2 }) {
            driftEigenvector(idim, is) = q(idim)
                    * sqrt(bosep1 / temp / A(idim));
        }
    }

    Eigen::MatrixXd D(3, 3);
    D = driftEigenvector * sMatrix.dot(driftEigenvector.transpose());
    D /= volume * numPoints;

    Eigen::Tensor<double, 3> tmpDriftEigvecs(3, 3, numStates);
    tmpDriftEigvecs.setZero();
    Eigen::MatrixXd W(3, 3);
    W.setZero();
    for (long is = firstState; is < bandStructure.getNumStates(); is++) {
        auto v = bandStructure.getGroupVelocity(is);
        for (auto i : { 0, 1, 2 }) {
            for (auto j : { 0, 1, 2 }) {
                tmpDriftEigvecs(i, j, is) = driftEigenvector(j, is) * v(i);
                W(i, j) += vector0.data(0, is) * v(i) * driftEigenvector(j, is);
            }
        }
    }
    W /= volume * numPoints;

    Eigen::Tensor<double, 3> w(3, 3, numStates);
    w.setZero();
    for (auto i : { 0, 1, 2 }) {
        for (auto j : { 0, 1, 2 }) {
            Eigen::VectorXd x(numStates);
            for (long is = 0; is < numStates; is++) {
                x(is) = tmpDriftEigvecs(i, j, is);
            }
            auto x2 = x.transpose() * eigenvectors;
            for (long is = 0; is < numStates; is++) {
                w(i, j, is) = x2(is) / volume / numPoints;
            }
        }
    }

    // Eq. 9, Simoncelli PRX (2019)
    tensordxdxdxd.setZero();
    for (long is = firstState; is < numStates; is++) {
        for (long i = 0; i < dimensionality; i++) {
            for (long j = 0; j < dimensionality; j++) {
                for (long k = 0; k < dimensionality; k++) {
                    for (long l = 0; l < dimensionality; l++) {
                        tensordxdxdxd(iCalc, i, j, k, l) += 0.5
                                * (w(i, j, is) * w(k, l, is)
                                        + w(i, l, is) * w(k, j, is)) * A(i)
                                * A(k) * relTimes.data(0, is);
                    }
                }
            }
        }
    }
}

void PhononViscosity::print() {
    std::string units;
    if (dimensionality == 1) {
        units = "Pa s / m^2";
    } else if (dimensionality == 2) {
        units = "Pa s / m";
    } else {
        units = "Pa s";
    }

    std::cout << "\n";
    std::cout << "Thermal Viscosity (" << units << ")\n";
    std::cout << "i, j, k, eta[i,j,k,1,0], eta[i,j,k,1], eta[i,j,k,2]\n";

    for (long iCalc = 0; iCalc < numCalcs; iCalc++) {

        double conversion = pow(hBarSi, 2) // momentum is hbar q
        / pow(distanceRyToSi, dimensionality) // volume conversion
        * rydbergSi / hBarSi // conversion time (q^2 v^2 tau = [time])
                / rydbergSi; // temperature conversion

        auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
        double temp = calcStat.temperature;

        std::cout << std::fixed;
        std::cout.precision(2);
        std::cout << "Temperature: " << temp * temperatureAuToSi << " (K)\n";
        std::cout.precision(5);
        std::cout << std::scientific;
        for (long i = 0; i < dimensionality; i++) {
            for (long j = 0; j < dimensionality; j++) {
                for (long k = 0; k < dimensionality; k++) {
                    std::cout << i << " " << j << " " << k;
                    for (long l = 0; l < dimensionality; l++) {
                        std::cout << " " << std::setw(12) << std::right
                                << tensordxdxdxd(iCalc, i, j, k, l)
                                        * conversion;
                    }
                    std::cout << "\n";
                }
            }
        }
        std::cout << "\n";
    }
}

int PhononViscosity::whichType() {
    return is4Tensor;
}
