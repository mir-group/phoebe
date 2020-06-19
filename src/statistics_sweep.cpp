#include <algorithm>
#include "utilities.h"
#include "statistics_sweep.h"
#include "bandstructure.h"
#include "context.h"
#include "constants.h"
#include <cmath>

StatisticsSweep::StatisticsSweep(Context &context,
        FullBandStructure *fullBandStructure) :
        particle(
                fullBandStructure != nullptr ?
                        fullBandStructure->getParticle() : Particle::phonon) {

    if (particle.isPhonon()) {
        Eigen::VectorXd temperatures = context.getTemperatures();
        nTemp = temperatures.size();
        nDop = 1;
        nChemPot = 1;
        numCalcs = nTemp * std::max(nChemPot, nDop);
        infoCalcs = Eigen::MatrixXd::Zero(numCalcs, 3);
        for (long it = 0; it < nTemp; it++) {
            double temp = temperatures(it);
            infoCalcs(it, 0) = temp;
            // note: the other two columns are set to zero above
        }

    } else { // isElectron()
        // in this case, the class is more complicated, as it is tasked with
        // the calculation of the chemical potential or the doping
        // so, we first load parameters to compute these quantities

        bool hasSpinOrbit = context.getHasSpinOrbit();
        if (hasSpinOrbit) {
            spinFactor = 1.;
        } else { // count spin degeneracy
            spinFactor = 2.;
        }

        volume =
                fullBandStructure->getPoints().getCrystal().getVolumeUnitCell();

        // flatten the energies (easier to work with)
        numPoints = fullBandStructure->getNumPoints();
        numBands = fullBandStructure->getNumBands();
        numStates = numBands * numPoints;
        energies = Eigen::VectorXd::Zero(numBands * numPoints);
        long i = 0;
        for (long ik = 0; ik < numPoints; ik++) {
            auto state = fullBandStructure->getState(ik);
            auto ens = state.getEnergies();
            for (long ib = 0; ib < numBands; ib++) {
                energies(i) = ens(ib);
                i++;
            }
        }

        // determine ground state properties
        // i.e. we define the number of filled bands and the fermi energy

        occupiedStates = context.getNumOccupiedStates();
        if (std::isnan(occupiedStates)) {
            // in this case we try to compute it from the Fermi-level
            fermiLevel = context.getFermiLevel();
            if (std::isnan(fermiLevel)) {
                Error e("Must provide either the Fermi level or the number of"
                        " occupied states");
            }
            occupiedStates = 0.;
            for (long i = 0; i < numStates; i++) {
                if (energies(i) < fermiLevel) {
                    occupiedStates += 1.;
                }
            }
            occupiedStates /= numPoints;
        } else {
            occupiedStates /= spinFactor;
            fermiLevel = energies.mean(); // first guess of Fermi level
            // refine this guess at zero temperature and zero doping
            fermiLevel = findChemicalPotentialFromDoping(0., 0.);
        }

        // load input sweep parameters

        Eigen::VectorXd temperatures = context.getTemperatures();
        Eigen::VectorXd dopings = context.getDopings();
        Eigen::VectorXd chemicalPotentials = context.getChemicalPotentials();

        nChemPot = chemicalPotentials.size();
        nDop = dopings.size();
        nTemp = temperatures.size();
        if ((nDop == 0) && (nChemPot == 0)) {
            Error e("Didn't find chemical potentials or doping in input");
        }

        // now have two cases:
        // if doping is passed in input, compute chemical potential
        // or viceversa

        Eigen::MatrixXd calcTable(nTemp * std::max(nChemPot, nDop), 3);
        calcTable.setZero();

        // in this case, I have dopings, and want to find chemical potentials
        if (nChemPot == 0) {
            nChemPot = nDop;

            for (long it = 0; it < nTemp; it++) {
                for (long id = 0; id < nDop; id++) {
                    double temp = temperatures(it);
                    double doping = dopings(id);
                    double chemPot = findChemicalPotentialFromDoping(doping,
                            temp);

                    long iCalc = compress2Indeces(it, id, nTemp, nDop);
                    calcTable(iCalc, 0) = temp;
                    calcTable(iCalc, 1) = chemPot;
                    calcTable(iCalc, 2) = doping;
                }
            }

            // in this case, I have chemical potentials
        } else if (nDop == 0) {
            nDop = nChemPot;

            for (long it = 0; it < nTemp; it++) {
                for (long imu = 0; imu < nChemPot; imu++) {
                    double temp = temperatures(it);
                    double chemPot = chemicalPotentials(imu);
                    double doping = findDopingFromChemicalPotential(chemPot,
                            temp);

                    long iCalc = compress2Indeces(it, imu, nTemp, nChemPot);
                    calcTable(iCalc, 0) = temp;
                    calcTable(iCalc, 1) = chemPot;
                    calcTable(iCalc, 2) = doping;
                }
            }
        }

        // clean memory
        energies = Eigen::VectorXd::Zero(1);
        // save potentials and dopings
        infoCalcs = calcTable;
        numCalcs = calcTable.rows();
    }
}

// copy constructor
StatisticsSweep::StatisticsSweep(const StatisticsSweep &that) :
        particle(that.particle) {
    numCalcs = that.numCalcs;
    infoCalcs = that.infoCalcs;
    nTemp = that.nTemp;
    nChemPot = that.nChemPot;
    nDop = that.nDop;
}

// copy assignment
StatisticsSweep& StatisticsSweep::operator =(const StatisticsSweep &that) {
    if (this != &that) {
        particle = that.particle;
        numCalcs = that.numCalcs;
        infoCalcs = that.infoCalcs;
        nTemp = that.nTemp;
        nChemPot = that.nChemPot;
        nDop = that.nDop;
    }
    return *this;
}

double StatisticsSweep::fPop(const double &chemPot, const double &temp) {
    // fPop = 1/NK \sum_\mu FermiDirac(\mu) - N
    // Note that I don`t normalize the integral, which is the same thing I did
    // for computing the particle number
    double fPop_ = 0.;
    for (long i = 0; i < numStates; i++) {
        fPop_ += particle.getPopulation(energies(i), temp, chemPot);
    }
    fPop_ /= numPoints;
    fPop_ = numElectronsDoped - fPop_;
    return fPop_;
}

double StatisticsSweep::findChemicalPotentialFromDoping(const double &doping,
        const double &temperature) {
    // given the carrier concentration, finds the fermi energy
    // temperature is set inside glob
    // To find fermi energy, I must find the root of \sum_s f(s) - N = 0
    // the root is found with a bisection algorithm
    // Might be numerically unstable for VERY small doping concentration

    // numElectronsDoped is the total number of electrons in the unit cell
    // numElectrons is the number of electrons in the unit cell before doping
    // doping > 0 means p-doping (fermi level in the valence band)

    numElectronsDoped = occupiedStates
            - doping * volume * pow(distanceBohrToCm, 3) / spinFactor;

    // bisection method: I need to find the root of N - \int fermi dirac = 0

    // initial guess
    double chemicalPotential = fermiLevel;

    // I choose the following  (generous) boundaries
    double aX = chemicalPotential - 0.25; // - 3.5eV
    double bX = chemicalPotential + 0.25; // + 3.5eV
    double aY = fPop(aX, temperature);
    double bY = fPop(bX, temperature);

    if (sgn(aY) == sgn(bY)) {
        Error e("I should revisit the boundary limits for bisection method");
    }

    for (long iter = 0; iter < maxIter; iter++) {
        if (iter == maxIter - 1) {
            Error e("Max iteration reached in finding mu");
        }
        double cX = (aX + bX) / 2.;
        double cY = fPop(cX, temperature);

        // exit condition: the guess is exact or didn't change much
        if ((cY == 0.) || (abs(bX - aX) < 1.0e-8)) {
            chemicalPotential = cX;
            break; // go out of the loop
        }

        // check the sign
        if (sgn(cY) == sgn(aY)) {
            aX = cX;
        } else {
            bX = cX;
        }
    }
    return chemicalPotential;
}

double StatisticsSweep::findDopingFromChemicalPotential(
        const double &chemicalPotential, const double &temperature) {
    double fPop = 0.;
    for (long i = 0; i < numStates; i++) {
        fPop += particle.getPopulation(energies(i), temperature,
                chemicalPotential);
    }
    fPop /= numPoints;
    double doping = (occupiedStates - fPop);
    doping *= spinFactor / volume / pow(distanceBohrToCm, 3);
    return doping;
}

struct CalcStatistics StatisticsSweep::getCalcStatistics(const long &index) {
    struct CalcStatistics sc;
    sc.temperature = infoCalcs(index, 0);
    sc.chemicalPotential = infoCalcs(index, 1);
    sc.doping = infoCalcs(index, 2);
    return sc;
}

struct CalcStatistics StatisticsSweep::getCalcStatistics(const TempIndex &iTemp,
        const ChemPotIndex &iChemPot) {
    long index = compress2Indeces(iTemp.get(), iChemPot.get(), nTemp, nChemPot);
    return getCalcStatistics(index);
}

long StatisticsSweep::getNumCalcs() {
    return numCalcs;
}

long StatisticsSweep::getNumChemicalPotentials() {
    return nChemPot;
}

long StatisticsSweep::getNumTemperatures() {
    return nTemp;
}
