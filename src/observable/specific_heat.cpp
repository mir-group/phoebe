#include "specific_heat.h"
#include "constants.h"

SpecificHeat::SpecificHeat(StatisticsSweep & statisticsSweep_,
		Crystal & crystal_, BaseBandStructure & bandStructure_) :
				Observable(statisticsSweep_, crystal_),
				bandStructure(bandStructure_) {
	scalar = Eigen::VectorXd::Zero(numCalcs);
};

// copy constructor
SpecificHeat::SpecificHeat(const SpecificHeat & that) : Observable(that),
				bandStructure(that.bandStructure) {
}

// copy assignment
SpecificHeat & SpecificHeat::operator = (const SpecificHeat & that) {
	Observable::operator=(that);
    if ( this != &that) {
    	bandStructure = that.bandStructure;
    }
    return *this;
}

void SpecificHeat::calc() {
	double norm = 1. / bandStructure.getNumPoints()
			/ crystal.getVolumeUnitCell(dimensionality);
	scalar.setZero();
	auto particle = bandStructure.getParticle();
	for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {
		auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
		double temp = calcStat.temperature;
		double chemPot = calcStat.chemicalPotential;
		for ( long is=0; is<bandStructure.getNumStates(); is++ ) {
			auto en = bandStructure.getEnergy(is);
			auto dndt = particle.getDndt(en, temp, chemPot);
			scalar(iCalc) += dndt * en * norm;
		}
	}
}

void SpecificHeat::print() {
	std::string units;
	if ( dimensionality == 1 ) {
		units = "J / K / m";
	} else if ( dimensionality == 2 ) {
		units = "J / K / m^2";
	} else {
		units = "J / K / m^3";
	}

	double conversion = kBoltzmannSi / pow(bohrRadiusSi,3);

	std::cout << "\n";
	std::cout << "Specific heat (" << units << ")\n";

	for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {

		auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
		double temp = calcStat.temperature;

		std::cout << std::fixed;
		std::cout.precision(2);
		std::cout << "Temperature: " << temp * temperatureAuToSi
				<< " (K), C = ";
		std::cout << std::scientific;
		std::cout.precision(5);
		std::cout << scalar(iCalc)*conversion;
		std::cout << "\n";
	}
}

int SpecificHeat::whichType() {
	return isScalar;
}

const double & SpecificHeat::get(const ChemPotIndex & imu, const TempIndex & it) {
	auto i = glob2Loc(imu,it);
	return scalar(i);
}




