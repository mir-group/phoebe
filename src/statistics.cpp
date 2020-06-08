#include "statistics.h"
#include "exceptions.h"

Statistics::Statistics(int particle_) {
	if ( particle_ == phonon ) {
		statistics = bose;
		particle = phonon;
	} else if ( particle_ == electron ) {
		statistics = fermi;
		particle = electron;
	} else {
		Error e("Wrong initialization of Statistics", 1);
	}
}

// copy constructor
Statistics::Statistics(const Statistics & obj) {
	statistics = obj.statistics;
	particle = obj.particle;
	temperatures = obj.temperatures;
	chemicalPotentials = obj.chemicalPotentials;
	dopings = obj.dopings;
}

// copy assignment operator
Statistics & Statistics::operator=(const Statistics & obj) {
	if ( this != &obj ) {
		statistics = obj.statistics;
		particle = obj.particle;
		temperatures = obj.temperatures;
		chemicalPotentials = obj.chemicalPotentials;
		dopings = obj.dopings;
	}
	return *this;
}

Statistics::~Statistics() {
}

bool Statistics::isFermi() {
	if ( statistics == fermi ) {
		return true;
	} else {
		return false;
	}
}

bool Statistics::isBose() {
	if ( statistics == bose ) {
		return true;
	} else {
		return false;
	}
}

bool Statistics::isElectron() {
	if ( particle == electron ) {
		return true;
	} else {
		return false;
	}
}

bool Statistics::isPhonon() {
	if ( particle == phonon ) {
		return true;
	} else {
		return false;
	}
}

double Statistics::getPopulation(const double & energy,
		const double & temperature, const double & chemicalPotential) {
	double population = 0.;
	double y = ( energy - chemicalPotential ) / temperature;

	if ( statistics == bose ) {
		population = 1. / ( exp(y) - 1 );
		if ( population < 0. ) {
			population = 0.;
		}
	} else {
		population = 1. / ( exp(y) + 1 );

		if ( population < 0. ) {
			population = 0.;
		}
		if ( population > 1. ) {
			population = 1.;
		}
	}
	return population;
}

double Statistics::getDndt(const double & energy, const double & temperature,
		const double & chemicalPotential) {
	double x = getPopPopPm1(energy, temperature, chemicalPotential);
	double y = energy - chemicalPotential;
	double dndt = x * y / temperature / temperature;
	return dndt;
}

double Statistics::getDnde(const double & energy, const double & temperature,
		const double & chemicalPotential) {
	double x = getPopPopPm1(energy, temperature, chemicalPotential);
	double dnde = - x / temperature;
	return dnde;
}

double Statistics::getPopPopPm1(const double & energy,
		const double & temperature, const double & chemicalPotential) {
	double x = energy - chemicalPotential;
	double arg = x / 2. / temperature;
	if ( statistics == bose ) {
		x = sinh(arg);
	} else {
		x = cosh(arg);
	}
	x = 0.25 / x / x;
	return x;
}

Eigen::VectorXd Statistics::getDndt(const Eigen::VectorXd & energies,
		const double & temperature, const double & chemicalPotential) {
	double x;
	long numBands = energies.size();
	Eigen::VectorXd dndt(numBands);
	for ( long ib=0; ib<numBands; ib++ ) {
		x = getDndt(energies(ib), temperature, chemicalPotential);
		dndt(ib) = x;
	}
	return dndt;
}

Eigen::VectorXd Statistics::getDnde(const Eigen::VectorXd & energies,
		const double & temperature, const double & chemicalPotential) {
	double x;
	long numBands = energies.size();
	Eigen::VectorXd dnde(numBands);
	for ( long ib=0; ib<numBands; ib++ ) {
		x = getDnde(energies(ib), temperature, chemicalPotential);
		dnde(ib) = x;
	}
	return dnde;
}
