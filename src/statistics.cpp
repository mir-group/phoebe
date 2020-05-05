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
Statistics::Statistics(const Statistics &obj) {
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

double Statistics::getPopulation(double energy, double temperature,
		double chemicalPotential) {
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

double Statistics::getDndt(double energy, double temperature,
		double chemicalPotential) {
	double arg, x, y;
	y = energy - chemicalPotential;
	arg = y / 2. / temperature;
	if ( statistics == bose ) {
		x = sinh(arg); // n(n+1) = 1/4/sinh^2
	} else {
		x = cosh(arg); // n(1-n) = 1/4/cosh^2
	}
	double dndt = y / temperature / temperature / 4. / x / x;
	return dndt;
}

double Statistics::getDnde(double energy, double temperature,
		double chemicalPotential) {
	double arg, x, y;
	y = energy - chemicalPotential;
	arg = y / 2. / temperature;
	if ( statistics == bose ) {
		x = sinh(arg);
	} else {
		x = cosh(arg);
	}
	double dnde = - 1. / temperature / 4. / x / x;
	return dnde;
}

Eigen::VectorXd Statistics::getDndt(Eigen::VectorXd energies,
		double temperature, double chemicalPotential) {
	double x;
	int numTemp = temperatures.size();
	int numEn = temperatures.size();
	Eigen::MatrixXd dndt(numEn, numTemp);
	for ( int ib=0; ib<numEn; ib++ ) {
		for ( int it=0; it<numTemp; it++ ) {
			x = getDndt(energies(ib), temperatures(it), chemicalPotential);
			dndt(ib,it) = x;
		}
	}
	return dndt;
}

Eigen::VectorXd Statistics::getDnde(Eigen::VectorXd energies,
		double temperature, double chemicalPotential) {
	double x;
	int numTemp = temperatures.size();
	int numEn = temperatures.size();
	Eigen::MatrixXd dnde(numEn, numTemp);
	for ( int ib=0; ib<numEn; ib++ ) {
		for ( int it=0; it<numTemp; it++ ) {
			x = getDnde(energies(ib), temperatures(it), chemicalPotential);
			dnde(ib,it) = x;
		}
	}
	return dnde;
}
