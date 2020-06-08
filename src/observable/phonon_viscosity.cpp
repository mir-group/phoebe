#include "phonon_viscosity.h"
#include "constants.h"
#include <iomanip>

PhononViscosity::PhononViscosity(StatisticsSweep & statisticsSweep_,
		Crystal & crystal_, FullBandStructure<FullPoints> & bandStructure_) :
				Observable(statisticsSweep_, crystal_),
				bandStructure(bandStructure_) {

	tensordxdxdxd = Eigen::Tensor<double,5>(numCalcs,dimensionality,
			dimensionality,dimensionality,dimensionality);
	tensordxdxdxd.setZero();
};

// copy constructor
PhononViscosity::PhononViscosity(const PhononViscosity & that) :
		Observable(that), bandStructure(that.bandStructure) {
}

// copy assigmnent
PhononViscosity & PhononViscosity::operator = (const PhononViscosity & that) {
	Observable::operator=(that);
    if ( this != &that) {
    	bandStructure = that.bandStructure;
    }
    return *this;
}

void PhononViscosity::calcRTA(VectorBTE & tau) {
	double norm = 1. / bandStructure.getNumPoints()
			/ crystal.getVolumeUnitCell(dimensionality);

	auto statistics = bandStructure.getStatistics();
	tensordxdxdxd.setZero();

	for ( long ik=0; ik<bandStructure.getNumPoints(); ik++ ) {
		auto s = bandStructure.getState(ik);
		auto ens = s.getEnergies();
		auto vel = s.getVelocities();
		auto q = s.getCoords(Points::cartesianCoords);

		for ( long ib=0; ib<ens.size(); ib++ ) {
			long is = bandStructure.getIndex(ik,ib);
			double en = ens(ib);

			// skip the acoustic phonons
			if ( q.norm()==0. && ib<3 ) continue;

			for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {

				auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
				double temperature = calcStat.temperature;
				double chemPot = calcStat.chemicalPotential;
				double bose = statistics.getPopulation(en,temperature,chemPot);

				for ( long i=0; i<dimensionality; i++ ) {
					for ( long j=0; j<dimensionality; j++ ) {
						for ( long k=0; k<dimensionality; k++ ) {
							for ( long l=0; l<dimensionality; l++ ) {
								tensordxdxdxd(iCalc,i,j,k,l) +=
										q(i)
										* vel(ib,ib,j).real()
										* q(k)
										* vel(ib,ib,l).real()
										* bose * ( 1. + bose )
										* tau.data(iCalc,is)
										/ temperature
										* norm;
							}
						}
					}
				}
			}
		}
	}
}

void PhononViscosity::print() {
	std::string units;
	if ( dimensionality == 1 ) {
		units = "Pa s / m^2";
	} else if ( dimensionality == 2 ) {
		units = "Pa s / m";
	} else {
		units = "Pa s";
	}

	std::cout << "\n";
	std::cout << "Thermal Viscosity (" << units << ")\n";
	std::cout << "i, j, k, eta[i,j,k,1,0], eta[i,j,k,1], eta[i,j,k,2]\n";

	for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {

		double conversion = pow(hBarSi,2) // momentum is hbar q
				/ pow(distanceRyToSi,dimensionality) // volume conversion
				* rydbergSi / hBarSi // conversion time (q^2 v^2 tau = [time])
				/ rydbergSi; // temperature conversion

		auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
		double temp = calcStat.temperature;

		std::cout.precision(5);
		std::cout << "Temperature: " << temp * temperatureAuToSi << " (K)\n";
		std::cout << std::scientific;
		for ( long i=0; i<dimensionality; i++ ) {
			for ( long j=0; j<dimensionality; j++ ) {
				for ( long k=0; k<dimensionality; k++ ) {
					std::cout << i << " " << j << " " << k;
					for ( long l=0; l<dimensionality; l++ ) {
						std::cout << " " << std::setw(12) << std::right <<
								tensordxdxdxd(iCalc,i,j,k,l) * conversion;
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
