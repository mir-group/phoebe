#include "phonon_thermal_cond.h"
#include "constants.h"
#include <time.h>
#include <iomanip>

PhononThermalConductivity::PhononThermalConductivity(
		StatisticsSweep & statisticsSweep_,
		Crystal & crystal_, FullBandStructure & bandStructure_) :
				Observable(statisticsSweep_, crystal_),
				bandStructure(bandStructure_) {

	tensordxd = Eigen::Tensor<double,3>(numCalcs,dimensionality,dimensionality);
	tensordxd.setZero();
};

// copy constructor
PhononThermalConductivity::PhononThermalConductivity(
		const PhononThermalConductivity & that) : Observable(that),
				bandStructure(that.bandStructure) {
}

// copy assigmnent
PhononThermalConductivity & PhononThermalConductivity::operator = (
		const PhononThermalConductivity & that) {
	Observable::operator=(that);
    if ( this != &that) {
    	bandStructure = that.bandStructure;
    }
    return *this;
}

PhononThermalConductivity PhononThermalConductivity::operator - (
		const PhononThermalConductivity & that) {
	PhononThermalConductivity newObservable(statisticsSweep, crystal,
			bandStructure);
	baseOperatorMinus(newObservable, that);
	return newObservable;
}

void PhononThermalConductivity::calcFromCanonicalPopulation(VectorBTE & f) {
	VectorBTE n = f;
	n.canonical2Population(); // n = bose (bose+1) f
	calcFromPopulation(n);
}

void PhononThermalConductivity::calcFromPopulation(VectorBTE & n) {
	double norm = 1. / bandStructure.getNumPoints()
			/ crystal.getVolumeUnitCell(dimensionality);

	tensordxd.setZero();

	for ( long ik=0; ik<bandStructure.getNumPoints(); ik++ ) {
		auto s = bandStructure.getState(ik);
		auto en = s.getEnergies();
		auto vel = s.getVelocities();
		for ( long ib=0; ib<en.size(); ib++ ) {
			long is = bandStructure.getIndex(WavevectorIndex(ik),
					BandIndex(ib));

			// skip the acoustic phonons
			if ( s.getCoords(Points::crystalCoords).norm()==0. && ib<3 ) continue;

			for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {
				auto [imu,it] = loc2Glob(iCalc);

				for ( long i=0; i<dimensionality; i++ ) {
					long icPop = n.glob2Loc(imu,it,DimIndex(i));
					for ( long j=0; j<dimensionality; j++ ) {
							tensordxd(iCalc,i,j) +=
									n.data(icPop,is) * vel(ib,ib,j).real()
									* en(ib) * norm;
					}
				}
			}
		}
	}
}

void PhononThermalConductivity::calcVariational(VectorBTE & af, VectorBTE & f,
		VectorBTE & scalingCG) {
	double norm = 1. / bandStructure.getNumPoints()
			/ crystal.getVolumeUnitCell(dimensionality);

	auto fUnscaled = f;
	fUnscaled = fUnscaled / scalingCG;

	calcFromCanonicalPopulation(fUnscaled);

	tensordxd *= tensordxd.constant(2.);

	for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {
		auto [imu,it] = loc2Glob(iCalc);

		for ( long i=0; i<dimensionality; i++ ) {
			long icPop1 = f.glob2Loc(imu,it,DimIndex(i));
			for ( long j=0; j<dimensionality; j++ ) {
				long icPop2 = f.glob2Loc(imu,it,DimIndex(j));
				for ( long is=0; is<bandStructure.getNumStates(); is++ ) {
					tensordxd(iCalc,i,j) -=
							f.data(icPop1,is) * af.data(icPop2,is) * norm;
				}
			}
		}
	}
}

void PhononThermalConductivity::calcFromRelaxons(SpecificHeat & specificHeat,
		VectorBTE & relaxonV, VectorBTE & relaxationTimes) {

	// we decide to skip relaxon states
	// 1) there is a relaxon with zero (or epsilon) eigenvalue -> infinite tau
	// 2) if we include (3) acoustic modes at gamma, we have 3 zero eigenvalues
	//    because we set some matrix rows/cols to zero
	bool hasAcousticGamma = false;
	auto s = bandStructure.getState(0);
	auto ens = s.getEnergies();
	if ( ens.size() == crystal.getNumAtoms()*3) hasAcousticGamma = true;
	long firstState;
	if ( hasAcousticGamma ) {
		firstState = 4;
	} else {
		firstState = 1;
	}

	tensordxd.setZero();

	for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {
		auto [imu,it] = loc2Glob(iCalc);
		double c = specificHeat.get(imu,it);
		// is = 3 is a temporary patch, I should discard the first three
		// rows/columns altogether + I skip the bose eigenvector
		for ( long is=firstState; is<relaxationTimes.numStates; is++ ) {

			if ( relaxationTimes.data(iCalc,is) <= 0. ) continue;
//			std::cout << iCalc << " " << is << " " << c << " "
//					<< relaxonV.data.col(is).transpose() << " "
//					<< relaxonV.data.col(is).transpose() << " "
//					<< relaxationTimes.data(iCalc,is) << "\n";

			for ( long i=0; i<dimensionality; i++ ) {
				for ( long j=0; j<dimensionality; j++ ) {
					auto i1 = relaxonV.glob2Loc(imu,it,DimIndex(i));
					auto j1 = relaxonV.glob2Loc(imu,it,DimIndex(j));

					tensordxd(iCalc,i,j) += c
							* relaxonV.data(i1,is)
							* relaxonV.data(j1,is)
							* relaxationTimes.data(iCalc,is);
				}
			}
		}
	}
}

void PhononThermalConductivity::print() {

	std::string units;
	if ( dimensionality == 1 ) {
		units = "W m / K";
	} else if ( dimensionality == 2 ) {
		units = "W / K";
	} else {
		units = "W / m / K";
	}

	std::cout << "\n";
	std::cout << "Thermal Conductivity (" << units << ")\n";

	for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {

		auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
		double temp = calcStat.temperature;

		std::cout << std::fixed;
		std::cout.precision(2);
		std::cout << "Temperature: " << temp * temperatureAuToSi << " (K)\n";
		std::cout.precision(5);
		for ( long i=0; i<dimensionality; i++ ) {
			std::cout << "  " << std::scientific;
			for ( long j=0; j<dimensionality; j++ ) {
				std::cout << " " << std::setw(13) << std::right;
				std::cout << tensordxd(iCalc,i,j)*thConductivityAuToSi;
			}
			std::cout << "\n";
		}
		std::cout << "\n";
	}
}

void PhononThermalConductivity::print(const int & iter) {
	// get the time
	time_t currentTime;
	currentTime = time(NULL);
	// and format the time nicely
	char s[200];
	struct tm * p = localtime(&currentTime);
	strftime(s, 200, "%F, %T", p);


	std::cout << "Iteration: " << iter << " | " << s << "\n";
	for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {
		auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
		double temp = calcStat.temperature;
		std::cout << std::fixed;
		std::cout.precision(2);
		std::cout << "T = " << temp * temperatureAuToSi << ", k = ";
		std::cout.precision(5);
		for ( long i=0; i<dimensionality; i++ ) {
			std::cout << std::scientific;
			std::cout << tensordxd(iCalc,i,i)*thConductivityAuToSi << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}


int PhononThermalConductivity::whichType() {
	return is2Tensor;
}
