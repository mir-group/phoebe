#include "phonon_thermal_cond.h"
#include "constants.h"

PhononThermalConductivity::PhononThermalConductivity(
		StatisticsSweep & statisticsSweep_,
		Crystal & crystal_, FullBandStructure<FullPoints> & bandStructure_) :
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
			long is = bandStructure.getIndex(ik,ib);

			// skip the acoustic phonons
			if ( s.getCoords(Points::crystalCoords).norm()==0. && ib<3 ) continue;

			for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {
				auto [imu,it] = loc2Glob(iCalc);

				for ( long i=0; i<dimensionality; i++ ) {
					long icPop = n.glob2Loc(imu,it,i);
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
			long icPop1 = f.glob2Loc(imu,it,i);
			for ( long j=0; j<dimensionality; j++ ) {
				long icPop2 = f.glob2Loc(imu,it,j);
				for ( long is=0; is<bandStructure.getNumStates(); is++ ) {
					tensordxd(iCalc,i,j) -=
							f.data(icPop1,is) * af.data(icPop2,is) * norm;
				}
			}
		}
	}
}


//void PhononThermalConductivity::calcVariational(VectorBTE & af,
//		VectorBTE & f, VectorBTE & b) {
//	Eigen::VectorXd temperatures = context.getTemperatures();
//	Eigen::VectorXd lambda = 1. / f.bandStructure.getNumPoints()
//			/ crystal.getVolumeUnitCell(dimensionality) / temperatures.array();
//	long numStates = f.numStates;
//
//	for ( long it=0; it<numTemps; it++ ) {
//		for ( long imu=0; imu<numChemPots; imu++ ) {
//			long icLoc = glob2Loc(imu, it);
//			for ( long i=0; i<dimensionality; i++ ) {
//				for ( long j=0; j<dimensionality; j++ ) {
//					long icPop1 = f.glob2Loc(imu,it,i);
//					long icPop2 = b.glob2Loc(imu,it,j);
//					for ( long istate=0; istate<numStates; istate++ ) {
//						tensordxd(icLoc,i,j) += 0.5 * f.data(icPop1,istate)
//								* af.data(icPop2,istate);
//						tensordxd(icLoc,i,j) -= f.data(icPop1,istate)
//								* b.data(icPop2,istate);
//					}
//					tensordxd(icLoc,i,j) /= lambda(it);
//				}
//			}
//		}
//	}
//}

void PhononThermalConductivity::print() {
	for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {

		auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
		double temp = calcStat.temperature;

		std::cout << "\n";
		std::cout << "Thermal Conductivity (W/m/K)\n";
		std::cout.precision(5);
		std::cout << "Temperature: " << temp * temperatureAuToSi << " (K)\n";
		std::cout << std::scientific;
		for ( long i=0; i<dimensionality; i++ ) {
			for ( long j=0; j<dimensionality; j++ ) {
				std::cout << tensordxd(iCalc,i,j)*thConductivityAuToSi << " ";
			}
			std::cout << "\n";
		}
		std::cout << "\n";
	}
}

void PhononThermalConductivity::print(const int & iter) {
	std::cout << "Iteration: " << iter << "\n";
	for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {
		auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
		double temp = calcStat.temperature;
		std::cout.precision(5);
		std::cout << "T = " << temp * temperatureAuToSi << ", k = ";
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
