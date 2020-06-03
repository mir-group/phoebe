#include "phonon_thermal_cond.h"
#include "constants.h"

PhononThermalConductivity::PhononThermalConductivity(Context & context_,
		Crystal & crystal_, FullBandStructure<FullPoints> & bandStructure_) :
				Observable(context_, crystal_),
				bandStructure(bandStructure_) {
};

void PhononThermalConductivity::calcFromPopulation(VectorBTE & n) {
	double norm = 1. / bandStructure.getNumPoints()
			/ crystal.getVolumeUnitCell(dimensionality);

	tensordxd = Eigen::Tensor<double,3>(numCalcs,dimensionality,dimensionality);
	tensordxd.setZero();

	long is = -1;
	for ( long ik=0; ik<bandStructure.getNumPoints(); ik++ ) {
		auto s = bandStructure.getState(ik);
		auto en = s.getEnergies();
		auto vel = s.getVelocities();
		for ( long ib=0; ib<en.size(); ib++ ) {

			is++;
			// skip the acoustic phonons
			if ( ik==0 && ib<3 ) continue;

			for ( long it=0; it<numTemps; it++ ) {
				for ( long imu=0; imu<numChemPots; imu++ ) {
					long icLoc = glob2Loc(imu, it);
					for ( long i=0; i<dimensionality; i++ ) {
						long icPop = n.glob2Loc(imu,it,i);

						for ( long j=0; j<dimensionality; j++ ) {
							tensordxd(icLoc,i,j) +=
									n.data(icPop,is) * vel(ib,ib,j).real()
									* en(ib) * norm;
						}
					}
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
	auto temperatures = context.getTemperatures();
	for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {
		auto [imu,it] = loc2Glob(iCalc);

		std::cout << "\n";
		std::cout << "Thermal Conductivity (W/m/K)\n";
		std::cout.precision(5);
		std::cout << "Temperature: " << temperatures(it)*temperatureAuToSi
				<< " (K)\n";
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
