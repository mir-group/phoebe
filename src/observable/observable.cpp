#include "observable.h"
#include <cmath>
#include "constants.h"

Observable::Observable(Context & context_, Crystal & crystal_) :
		context(context_), crystal(crystal_) {
	Eigen::VectorXd temperatures = context.getTemperatures();
	Eigen::VectorXd chemicalPotentials = context.getChemicalPotentials();
	numChemPots = chemicalPotentials.size();
	numTemps = temperatures.size();
	numCalcs = numTemps * numChemPots;
	dimensionality = context.getDimensionality();
	//note: each subclass should define type
}

// copy constructor
Observable::Observable(const Observable & that) :
		context(that.context), crystal(that.crystal) {
	numChemPots = that.numChemPots;
	numTemps = that.numTemps;
	dimensionality = that.dimensionality;
	numCalcs = that.numCalcs;
	type = that.type;
	scalar = that.scalar;
	vectord = that.vectord;
	tensordxd = that.tensordxd;
	tensordxdxdxd = that.tensordxdxdxd;
}

// copy assigmnent
Observable & Observable::operator = (const Observable & that) {
    if ( this != &that) {
    	context = that.context;
    	crystal = that.crystal;
		numChemPots = that.numChemPots;
		numTemps = that.numTemps;
		numCalcs = that.numCalcs;
		type = that.type;
		scalar = that.scalar;
		vectord = that.vectord;
		tensordxd = that.tensordxd;
		tensordxdxdxd = that.tensordxdxdxd;
    }
    return *this;
}

//void Observable::add(long & it, long & imu, double & term) {
//	long is = glob2loc(imu,it);
//	scalar(is) += term;
//}
//
//void Observable::add(long & it, long & imu, double & term, long & i) {
//	long is = glob2loc(imu,it);
//	vectord(is,i) += term;
//}
//
//void Observable::add(long & it, long & imu, double & term, long & i, long & j) {
//	long is = glob2loc(imu,it);
//	tensordxd(is,i,j) += term;
//}
//
//void Observable::add(long & it, long & imu, double & term,
//		int & i, int & j, int & k, int & l) {
//	long is = glob2loc(imu,it);
//	tensordxdxdxd(is,i,j,k,l) += term;
//}
//
//void Observable::rescale(long & it, long & imu, double & term) {
//	long is = glob2loc(imu,it);
//	if ( type == isScalar ) {
//		vectord(is) *= term;
//	} else if ( type == isVector ) {
//	 	for ( int i=0; i<dimensionality; i++ ) {
//	 		vectord(is,i) *= term;
//		}
//	} else if ( type == is2Tensor ) {
//	 	for ( int i=0; i<dimensionality; i++ ) {
//		 	for ( int j=0; j<dimensionality; j++ ) {
//				 		tensordxd(is,i,j) *= term;
//		 	}
//	 	}
//	} else if ( type == is4Tensor ) {
//	 	for ( int i=0; i<dimensionality; i++ ) {
//		 	for ( int j=0; j<dimensionality; j++ ) {
//			 	for ( int k=0; k<dimensionality; k++ ) {
//				 	for ( int l=0; l<dimensionality; l++ ) {
//				 		tensordxdxdxd(is,i,j,k,l) *= term;
//				 	}
//			 	}
//		 	}
//	 	}
//	} else {
//		Error e("type not defined in observable",1);
//	}
//}

long Observable::glob2Loc(const long & imu, const long & it) {
	return compress2Indeces(imu, it, numChemPots, numTemps);
}

std::tuple<long,long> Observable::loc2Glob(const long & i) {
	return decompress2Indeces(i, numChemPots, numTemps);
}

Observable Observable::operator - (const Observable & that) {
	Observable newObservable(context, crystal);
	if ( type == isScalar ) {
		for ( long is=0; is<numCalcs; is++ ) {
			newObservable.scalar(is) = scalar(is) - that.scalar(is);
		}
	} else if ( type == isVector ) {
		for ( long is=0; is<numCalcs; is++ ) {
			for ( int i=0; i<dimensionality; i++ ) {
				newObservable.vectord(is,i) = vectord(is,i) - that.vectord(is,i);
			}
		}
	} else if ( type == is2Tensor ) {
		for ( long is=0; is<numCalcs; is++ ) {
			for ( int i=0; i<dimensionality; i++ ) {
				for ( int j=0; j<dimensionality; j++ ) {
					newObservable.tensordxd(is,i,j) = tensordxd(is,i,j)
		 						- that.tensordxd(is,i,j);
				}
		 	}
	 	}
	} else if ( type == is4Tensor ) {
		for ( long is=0; is<numCalcs; is++ ) {
			for ( int i=0; i<dimensionality; i++ ) {
				for ( int j=0; j<dimensionality; j++ ) {
					for ( int k=0; k<dimensionality; k++ ) {
						for ( int l=0; l<dimensionality; l++ ) {
							newObservable.tensordxdxdxd(is,i,j,k,l) =
									tensordxdxdxd(is,i,j,k,l)
									- that.tensordxdxdxd(is,i,j,k,l);
						}
				 	}
			 	}
		 	}
	 	}
	}
	return newObservable;
}

Eigen::VectorXd Observable::getNorm() {
	Eigen::VectorXd norm(numCalcs);
	norm.setZero();
	if ( type == isScalar ) {
		for ( long is=0; is<numCalcs; is++ ) {
			norm(is) = abs(scalar(is));
		}
	} else if ( type == isVector ) {
		for ( long is=0; is<numCalcs; is++ ) {
			for ( int i=0; i<dimensionality; i++ ) {
				norm(is) += vectord(is,i) * vectord(is,i);
			}
			norm(is) = sqrt(norm(is)) / double(dimensionality);
		}
	} else if ( type == is2Tensor ) {
		for ( long is=0; is<numCalcs; is++ ) {
			for ( int i=0; i<dimensionality; i++ ) {
				for ( int j=0; j<dimensionality; j++ ) {
					norm(is) += tensordxd(is,i,j) * tensordxd(is,i,j);
				}
		 	}
			norm(is) = sqrt(norm(is)) / double(dimensionality*dimensionality);
	 	}
	} else if ( type == is4Tensor ) {
		for ( long is=0; is<numCalcs; is++ ) {
			for ( int i=0; i<dimensionality; i++ ) {
				for ( int j=0; j<dimensionality; j++ ) {
					for ( int k=0; k<dimensionality; k++ ) {
						for ( int l=0; l<dimensionality; l++ ) {
							norm(is) += tensordxdxdxd(is,i,j,k,l)
									* tensordxdxdxd(is,i,j,k,l);
						}
				 	}
			 	}
		 	}
			norm(is) = sqrt(norm(is)) / double(pow(dimensionality,4));
	 	}
	}
	return norm;
}

PhononThermalConductivity::PhononThermalConductivity(Context & context_,
		Crystal & crystal_) : Observable(context_, crystal_) {
};

void PhononThermalConductivity::calcFromPopulation(VectorBTE & f,
		VectorBTE & b) {
	Eigen::VectorXd temperatures = context.getTemperatures();
	Eigen::VectorXd lambda(numTemps);
	lambda << 1. / f.bandStructure.getNumPoints()
			/ crystal.getVolumeUnitCell(dimensionality) / temperatures.array();
	long numStates = f.numStates;

	tensordxd = Eigen::Tensor<double,3>(numCalcs,dimensionality,dimensionality);
	tensordxd.setZero();

	for ( long it=0; it<numTemps; it++ ) {
		for ( long imu=0; imu<numChemPots; imu++ ) {
			long icLoc = glob2Loc(imu, it);
			for ( long i=0; i<dimensionality; i++ ) {
				for ( long j=0; j<dimensionality; j++ ) {
					long icPop1 = f.glob2Loc(imu,it,i);
					long icPop2 = b.glob2Loc(imu,it,j);
					//TODO: I'm skipping the acoustic phonons
					// we should clean this offset
					for ( long istate=3; istate<numStates; istate++ ) {
						tensordxd(icLoc,i,j) += f.data(icPop1,istate)
								* b.data(icPop2,istate);
					}
					tensordxd(icLoc,i,j) /= lambda(it);
				}
			}
		}
	}
}

void PhononThermalConductivity::calcVariational(VectorBTE & af,
		VectorBTE & f, VectorBTE & b) {
	Eigen::VectorXd temperatures = context.getTemperatures();
	Eigen::VectorXd lambda = 1. / f.bandStructure.getNumPoints()
			/ crystal.getVolumeUnitCell(dimensionality) / temperatures.array();
	long numStates = f.numStates;

	for ( long it=0; it<numTemps; it++ ) {
		for ( long imu=0; imu<numChemPots; imu++ ) {
			long icLoc = glob2Loc(imu, it);
			for ( long i=0; i<dimensionality; i++ ) {
				for ( long j=0; j<dimensionality; j++ ) {
					long icPop1 = f.glob2Loc(imu,it,i);
					long icPop2 = b.glob2Loc(imu,it,j);
					for ( long istate=0; istate<numStates; istate++ ) {
						tensordxd(icLoc,i,j) += 0.5 * f.data(icPop1,istate)
								* af.data(icPop2,istate);
						tensordxd(icLoc,i,j) -= f.data(icPop1,istate)
								* b.data(icPop2,istate);
					}
					tensordxd(icLoc,i,j) /= lambda(it);
				}
			}
		}
	}
}

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
