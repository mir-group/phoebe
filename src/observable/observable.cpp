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
