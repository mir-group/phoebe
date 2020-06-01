#include "vector_bte.h"

// default constructor
VectorBTE::VectorBTE(Context & context_,
		FullBandStructure<FullPoints> & bandStructure_,
		const long & dimensionality_) :
		context(context_), bandStructure(bandStructure_) {

	Eigen::VectorXd temperatures = context.getTemperatures();
	Eigen::VectorXd chemicalPotentials = context.getChemicalPotentials();
	numChemPots = chemicalPotentials.size();
	numTemps = temperatures.size();

	if ( dimensionality_ == 0 ) {
		dimensionality = context.getDimensionality();
	} else if ( dimensionality_ < 0 ) {
		Error e("VectorBTE doesn't accept negative dimensions");
		dimensionality = 0;
	} else {
		dimensionality = dimensionality_;
	}

	numCalcs = numTemps * dimensionality;
	if ( numChemPots > 0 ) numCalcs *= numChemPots;
	numStates = bandStructure.getNumStates();
	data = Eigen::MatrixXd::Zero(numCalcs,numStates);
}

// copy constructor
VectorBTE::VectorBTE(const VectorBTE & that)
		: context(that.context), bandStructure(that.bandStructure) {
	numChemPots = that.numChemPots;
	numTemps = that.numTemps;
	numCalcs = that.numCalcs;
	numStates = that.numStates;
	dimensionality = that.dimensionality;
	data = that.data;
}

// copy assignment
VectorBTE & VectorBTE::operator = (const VectorBTE & that) {
    if ( this != &that) {
    	bandStructure = that.bandStructure;
		context = that.context;
    	numChemPots = that.numChemPots;
		numTemps = that.numTemps;
		dimensionality = that.dimensionality;
		numCalcs = that.numCalcs;
		numStates = that.numStates;
		data = that.data;
    }
    return *this;
}

// product operator overload
Eigen::VectorXd VectorBTE::operator * (const VectorBTE & that) {
	Eigen::VectorXd result(numCalcs);
	for ( long i=0; i<numCalcs; i++ ) {
		result.row(i) = this->data.row(i).transpose() * that.data.row(i);
	}
	return result;
}

// product operator overload
VectorBTE VectorBTE::operator * (const double & scalar) {
	VectorBTE newPopulation(context, bandStructure);
	for ( long i=0; i<numCalcs; i++ ) {
		newPopulation.data.row(i) = this->data.row(i) * scalar;
	}
	return newPopulation;
}

// product operator overload
VectorBTE VectorBTE::operator * (const Eigen::VectorXd & vector) {
	VectorBTE newPopulation(context, bandStructure);
	for ( long i=0; i<numCalcs; i++ ) {
		newPopulation.data.row(i) = this->data.row(i) * vector(i);
	}
	return newPopulation;
}

// product operator overload
VectorBTE VectorBTE::operator + (const VectorBTE & that) {
	VectorBTE newPopulation(context, bandStructure);
	newPopulation.data = this->data + that.data;
	return newPopulation;
}

// product operator overload
VectorBTE VectorBTE::operator - (const VectorBTE & that) {
	VectorBTE newPopulation(context, bandStructure);
	newPopulation.data = this->data - that.data;
	return newPopulation;
}

// difference operator overload
VectorBTE VectorBTE::operator - () {
	VectorBTE newPopulation(context, bandStructure);
	newPopulation.data = - this->data;
	return newPopulation;
}

// division operator overload
VectorBTE VectorBTE::operator / (const VectorBTE & that) {
	VectorBTE newPopulation(context, bandStructure);
	//
	if ( dimensionality == that.dimensionality ) {
		newPopulation.data << this->data.array() / that.data.array();
	} else if ( that.dimensionality == 1 ) {
		long i1, i2;
		for ( long imu=0; imu<numChemPots; imu++ ) {
			for ( long it=0; it<numTemps; it++ ) {
				i2 = glob2Loc(imu,it,0);
				for ( long idim=0; idim<dimensionality; idim++ ) {
					i1 = glob2Loc(imu,it,idim);
					newPopulation.data.row(i1) = this->data.row(i1).array()
							/ that.data.row(i2).array();
				}
			}
		}
	}
	return newPopulation;
}

VectorBTE VectorBTE::sqrt() {
	VectorBTE newPopulation(context, bandStructure);
	newPopulation.data << this->data.array().sqrt();
	return newPopulation;
}

void VectorBTE::setConst(const double & constant) {
	data.setConstant(constant);
}

long VectorBTE::glob2Loc(const long & imu, const long & it, const long & idim){
	long i = compress3Indeces(imu,it,idim,numChemPots,numTemps,dimensionality);
	return i;
}

std::tuple<long,long,long> VectorBTE::loc2Glob(const long & i) {
	auto [imu, it, idim] = decompress3Indeces(i, numChemPots, numTemps,
			dimensionality);
	return {imu,it,idim};
}
