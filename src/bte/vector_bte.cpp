#include "vector_bte.h"

// default constructor
VectorBTE::VectorBTE(Context & context_,
		ActiveBandStructure & activeBandStructure_) :
		context(context_), activeBandStructure(activeBandStructure_) {
	Eigen::VectorXd temperatures = context.getTemperatures();
	Eigen::VectorXd chemicalPotentials = context.getChemicalPotentials();
	numChemPots = chemicalPotentials.size();
	numTemps = temperatures.size();
	dimensionality = context.getDimensionality();
	numCalcs = numTemps * numChemPots * dimensionality;
	numStates = activeBandStructure.getNumStates();
	data = Eigen::MatrixXd::Zero(numCalcs,numStates);
}

// copy constructor
VectorBTE::VectorBTE(const VectorBTE & that)
		: VectorBTE(that.context, that.activeBandStructure) {
	data = that.data;
}

// copy assigmnent
VectorBTE & VectorBTE::operator = (const VectorBTE & that) {
    if ( this != &that) {
    	activeBandStructure = that.activeBandStructure;
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
	VectorBTE newPopulation(context, activeBandStructure);
	for ( long i=0; i<numCalcs; i++ ) {
		newPopulation.data.row(i) = this->data.row(i) * scalar;
	}
	return newPopulation;
}

// product operator overload
VectorBTE VectorBTE::operator * (const Eigen::VectorXd & vector) {
	VectorBTE newPopulation(context, activeBandStructure);
	for ( long i=0; i<numCalcs; i++ ) {
		newPopulation.data.row(i) = this->data.row(i) * vector(i);
	}
	return newPopulation;
}

// product operator overload
VectorBTE VectorBTE::operator + (const VectorBTE & that) {
	VectorBTE newPopulation(context, activeBandStructure);
	newPopulation.data = this->data + that.data;
	return newPopulation;
}

// product operator overload
VectorBTE VectorBTE::operator - (const VectorBTE & that) {
	VectorBTE newPopulation(context, activeBandStructure);
	newPopulation.data = this->data - that.data;
	return newPopulation;
}

// product operator overload
VectorBTE VectorBTE::operator - () {
	VectorBTE newPopulation(context, activeBandStructure);
	newPopulation.data = - this->data;
	return newPopulation;
}

// product operator overload
VectorBTE VectorBTE::operator / (const VectorBTE & that) {
	VectorBTE newPopulation(context, activeBandStructure);
	newPopulation.data << this->data.array() / that.data.array();
	return newPopulation;
}

VectorBTE VectorBTE::sqrt() {
	VectorBTE newPopulation(context, activeBandStructure);
	newPopulation.data << this->data.array().sqrt();
	return newPopulation;
}

long VectorBTE::glob2Loc(long & imu, long & it, long & idim) {
	long i = imu * numTemps * dimensionality + it * dimensionality + idim;
	return i;
}

std::tuple<long,long,long> VectorBTE::loc2Glob(long & i) {
	long imu = i / (numTemps * dimensionality);
	long remainder = i - imu * numTemps * dimensionality;
	long it = remainder / dimensionality;
	remainder -= it * dimensionality;
	long idim = remainder;
	return {imu,it,idim};
}
