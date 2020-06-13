#include "vector_bte.h"
#include "constants.h"

// default constructor
VectorBTE::VectorBTE(StatisticsSweep & statisticsSweep_,
		FullBandStructure & bandStructure_,
		const long & dimensionality_) :
		statisticsSweep(statisticsSweep_), bandStructure(bandStructure_) {

	if ( dimensionality_ < 0 ) {
		Error e("VectorBTE doesn't accept negative dimensions");
		dimensionality = 0;
	} else {
		dimensionality = dimensionality_;
	}

	numCalcs = statisticsSweep.getNumCalcs();
	numCalcs *= dimensionality;

	numChemPots = statisticsSweep.getNumChemicalPotentials();
	numTemps = statisticsSweep.getNumTemperatures();

	numStates = bandStructure.getNumStates();
	data = Eigen::MatrixXd::Zero(numCalcs,numStates);

	for ( long is=0; is<numStates; is++ ) {
		double en = bandStructure.getEnergy(is);
		if ( en < 0.1 / ryToCmm1 ) { // cutoff at 0.1 cm^-1
			excludeIndeces.push_back(is);
		}
	}
}

// copy constructor
VectorBTE::VectorBTE(const VectorBTE & that)
		: statisticsSweep(that.statisticsSweep),
		  bandStructure(that.bandStructure) {
	numCalcs = that.numCalcs;
	numStates = that.numStates;
	numChemPots = that.numChemPots;
	numTemps = that.numTemps;
	dimensionality = that.dimensionality;
	data = that.data;
	excludeIndeces = that.excludeIndeces;
}

// copy assignment
VectorBTE & VectorBTE::operator = (const VectorBTE & that) {
    if ( this != &that) {
    	statisticsSweep = that.statisticsSweep;
    	bandStructure = that.bandStructure;
		numCalcs = that.numCalcs;
		numStates = that.numStates;
		numChemPots = that.numChemPots;
		numTemps = that.numTemps;
		dimensionality = that.dimensionality;
		data = that.data;
		excludeIndeces = that.excludeIndeces;
    }
    return *this;
}

// product operator overload
Eigen::VectorXd VectorBTE::dot(const VectorBTE & that) {
	Eigen::VectorXd result(numCalcs);
	result.setZero();
	for ( long i=0; i<numCalcs; i++ ) {
		for ( long is=0; is<numStates; is++ ) {
			result(i) += this->data(i,is) * that.data(i,is);
		}
	}
	return result;
}


VectorBTE VectorBTE::baseOperator(VectorBTE & that,
		const int & operatorType) {
	VectorBTE newPopulation(statisticsSweep, bandStructure, dimensionality);

	if ( dimensionality == that.dimensionality ) {

		if ( operatorType == operatorSums ) {
			newPopulation.data << this->data.array() + that.data.array();
		} else if ( operatorType == operatorDivs ) {
			newPopulation.data << this->data.array() / that.data.array();
		} else if ( operatorType == operatorProd ) {
			newPopulation.data << this->data.array() * that.data.array();
		} else if ( operatorType == operatorDiff ) {
			newPopulation.data << this->data.array() - that.data.array();
		} else {
			Error e("Operator type for VectorBTE not recognized");
		}

	} else if ( that.dimensionality == 1 ) {

		for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {
			auto [imu,it,idim] = loc2Glob(iCalc);
			auto i2 = that.glob2Loc(imu,it,DimIndex(0));

			if ( operatorType == operatorSums ) {
				newPopulation.data.row(iCalc) = this->data.row(iCalc).array()
						+ that.data.row(i2).array();
			} else if ( operatorType == operatorDivs ) {
				newPopulation.data.row(iCalc) = this->data.row(iCalc).array()
						/ that.data.row(i2).array();
			} else if ( operatorType == operatorProd ) {
				newPopulation.data.row(iCalc) = this->data.row(iCalc).array()
						* that.data.row(i2).array();
			} else if ( operatorType == operatorDiff ) {
				newPopulation.data.row(iCalc) = this->data.row(iCalc).array()
						- that.data.row(i2).array();
			} else {
				Error e("Operator type for VectorBTE not recognized");
			}
		}
	} else {
		Error e("VectorBTE can't handle dimensionality for this case");
	}
	for ( auto is : excludeIndeces ) {
		newPopulation.data.col(is).setZero();
	}
	return newPopulation;
}

// product operator overload
VectorBTE VectorBTE::operator * (VectorBTE & that) {
	return baseOperator(that, operatorProd);
}

// product operator overload
VectorBTE VectorBTE::operator * (const double & scalar) {
	VectorBTE newPopulation(statisticsSweep, bandStructure, dimensionality);
	for ( long i=0; i<numCalcs; i++ ) {
		newPopulation.data.row(i) = this->data.row(i) * scalar;
	}
	return newPopulation;
}

// product operator overload
VectorBTE VectorBTE::operator * (const Eigen::VectorXd & vector) {
	VectorBTE newPopulation(statisticsSweep, bandStructure, dimensionality);
	for ( long i=0; i<numCalcs; i++ ) {
		newPopulation.data.row(i) = this->data.row(i) * vector(i);
	}
	return newPopulation;
}

// product operator overload
VectorBTE VectorBTE::operator + (VectorBTE & that) {
	return baseOperator(that, operatorSums);
}

// product operator overload
VectorBTE VectorBTE::operator - (VectorBTE & that) {
	return baseOperator(that, operatorDiff);
}

// difference operator overload
VectorBTE VectorBTE::operator - () {
	VectorBTE newPopulation(statisticsSweep, bandStructure, dimensionality);
	newPopulation.data = - this->data;
	return newPopulation;
}

// division operator overload
VectorBTE VectorBTE::operator / (VectorBTE & that) {
	return baseOperator(that, operatorDivs);
}

VectorBTE VectorBTE::sqrt() {
	VectorBTE newPopulation(statisticsSweep, bandStructure, dimensionality);
	newPopulation.data << this->data.array().sqrt();
	return newPopulation;
}

VectorBTE VectorBTE::reciprocal() {
	VectorBTE newPopulation(statisticsSweep, bandStructure, dimensionality);
	newPopulation.data << 1. / this->data.array();
	return newPopulation;
}

void VectorBTE::setConst(const double & constant) {
	data.setConstant(constant);
}

long VectorBTE::glob2Loc(const ChemPotIndex & imu, const TempIndex & it,
		const DimIndex & idim){
	long i = compress3Indeces(imu.get(),it.get(),idim.get(),
			numChemPots,numTemps,dimensionality);
	return i;
}

std::tuple<ChemPotIndex,TempIndex,DimIndex> VectorBTE::loc2Glob(
		const long & i) {
	auto [imu, it, idim] = decompress3Indeces(i, numChemPots, numTemps,
			dimensionality);
	return {ChemPotIndex(imu),TempIndex(it),DimIndex(idim)};
}

void VectorBTE::canonical2Population() {
	auto particle = bandStructure.getParticle();
	for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {
		auto [imu, it, idim] = loc2Glob(iCalc);
		auto calcStatistics = statisticsSweep.getCalcStatistics(it, imu);
		auto temp = calcStatistics.temperature;
		auto chemPot = calcStatistics.chemicalPotential;
		for ( long is=0; is<numStates; is++ ) {
			double en = bandStructure.getEnergy(is);
			double term = particle.getPopPopPm1(en, temp, chemPot);
			data(iCalc,is) *= term;
		}
	}
}

void VectorBTE::population2Canonical() {
	auto particle = bandStructure.getParticle();
	for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {
		auto [imu, it, idim] = loc2Glob(iCalc);
		auto calcStatistics = statisticsSweep.getCalcStatistics(it, imu);
		auto temp = calcStatistics.temperature;
		auto chemPot = calcStatistics.chemicalPotential;
		for ( long is=0; is<numStates; is++ ) {
			double en = bandStructure.getEnergy(is);
			double term = particle.getPopPopPm1(en, temp, chemPot);
			data(iCalc,is) /= term;
		}
	}
}
