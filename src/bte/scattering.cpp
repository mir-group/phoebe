#include "scattering.h"

ScatteringMatrix::ScatteringMatrix(Context & context_,
		StatisticsSweep & statisticsSweep_,
		FullBandStructure & innerBandStructure_,
		FullBandStructure * outerBandStructure_) :
		context(context_), statisticsSweep(statisticsSweep_),
		innerBandStructure(innerBandStructure_),
		outerBandStructure(outerBandStructure_),
		internalDiagonal(context,bandStructure,1) {
	if ( constantRTA ) return;

	if ( outerBandStructure == nullptr ) {
		outerBandStructure = innerBandStructure;
	}

	deltaFunctionSelection = context.getSmearingMethod();

	numStates = outerBandStructure->getNumStates();
	numPoints = outerBandStructure->getNumPoints();

	// the difference arises from the fact that vectors in the BTE have
	// numCalcs set to nTemp * 3, whereas the matrix only depends on nTemp.
	// so, when we compute the action of the scattering matrix, we must take
	// into account for this misalignment

	if ( highMemory ) {
		numCalcs = statisticsSweep.getNumCalcs();
	} else {
		numCalcs = statisticsSweep.getNumCalcs() * context.getDimensionality();
	}

	if ( highMemory ) {

		if ( numCalcs > 1 ) {
			// note: one could write code around this
			// but the methods are very memory intensive for production runs
			Error e("High memory BTE methods can only work with one "
					"temperature and/or chemical potential in a single run");
		}
		theMatrix = Eigen::MatrixXd::Zero(numStates,numStates);
		builder(theMatrix,internalDiagonal,,); // calc matrix and linew.
	} else {
		builder(,internalDiagonal,,); // calc linewidths only
	}
}

ScatteringMatrix(const ScatteringMatrix & that) :  // copy constructor
	context(that.context),
	statisticsSweep(that.statisticsSweep),
	innerBandStructure(that.innerBandStructure),
	outerBandStructure(that.outerBandStructure),
	constantRTA(that.constantRTA),
	highMemory(that.highMemory),
	hasCGScaling(that.hasCGScaling),
	internalDiagonal(that.internalDiagonal),
	theMatrix(that.theMatrix),
	numStates(that.numStates),
	numPoints(that.numPoints),
	numCalcs(that.numCalcs),
	deltaFunctionSelection(that.deltaFunctionSelection) {
}

ScatteringMatrix & operator=(const ScatteringMatrix & that) {//assignment op
	if ( this != &that ) {
		context = that.context;
		statisticsSweep = that.statisticsSweep;
		innerBandStructure = that.innerBandStructure;
		outerBandStructure = that.outerBandStructure;
		constantRTA = that.constantRTA;
		highMemory = that.highMemory;
		hasCGScaling = that.hasCGScaling;
		internalDiagonal = that.internalDiagonal;
		theMatrix = that.theMatrix;
		numStates = that.numStates;
		numPoints = that.numPoints;
		numCalcs = that.numCalcs;
		deltaFunctionSelection = that.deltaFunctionSelection;
	}
	return *this;
}


VectorBTE ScatteringMatrix::diagonal() {
	if ( constantRTA ) {
		double crt = context.getConstantRelaxationTime();
		VectorBTE diag(context,&outerBandStructure,1);
		diag.setConst(1./crt);
		return vector;
	} else {
		return internalDiagonal;
	}
}

VectorBTE ScatteringMatrix::offDiagonalDot(VectorBTE & inPopulation) {
	if ( highMemory ) {
		// it's just the full matrix product, minus the diagonal contribution
		VectorBTE outPopulation = theMatrix * inPopulation;
		outPopulation -= internalDiagonal * inPopulation;
		return outPopulation;
	} else {
		VectorBTE outPopulation(context, &outerBandStructure);
		builder(,,&inPopulation,&outPopulation);
		outPopulation -= internalDiagonal * inPopulation;
		return outPopulation;
	}
}

VectorBTE ScatteringMatrix::dot(VectorBTE & inPopulation) {
	if ( highMemory ) {
		VectorBTE outPopulation = theMatrix * inPopulation;
		if ( hasCGScaling ) {
			outPopulation = outPopulation + inPopulation
					- internalDiagonal * inPopulation;
		}
		return outPopulation;
	} else {
		VectorBTE outPopulation(context,&outerBandStructure);
		builder(,,&inPopulation,&outPopulation);
		if ( hasCGScaling ) {
			outPopulation = outPopulation + inPopulation
					- internalDiagonal * inPopulation;
		}
		return outPopulation;
	}
}

void ScatteringMatrix::setCGScaling() {
	hasCGScaling = true;
}

void ScatteringMatrix::unsetCGScaling() {
	hasCGScaling = false;
}
