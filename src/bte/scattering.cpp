#include "scattering.h"

ScatteringMatrix::ScatteringMatrix(Context & context_,
		StatisticsSweep * statisticsSweep_,
		FullBandStructure<FullPoints> & innerBandStructure_,
		FullBandStructure<FullPoints> & outerBandStructure_) :
		context(context_), statisticsSweep(statisticsSweep_),
		innerBandStructure(innerBandStructure_),
		outerBandStructure(outerBandStructure_),
		internalDiagonal(context,outerBandStructure,1) {

	numStates = outerBandStructure.getNumStates();
	numPoints = outerBandStructure.getNumPoints();

	if ( constantRTA ) return;

	smearing = DeltaFunction::smearingFactory(context,innerBandStructure);

	statisticsSweep = statisticsSweep;

	// the difference arises from the fact that vectors in the BTE have
	// numCalcs set to nTemp * 3, whereas the matrix only depends on nTemp.
	// so, when we compute the action of the scattering matrix, we must take
	// into account for this misalignment
	if ( highMemory ) {
		numCalcs = statisticsSweep->getNumCalcs();
	} else {
		numCalcs = statisticsSweep->getNumCalcs() *context.getDimensionality();
	}
}

// copy constructor
ScatteringMatrix::ScatteringMatrix(const ScatteringMatrix & that) :
	context(that.context),
	statisticsSweep(that.statisticsSweep),
	smearing(that.smearing),
	innerBandStructure(that.innerBandStructure),
	outerBandStructure(that.outerBandStructure),
	constantRTA(that.constantRTA),
	highMemory(that.highMemory),
	hasCGScaling(that.hasCGScaling),
	internalDiagonal(that.internalDiagonal),
	theMatrix(that.theMatrix),
	numStates(that.numStates),
	numPoints(that.numPoints),
	numCalcs(that.numCalcs) {
}

//assignment operator
ScatteringMatrix & ScatteringMatrix::operator=(const ScatteringMatrix & that) {
	if ( this != &that ) {
		context = that.context;
		statisticsSweep = that.statisticsSweep;
		smearing = that.smearing;
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
	}
	return *this;
}


VectorBTE ScatteringMatrix::diagonal() {
	if ( constantRTA ) {
		double crt = context.getConstantRelaxationTime();
		VectorBTE diag(context,outerBandStructure,1);
		diag.setConst(1./crt);
		return diag;
	} else {
		return internalDiagonal;
	}
}

VectorBTE ScatteringMatrix::offDiagonalDot(VectorBTE & inPopulation) {
	if ( highMemory ) {
		// it's just the full matrix product, minus the diagonal contribution
		VectorBTE outPopulation = dot(inPopulation);
//		outPopulation -= internalDiagonal * inPopulation;
		// must resolve the different cartesian indices
		for ( long i=0; i<outPopulation.numCalcs; i++ ) {
			outPopulation.data(0,i) -=
				- internalDiagonal.data(0,i/internalDiagonal.dimensionality)
				* inPopulation.data(0,i);
		}
		return outPopulation;
	} else {
		VectorBTE outPopulation(context, outerBandStructure);
		builder(nullptr,nullptr,&inPopulation,&outPopulation);
		// outPopulation = outPopulation - internalDiagonal * inPopulation;
		for ( long i=0; i<outPopulation.numCalcs; i++ ) {
			outPopulation.data(0,i) -=
				- internalDiagonal.data(0,i/internalDiagonal.dimensionality)
				* inPopulation.data(0,i);
		}
		return outPopulation;
	}
}

VectorBTE ScatteringMatrix::dot(VectorBTE & inPopulation) {
	if ( highMemory ) {
		VectorBTE outPopulation(context, outerBandStructure);
		// note: we are assuming that ScatteringMatrix has numCalcs = 1
		for ( auto i=0; i<numStates; i++ ) {
			for ( auto j=0; j<numStates; j++ ) {
				outPopulation.data(0,i) +=
						theMatrix(i,j) * inPopulation.data(0,j);
			}
		}
		if ( hasCGScaling ) {
			// notice that the linewidths are "scalar",
			// populations instead have cartesian indeces
			for ( long i=0; i<outPopulation.numCalcs; i++ ) {
				outPopulation.data(0,i) += inPopulation.data(0,i)
					-internalDiagonal.data(0,i/internalDiagonal.dimensionality)
					* inPopulation.data(0,i);
			}
		}
		return outPopulation;
	} else {
		VectorBTE outPopulation(context, outerBandStructure);
		builder(nullptr,nullptr,&inPopulation,&outPopulation);
		if ( hasCGScaling ) {
			for ( long i=0; i<outPopulation.numCalcs; i++ ) {
				outPopulation.data(0,i) += inPopulation.data(0,i)
					-internalDiagonal.data(0,i/internalDiagonal.dimensionality)
					* inPopulation.data(0,i);
			}
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
