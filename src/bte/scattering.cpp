#include "scattering.h"

ScatteringMatrix::ScatteringMatrix(Context & context_,
		ActiveBandStructure & bandStructure_,
		StatisticsSweep & statisticsSweep_) :
		context(context_), bandStructure(bandStructure_),
		internalDiagonal(context,bandStructure),
		statisticsSweep(statisticsSweep_) {
	if ( constantRTA ) return;

	deltaFunctionSelection = context.getSmearingMethod();

	numStates = bandStructure.getNumStates();
	numPoints = bandStructure.getNumPoints();

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
		internalDiagonal = builderManager();
	} else {
		setInternalDiagonal(); // these we can store also with low mem
	}
}

ScatteringMatrix(const ScatteringMatrix & that) :  // copy constructor
	context(that.context),
	bandStructure(that.bandStructure),
	statisticsSweep(that.statisticsSweep),
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
		bandStructure = that.bandStructure;
		statisticsSweep = that.statisticsSweep;
		constantRTA = that.constantRTA;
		highMemory = that.highMemory;
		hasCGScaling = thaht.hasCGScaling;
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
		VectorBTE diag(context,bandStructure);
		diag.setConst(1./crt);
		return vector;
	} else {
		if ( ! hasCGScaling ) {
			return internalDiagonal;
		} else {
			VectorBTE outVector(context,bandStructure);
			outVector.setConst(1.);
			return outVector;
		}
	}
}

void ScatteringMatrix::setInternalDiagonal() {
	internalDiagonal = builderManager();
}

VectorBTE ScatteringMatrix::offDiagonalDot(VectorBTE & inPopulation) {
	if ( highMemory ) {
		// it's just the full matrix product, minus the diagonal contribution
		VectorBTE outPopulation = theMatrix * inPopulation;
		outPopulation -= internalDiagonal * inPopulation;
		return outPopulation;
	} else {
		VectorBTE outPopulation = builderManager(&inPopulation);
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
		VectorBTE outPopulation = builderManager(&inPopulation);
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









VectorBTE PhScatteringMatrix::builderManager(VectorBTE * inPopulation) {
	VectorBTE outVector(context, bandStructure);
	outVector.setConst(0.);
	if ( inPopulation == nullptr ) {
		// case in which we want to compute the diagonal only
		if ( couplingPh3 != nullptr ) outVector += internal3Ph(,);
	} else {
		if ( highMemory ) {
			if ( couplingPh3 != nullptr ) {
				outVector += builder3Ph(*theMatrix,);
			}
		} else { // compute on the fly
			if ( couplingPh3 != nullptr ) {
				outVector += builder3Ph(,*inPopulation);
			}
		}
	}
	return outVector;
}


//VectorBTE builderPhIsotope(Eigen::MatrixXd * theMatrix=nullptr,
//		VectorBTE * inPopulation=nullptr) {
//}
//
//VectorBTE builderPhBoundary(Eigen::MatrixXd * theMatrix=nullptr,
//		VectorBTE * inPopulation=nullptr) {
//}

// 3 cases:
// theMatrix is passed: we compute and store in memory the scatt matrix
//                      we return the diagonal
// inPopulation is passed: compute sMatrix * vector, matrix not kept in memory
//                      we return outVec = sMatrix*vector
// neither is passed: we compute and return the diagonal of the scatt matrix
VectorBTE PhScatteringMatrix::builder3Ph(Eigen::MatrixXd * matrix,
		VectorBTE * inPopulation) {

	VectorBTE outVector(context, bandStructure);

	for ( i ) {
		for ( j ) {
			if ( inPopulation != nullptr ) {
				outVector(i) = term * inPopulation(j);
			if ( matrix != nullptr ) {
				*matrix(i,j) += term;
				outVector(i) += term; // cumulate the diagonal
			} else {
				outVector(i) += term; // cumulate the diagonal
			}
		}
	}

	return outVector;
}
