#include "scattering.h"
#include "constants.h"
#include <algorithm>

ScatteringMatrix::ScatteringMatrix(Context & context_,
		StatisticsSweep & statisticsSweep_,
		FullBandStructure<FullPoints> & innerBandStructure_,
		FullBandStructure<FullPoints> & outerBandStructure_) :
		context(context_), statisticsSweep(statisticsSweep_),
		innerBandStructure(innerBandStructure_),
		outerBandStructure(outerBandStructure_),
		internalDiagonal(statisticsSweep, outerBandStructure, 1) {

	numStates = outerBandStructure.getNumStates();
	numPoints = outerBandStructure.getNumPoints();

	double constantRelaxationTime = context.getConstantRelaxationTime();
	if ( constantRelaxationTime > 0. ) {
		constantRTA = true;
		return;
	}

	highMemory = context.getScatteringMatrixInMemory();

	smearing = DeltaFunction::smearingFactory(context,innerBandStructure);

	// the difference arises from the fact that vectors in the BTE have
	// numCalcs set to nTemp * 3, whereas the matrix only depends on nTemp.
	// so, when we compute the action of the scattering matrix, we must take
	// into account for this misalignment
	if ( highMemory ) {
		numCalcs = statisticsSweep.getNumCalcs();
	} else {
		numCalcs = statisticsSweep.getNumCalcs() * context.getDimensionality();
	}

	// we want to know the state index of acoustic modes at gamma,
	// so that we can set their populations to zero
	for ( long is=0; is<numStates; is++ ) {
		double en = outerBandStructure.getEnergy(is);
		if ( en < 0.1 / ryToCmm1 ) { // cutoff at 0.1 cm^-1
			excludeIndeces.push_back(is);
		}
	}
}

ScatteringMatrix::~ScatteringMatrix() {
	delete smearing;
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
	internalDiagonal(that.internalDiagonal),
	theMatrix(that.theMatrix),
	numStates(that.numStates),
	numPoints(that.numPoints),
	numCalcs(that.numCalcs),
	excludeIndeces(that.excludeIndeces) {
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
		internalDiagonal = that.internalDiagonal;
		theMatrix = that.theMatrix;
		numStates = that.numStates;
		numPoints = that.numPoints;
		numCalcs = that.numCalcs;
		excludeIndeces = that.excludeIndeces;
	}
	return *this;
}

void ScatteringMatrix::setup() {
	// note: here we want to build the matrix or its diagonal
	// builder is a pure virtual function, which is implemented in subclasses
	// c++ discourages calls to pure virtual functions in the constructor

	if ( constantRTA ) return; // nothing to construct

	if ( highMemory ) {
		if ( numCalcs > 1 ) {
			// note: one could write code around this
			// but the methods are very memory intensive for production runs
			Error e("High memory BTE methods can only work with one "
					"temperature and/or chemical potential in a single run");
		}
		theMatrix = Eigen::MatrixXd::Zero(numStates,numStates);
		// calc matrix and linew.
		builder(theMatrix, &internalDiagonal, nullptr, nullptr);
	} else {
		// calc linewidths only
		builder(theMatrix, &internalDiagonal, nullptr, nullptr);
	}
}

VectorBTE ScatteringMatrix::diagonal() {
	if ( constantRTA ) {
		double crt = context.getConstantRelaxationTime();
		VectorBTE diag(statisticsSweep,outerBandStructure,1);
		diag.setConst(1./crt);
		return diag;
	} else {
		return internalDiagonal;
	}
}

Eigen::MatrixXd ScatteringMatrix::dot(const Eigen::MatrixXd & otherMatrix) {
	return theMatrix * otherMatrix;
}

VectorBTE ScatteringMatrix::offDiagonalDot(VectorBTE & inPopulation) {
	if ( highMemory ) {
		// it's just the full matrix product, minus the diagonal contribution
		VectorBTE outPopulation = dot(inPopulation);

		// outPopulation -= internalDiagonal * inPopulation;
		for ( long i=0; i<outPopulation.numCalcs; i++ ) {
			auto [imu,it,idim] = inPopulation.loc2Glob(i);
			auto j = internalDiagonal.glob2Loc(imu,it,DimIndex(0));
			for ( long is=0; is<numStates; is++ ) {
				outPopulation.data(i,is) -= internalDiagonal.data(j,is)
						* inPopulation.data(i,is);
			}
		}
		return outPopulation;
	} else {
		VectorBTE outPopulation(statisticsSweep, outerBandStructure,
				inPopulation.dimensionality);
		builder(theMatrix,nullptr,&inPopulation,&outPopulation);
		// outPopulation = outPopulation - internalDiagonal * inPopulation;
		for ( long i=0; i<outPopulation.numCalcs; i++ ) {
			auto [imu,it,idim] = inPopulation.loc2Glob(i);
			auto j = internalDiagonal.glob2Loc(imu,it,DimIndex(0));
			for ( long is=0; is<numStates; is++ ) {
				outPopulation.data(i,is) -= internalDiagonal.data(j,is)
							* inPopulation.data(i,is);
			}
		}
		return outPopulation;
	}
}

VectorBTE ScatteringMatrix::dot(VectorBTE & inPopulation) {
	if ( highMemory ) {
		VectorBTE outPopulation(statisticsSweep, outerBandStructure,
				inPopulation.dimensionality);
		outPopulation.data.setZero();
		// note: we are assuming that ScatteringMatrix has numCalcs = 1
		for ( auto i1=0; i1<numStates; i1++ ) {
			for ( auto j1=0; j1<numStates; j1++ ) {
				for ( int idim=0; idim<inPopulation.dimensionality; idim++ ) {
					outPopulation.data(idim,i1) +=
							theMatrix(i1,j1) * inPopulation.data(idim,j1);
				}
			}
		}
		// normalization
//		outPopulation.data /= numPoints;
		return outPopulation;
	} else {
		VectorBTE outPopulation(statisticsSweep, outerBandStructure,
				inPopulation.dimensionality);
		builder(theMatrix,nullptr,&inPopulation,&outPopulation);
		return outPopulation;
	}
}

// set,unset the scaling of omega = A/sqrt(bose1*bose1+1)/sqrt(bose2*bose2+1)
void ScatteringMatrix::a2Omega() {

	if ( ! highMemory ) {
		Error e("a2Omega only works if the matrix is stored in memory");
	}

	if ( theMatrix.rows() == 0 ) {
		Error e("The scattering matrix hasn't been built yet");
	}

	if ( isMatrixOmega ) { // it's already with the scaling of omega
		return;
	}

	long iCalc = 0; // as there can only be one temperature

	auto statistics = outerBandStructure.getStatistics();
	auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
	double temp = calcStatistics.temperature;
	double chemPot = calcStatistics.chemicalPotential;

	for ( long ind1=0; ind1<numStates; ind1++ ) {

		if ( std::find(excludeIndeces.begin(), excludeIndeces.end(), ind1)
				!= excludeIndeces.end()) continue;

		double en1 = outerBandStructure.getEnergy(ind1);
		double pop1 = statistics.getPopulation(en1, temp, chemPot);

		double term1;
		if ( statistics.isFermi() ) {
			term1 = pop1 * (1. - pop1);
		} else {
			term1 = pop1 * (1. + pop1);
		}
		internalDiagonal.data(iCalc,ind1) /= term1;

		for ( long ind2=0; ind2<numStates; ind2++ ) {

			if ( std::find(excludeIndeces.begin(), excludeIndeces.end(), ind2)
					!= excludeIndeces.end()) continue;

			double en2 = outerBandStructure.getEnergy(ind2);
			double pop2 = statistics.getPopulation(en2, temp, chemPot);

			double term2;
			if ( statistics.isFermi() ) {
				term2 = pop2 * (1. - pop2);
			} else {
				term2 = pop2 * (1. + pop2);
			}
			theMatrix(ind1,ind2) /= sqrt( term1 * term2 );
		}
	}
	isMatrixOmega = true;
}



// add a flag to remember if we have A or Omega
void ScatteringMatrix::omega2A() {

	if ( ! highMemory ) {
		Error e("a2Omega only works if the matrix is stored in memory");
	}

	if ( theMatrix.rows() == 0 ) {
		Error e("The scattering matrix hasn't been built yet");
	}

	if ( ! isMatrixOmega ) { // it's already with the scaling of A
		return;
	}

	long iCalc = 0; // as there can only be one temperature

	auto statistics = outerBandStructure.getStatistics();
	auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
	double temp = calcStatistics.temperature;
	double chemPot = calcStatistics.chemicalPotential;

	for ( long ind1=0; ind1<numStates; ind1++ ) {

		if ( std::find(excludeIndeces.begin(), excludeIndeces.end(), ind1)
				!= excludeIndeces.end()) continue;

		double en1 = outerBandStructure.getEnergy(ind1);
		double pop1 = statistics.getPopulation(en1, temp, chemPot);

		double term1;
		if ( statistics.isFermi() ) {
			term1 = pop1 * (1. - pop1);
		} else {
			term1 = pop1 * (1. + pop1);
		}
		internalDiagonal.data(iCalc,ind1) *= term1;

		for ( long ind2=0; ind2<numStates; ind2++ ) {

			if ( std::find(excludeIndeces.begin(), excludeIndeces.end(), ind2)
					!= excludeIndeces.end()) continue;

			double en2 = outerBandStructure.getEnergy(ind2);
			double pop2 = statistics.getPopulation(en2, temp, chemPot);

			double term2;
			if ( statistics.isFermi() ) {
				term2 = pop2 * (1. - pop2);
			} else {
				term2 = pop2 * (1. + pop2);
			}
			theMatrix(ind1,ind2) *= sqrt( term1 * term2 );
		}
	}
	isMatrixOmega = false;
}

// to compute the RTA, get the single mode relaxation times
VectorBTE ScatteringMatrix::getSingleModeTimes() {
	if ( constantRTA ) {
		double crt = context.getConstantRelaxationTime();
		VectorBTE diag(statisticsSweep,outerBandStructure,1);
		diag.setConst(crt);
		return diag;
	} else {
		VectorBTE times = internalDiagonal;
		if ( isMatrixOmega ) {
			for ( long iCalc=0; iCalc<internalDiagonal.numCalcs; iCalc++ ) {
				for ( long is=0; is<internalDiagonal.numStates; is++ ) {
					times.data(iCalc,is) = 1. / times.data(iCalc,is);
				}
			}
			return times;
		} else { // A_nu,nu = N(1+-N) / tau
			auto statistics = outerBandStructure.getStatistics();
			for ( long iCalc=0; iCalc<internalDiagonal.numCalcs; iCalc++ ) {
				auto calcStatistics = statisticsSweep.getCalcStatistics(iCalc);
				double temp = calcStatistics.temperature;
				double chemPot = calcStatistics.chemicalPotential;

				for ( long is=0; is<internalDiagonal.numStates; is++ ) {
					double en = outerBandStructure.getEnergy(is);
					double population = statistics.getPopulation(en, temp,
							chemPot);
					if ( statistics.isFermi() ) {
						times.data(iCalc,is) = population * ( 1. - population )
								/ times.data(iCalc,is);
					} else {
						times.data(iCalc,is) = population * ( 1. + population )
								/ times.data(iCalc,is);
					}
				}
			}
			return times;
		}
	}
}

std::tuple<VectorBTE,Eigen::MatrixXd> ScatteringMatrix::diagonalize() {
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(theMatrix);
	Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
	Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();
	// place eigenvalues in an VectorBTE object
	VectorBTE eigvals(statisticsSweep,outerBandStructure,1);
	eigvals.data.row(0) = eigenvalues;

	// correct normalization of eigenvectors
	double volume = outerBandStructure.getPoints().getCrystal(
			).getVolumeUnitCell(context.getDimensionality());
	eigenvectors *= sqrt(numPoints * volume);

	return {eigvals, eigenvectors};
}
