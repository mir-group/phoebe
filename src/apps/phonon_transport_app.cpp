#include <string>
#include "phonon_transport_app.h"
#include "context.h"
#include "constants.h"
#include "exceptions.h"
#include "statistics.h"
#include "io.h"
#include "vector_bte.h"
#include "drift.h"
#include "scattering.h"
#include "observable.h"

void PhononTransportApp::run(Context & context) {

	// Read the necessary input files

	auto [crystal, phononH0] = setupPhononH0(context);

	// first we make compute the band structure on the fine grid

	// here we branch into two versions:

	FullPoints fullQPoints(crystal, context.getQMesh());
	bool withVelocities = true;
	bool withEigenvectors = true;
	Statistics statistics = phononH0.getStatistics();
	FullBandStructure fullPhBandStructure(phononH0.getNumBands(), statistics,
			withVelocities, withEigenvectors, &fullQPoints);
	fullPhBandStructure.populate(phononH0);

	auto [fullPoints,fullBandStructure] = buildFullBandStructure(crystal,
			context.getQMesh(), phononH0, withVelocities, withEigenvectors);

	// then we apply a filter to retain only useful energies
	auto [activePoints, activeBandStructure] = restrictBandStructure(context,
			fullBandStructure);

	// free up some memory
	fullBandStructure.~FullBandStructure();
	phononH0.~PhononH0();

	// set the chemical potentials to zero, load temperatures
	context.getTemperatures();
	Eigen::VectorXd chemicalPotentials(1);
	chemicalPotentials.setZero();

	// BTE solver

	// Let's first implement this here. then we will move it into a class
	// that can be shared for both electron and phonon (separate) BTE's
	// but just with the operators specialized to various cases

	// initialize populations

	BulkTDrift drift(context, activeBandStructure);

	std::vector<std::string> solverBTE = context.getSolverBTE();

	ScatteringMatrix scatteringMatrix(context);

	PhononThermalConductivity phTCond(context, crystal);
	PhononThermalConductivity phTCondOld(context, crystal);

	VectorBTE sMatrixDiagonal = scatteringMatrix.diagonal();

	VectorBTE popRTA = drift / sMatrixDiagonal;
	phTCond.calcFromPopulation(popRTA, drift);

	if ( std::find(solverBTE.begin(),solverBTE.end(),"iterative") !=
			solverBTE.end() ) {
		// recompute scatteringoperator everytime
		VectorBTE popNext(context, activeBandStructure);
		VectorBTE popOld(context, activeBandStructure);

		popOld = popRTA;

		double threshold = context.getConvergenceThresholdBTE();

		for (  long iter=0; iter<context.getMaxIterationsBTE(); iter++ ) {
			VectorBTE popNext = popRTA
					- scatteringMatrix.offDiagonalDot(popOld)
					/ sMatrixDiagonal;

			phTCond.calcFromPopulation(popNext, drift);

			// this exit condition must be improved
			// different temperatures might converge differently
			auto diff = phTCond - phTCondOld;
			if ( diff.getNorm().maxCoeff() < threshold ) {
				break;
			} else {
				phTCondOld = phTCond;
				popOld = popNext;
			}
		}
	}

	if ( std::find(solverBTE.begin(),solverBTE.end(),"minimization") !=
				solverBTE.end() ) {
		VectorBTE popNext(context, activeBandStructure);
		VectorBTE popOld(context, activeBandStructure);

		// note that we use the rescaled conjugate gradient method
		VectorBTE bTilde = drift / sMatrixDiagonal.sqrt();
		VectorBTE fRTA = bTilde;

		scatteringMatrix.setCGScaling();

		VectorBTE fOld = fRTA;
		VectorBTE gOld = scatteringMatrix.dot(fRTA) - fOld;
		VectorBTE hOld = - gOld;

		double threshold = context.getConvergenceThresholdBTE();

		for ( long iter=0; iter<context.getMaxIterationsBTE(); iter++ ) {
			VectorBTE tOld = scatteringMatrix.dot(hOld);

			Eigen::VectorXd alpha = ( gOld * hOld ).array() / ( hOld * tOld ).array();

			VectorBTE fNew = fOld - hOld * alpha;
			VectorBTE gNew = gOld - tOld * alpha;

			Eigen::VectorXd beta = ( gNew * gNew ).array() / ( gOld * gOld ).array();
			VectorBTE hNew = - gNew + hOld * beta;

			// action of scattering operator on vectors of populations:
			VectorBTE AfNew = gNew + bTilde;

			phTCond.calcVariational(AfNew, fNew, bTilde);

			auto diff = phTCond - phTCondOld;
			if ( diff.getNorm().maxCoeff() < threshold ) {
				break;
			} else {
				phTCondOld = phTCond;
				fOld = fNew;
			    gNew = gNew;
			    hOld = hNew;
			}
		}
	}
};
