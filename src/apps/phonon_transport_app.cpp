#include <string>
#include "phonon_transport_app.h"
#include "context.h"
#include "constants.h"
#include "exceptions.h"
#include "statistics.h"
#include "io.h"
#include "vector_bte.h"
#include "drift.h"
#include "ph_scattering.h"
#include "observable.h"
#include "bandstructure.h"
#include "ifc3_parser.h"

void PhononTransportApp::run(Context & context) {

	// Read the necessary input files

	auto [crystal, phononH0] = parser.parsePhHarmonic(context);
	phononH0.setAcousticSumRule(context.getSumRuleD2());

	// first we make compute the band structure on the fine grid

	FullPoints fullPoints(crystal, context.getQMesh());
	bool withVelocities = true;
	bool withEigenvectors = true;
	FullBandStructure<FullPoints> bandStructure = phononH0.populate(
			fullPoints, withVelocities, withEigenvectors);

	// set the chemical potentials to zero, load temperatures
	StatisticsSweep statisticsSweep(context);

	// load the 3phonon coupling
	auto coupling3Ph = IFC3Parser::parse(context, crystal);

	// initialize populations
	BulkTDrift drift(context, bandStructure);

	// build/initialize the scattering matrix and the smearing
	PhScatteringMatrix scatteringMatrix(context, statisticsSweep,
			bandStructure, bandStructure, &coupling3Ph);
	scatteringMatrix.setup();

	// solve the BTE at the relaxation time approximation level
	// we always do this, as it's the cheapest solver and is required to know
	// the diagonal for the exact method.

	std::cout << "\n";
	std::cout << "Solving BTE within the relaxation time approximation.\n";

	VectorBTE sMatrixDiagonal = scatteringMatrix.diagonal();

	// compute the populations
	VectorBTE popRTA = drift / sMatrixDiagonal;

	// compute the thermal conductivity
	PhononThermalConductivity phTCond(context, crystal);
	phTCond.calcFromPopulation(popRTA, drift);
	phTCond.print();

	// if needed, we solve the BTE exactly

	std::vector<std::string> solverBTE = context.getSolverBTE();
	if ( std::find(solverBTE.begin(),solverBTE.end(),"iterative") !=
			solverBTE.end() ) {

		// initialize the (old) thermal conductivity
		PhononThermalConductivity phTCondOld(context, crystal);

		// recompute scatteringoperator everytime
		VectorBTE popNext(context, bandStructure);
		VectorBTE popOld(context, bandStructure);

		popOld = popRTA;

		auto threshold = context.getConvergenceThresholdBTE();

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
		// initialize the (old) thermal conductivity
		PhononThermalConductivity phTCondOld(context, crystal);

		// initialize two vectors of populations at current and former step
		VectorBTE popNext(context, bandStructure);
		VectorBTE popOld(context, bandStructure);

		// note that we use the rescaled conjugate gradient method
		auto bTilde = drift / sMatrixDiagonal.sqrt();
		auto fRTA = bTilde;

		scatteringMatrix.setCGScaling();

		// do the conjugate gradient method for thermal conductivity.
		auto fOld = fRTA;
		auto gOld = scatteringMatrix.dot(fRTA) - fOld;
		auto hOld = - gOld;

		double threshold = context.getConvergenceThresholdBTE();

		for ( long iter=0; iter<context.getMaxIterationsBTE(); iter++ ) {
			auto tOld = scatteringMatrix.dot(hOld);

			Eigen::VectorXd alpha = ( gOld * hOld ).array()
					/ ( hOld * tOld ).array();

			auto fNew = fOld - hOld * alpha;
			auto gNew = gOld - tOld * alpha;

			Eigen::VectorXd beta = ( gNew * gNew ).array()
					/ ( gOld * gOld ).array();
			auto hNew = - gNew + hOld * beta;

			// action of scattering operator on vectors of populations:
			auto AfNew = gNew + bTilde;

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
