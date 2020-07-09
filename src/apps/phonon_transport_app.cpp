#include "phonon_transport_app.h"
#include "bandstructure.h"
#include "constants.h"
#include "context.h"
#include "drift.h"
#include "exceptions.h"
#include "ifc3_parser.h"
#include "io.h"
#include "observable.h"
#include "particle.h"
#include "ph_scattering.h"
#include "phonon_thermal_cond.h"
#include "wigner_phonon_thermal_cond.h"
#include "phonon_viscosity.h"
#include "qe_input_parser.h"
#include "specific_heat.h"
#include <iomanip>

void PhononTransportApp::run(Context &context) {

  // Read the necessary input files

  auto [crystal, phononH0] = QEParser::parsePhHarmonic(context);

  // first we make compute the band structure on the fine grid

  FullPoints fullPoints(crystal, context.getQMesh());

  //	bool withVelocities = true;
  //	bool withEigenvectors = true;
  //	FullBandStructure bandStructure = phononH0.populate(
  //			fullPoints, withVelocities, withEigenvectors);
  //	// set the chemical potentials to zero, load temperatures
  //	StatisticsSweep statisticsSweep(context);

  auto [bandStructure, statisticsSweep] =
      ActiveBandStructure::builder(context, phononH0, fullPoints);

  // load the 3phonon coupling
  auto coupling3Ph = IFC3Parser::parse(context, crystal);

  // build/initialize the scattering matrix and the smearing
  PhScatteringMatrix scatteringMatrix(context, statisticsSweep, bandStructure,
                                      bandStructure, &coupling3Ph, &phononH0);
  scatteringMatrix.setup();

  // solve the BTE at the relaxation time approximation level
  // we always do this, as it's the cheapest solver and is required to know
  // the diagonal for the exact method.

  if ( mpi->mpiHead()) {
    std::cout << "\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << "\n";
    std::cout << "Solving BTE within the relaxation time approximation.\n";
  }

  // compute the phonon populations in the relaxation time approximation.
  // Note: this is the total phonon population n (n != f(1+f) Delta n)

  auto dimensionality = context.getDimensionality();
  BulkTDrift drift(statisticsSweep, bandStructure, dimensionality);
  VectorBTE phononRelTimes = scatteringMatrix.getSingleModeTimes();
  VectorBTE popRTA = drift * phononRelTimes;

  // compute the thermal conductivity
  PhononThermalConductivity phTCond(statisticsSweep, crystal, bandStructure);
  phTCond.calcFromPopulation(popRTA);
  phTCond.print();

  // compute the Wigner thermal conductivity
  WignerPhononThermalConductivity phTCondWigner(statisticsSweep, crystal,
                                                bandStructure, phononRelTimes);
  phTCondWigner.calcFromPopulation(popRTA);
  phTCondWigner.print();

  // compute the thermal conductivity
  PhononViscosity phViscosity(statisticsSweep, crystal, bandStructure);
  phViscosity.calcRTA(phononRelTimes);
  phViscosity.print();

  // compute the specific heat
  SpecificHeat specificHeat(statisticsSweep, crystal, bandStructure);
  specificHeat.calc();
  specificHeat.print();

  if ( mpi->mpiHead()) {
    std::cout << "\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << "\n";
  }

  // if requested, we solve the BTE exactly

  std::vector<std::string> solverBTE = context.getSolverBTE();

  bool doIterative = false;
  bool doVariational = false;
  bool doRelaxons = false;
  for (auto s : solverBTE) {
    if (s.compare("iterative") == 0)
      doIterative = true;
    if (s.compare("variational") == 0)
      doVariational = true;
    if (s.compare("relaxons") == 0)
      doRelaxons = true;
  }

  // here we do validation of the input, to check for consistency
  if (doRelaxons && !context.getScatteringMatrixInMemory()) {
    Error e("Relaxons require matrix kept in memory");
  }
  if (context.getScatteringMatrixInMemory() &&
      statisticsSweep.getNumCalcs() != 1) {
    Error e("If scattering matrix is kept in memory, only one "
            "temperature/chemical potential is allowed in a run");
  }

  mpi->barrier();

  if (doIterative) {

    if ( mpi->mpiHead()) {
      std::cout << "Starting Omini Sparavigna BTE solver\n";
      std::cout << "\n";
    }

    // initialize the (old) thermal conductivity
    PhononThermalConductivity phTCondOld = phTCond;

    VectorBTE fNext(statisticsSweep, bandStructure, dimensionality);
    VectorBTE sMatrixDiagonal = scatteringMatrix.diagonal();

    // from n, we get f, such that n = bose(bose+1)f
    VectorBTE fRTA = popRTA;
    fRTA.population2Canonical();
    VectorBTE fOld = fRTA;

    auto threshold = context.getConvergenceThresholdBTE();

    for (long iter = 0; iter < context.getMaxIterationsBTE(); iter++) {

      fNext = scatteringMatrix.offDiagonalDot(fOld) / sMatrixDiagonal;
      fNext = fRTA - fNext;

      phTCond.calcFromCanonicalPopulation(fNext);
      phTCond.print(iter);

      // this exit condition must be improved
      // different temperatures might converge differently
      auto diff = phTCond - phTCondOld;
      if (diff.getNorm().maxCoeff() < threshold) {
        break;
      } else {
        phTCondOld = phTCond;
        fOld = fNext;
      }

      if (iter == context.getMaxIterationsBTE() - 1) {
        Error e("Reached max BTE iterations without convergence");
      }
    }
    phTCond.print();
    if ( mpi->mpiHead()) {
      std::cout << "Finished Omini Sparavigna BTE solver\n";
      std::cout << "\n";
      std::cout << std::string(80, '-') << "\n";
      std::cout << "\n";
    }
  }

  if (doVariational) {
    if ( mpi->mpiHead()) {
      std::cout << "Starting variational BTE solver\n";
      std::cout << "\n";
    }

    // note: each iteration should take approximately twice as long as
    // the iterative method above (in the way it's written here.

    // initialize the (old) thermal conductivity
    PhononThermalConductivity phTCondOld = phTCond;

    // load the conjugate gradient rescaling factor
    VectorBTE sMatrixDiagonalSqrt = scatteringMatrix.diagonal().sqrt();
    VectorBTE sMatrixDiagonal = scatteringMatrix.diagonal();

    // set the initial guess to the RTA solution
    VectorBTE fNew = popRTA;
    // from n, we get f, such that n = bose(bose+1)f
    fNew.population2Canonical();
    // CG rescaling
    fNew = fNew * sMatrixDiagonalSqrt; // CG scaling

    // save the population of the previous step
    VectorBTE fOld = fNew;

    // do the conjugate gradient method for thermal conductivity.
    //		auto gOld = scatteringMatrix.dot(fNew) - fOld;
    auto gOld = scatteringMatrix.dot(fNew);
    gOld = gOld / sMatrixDiagonal; // CG scaling
    gOld = gOld - fOld;
    auto hOld = -gOld;

    double threshold = context.getConvergenceThresholdBTE();

    for (long iter = 0; iter < context.getMaxIterationsBTE(); iter++) {
      // execute CG step, as in

      auto tOld = scatteringMatrix.dot(hOld);
      tOld = tOld / sMatrixDiagonal; // CG scaling

      Eigen::VectorXd alpha =
          (gOld.dot(hOld)).array() / (hOld.dot(tOld)).array();

      fNew = hOld * alpha;
      fNew = fOld - fNew;

      auto gNew = tOld * alpha;
      gNew = gOld - gNew;

      Eigen::VectorXd beta =
          (gNew.dot(gNew)).array() / (gOld.dot(gOld)).array();
      auto hNew = hOld * beta;
      hNew = -gNew + hNew;

      // calculate the thermal conductivity
      auto tmpF = scatteringMatrix.dot(fNew);
      phTCond.calcVariational(tmpF, fNew, sMatrixDiagonalSqrt);
      phTCond.print(iter);

      // decide whether to exit or run the next iteration
      auto diff = phTCond - phTCondOld;
      if (diff.getNorm().maxCoeff() < threshold) {
        break;
      } else {
        phTCondOld = phTCond;
        fOld = fNew;
        gOld = gNew;
        hOld = hNew;
      }

      if (iter == context.getMaxIterationsBTE() - 1) {
        Error e("Reached max BTE iterations without convergence");
      }
    }

    // nice formatting of the thermal conductivity at the last step
    phTCond.print();

    if ( mpi->mpiHead()) {
      std::cout << "Finished variational BTE solver\n";
      std::cout << "\n";
      std::cout << std::string(80, '-') << "\n";
      std::cout << "\n";
    }
  }

  if (doRelaxons) {
    if ( mpi->mpiHead()) {
      std::cout << "Starting relaxons BTE solver\n";
    }
    scatteringMatrix.a2Omega();
    auto [eigenvalues, eigenvectors] = scatteringMatrix.diagonalize();
    // EV such that Omega = V D V^-1
    // eigenvectors(phonon index, eigenvalue index)

    Vector0 boseEigenvector(statisticsSweep, bandStructure, specificHeat);

    VectorBTE relaxonV(statisticsSweep, bandStructure, 3);
    for (long iCalc = 0; iCalc < relaxonV.numCalcs; iCalc++) {

      double norm = 1. /
                    crystal.getVolumeUnitCell(context.getDimensionality()) /
                    bandStructure.getNumPoints(true);

      auto [imu, it, idim] = relaxonV.loc2Glob(iCalc);
      int idimIndex = idim.get();
      auto jCalc = boseEigenvector.glob2Loc(imu, it, DimIndex(0));
      for (long is = 0; is < bandStructure.getNumStates(); is++) {
        Eigen::Vector3d v = bandStructure.getGroupVelocity(is);
        relaxonV.data(iCalc, is) =
            boseEigenvector.data(jCalc, is) * norm * v(idimIndex);
      }
    }
    relaxonV = relaxonV * eigenvectors;

    VectorBTE relaxationTimes = eigenvalues.reciprocal();
    phTCond.calcFromRelaxons(specificHeat, relaxonV, relaxationTimes);
    phTCond.print();

    phViscosity.calcFromRelaxons(boseEigenvector, relaxationTimes,
                                 scatteringMatrix, eigenvectors);
    phViscosity.print();

    if ( mpi->mpiHead()) {
      std::cout << "Finished relaxons BTE solver\n";
      std::cout << "\n";
      std::cout << std::string(80, '-') << "\n";
      std::cout << "\n";
    }
  }
  mpi->barrier();
}

void PhononTransportApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getPhD2FileName(), "PhD2FileName");
  throwErrorIfUnset(context.getQMesh(), "qMesh");
  throwWarningIfUnset(context.getSumRuleD2(), "sumRuleD2");
  throwErrorIfUnset(context.getPhD3FileName(), "PhD3FileName");
  throwErrorIfUnset(context.getTemperatures(), "temperatures");
  throwErrorIfUnset(context.getSmearingMethod(), "smearingMethod");
  throwErrorIfUnset(context.getSmearingWidth(), "smearingWidth");
}
