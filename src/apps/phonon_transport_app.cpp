#include "phonon_transport_app.h"
#include "bandstructure.h"
#include "context.h"
#include "drift.h"
#include "exceptions.h"
#include "ifc3_parser.h"
#include "observable.h"
#include "parser.h"
#include "ph_scattering.h"
#include "phonon_thermal_cond.h"
#include "phonon_viscosity.h"
#include "points.h"
#include "specific_heat.h"
#include "wigner_phonon_thermal_cond.h"
#include <iomanip>

void PhononTransportApp::run(Context &context) {

  // Read the necessary input files

  auto tup = Parser::parsePhHarmonic(context);
  auto crystal = std::get<0>(tup);
  auto phononH0 = std::get<1>(tup);

  // first we make compute the band structure on the fine grid

  Points fullPoints(crystal, context.getQMesh());

  if (mpi->mpiHead()) {
    std::cout << "\nComputing phonon band structure." << std::endl;
  }
  auto tup1 = ActiveBandStructure::builder(context, phononH0, fullPoints);
  auto bandStructure = std::get<0>(tup1);
  auto statisticsSweep = std::get<1>(tup1);

  // print some info about state number reduction
  if (mpi->mpiHead()) {
    if(bandStructure.hasWindow() != 0) {
        std::cout << "Window selection reduced phonon band structure from "
                << fullPoints.getNumPoints()*phononH0.getNumBands() << " to "
                << bandStructure.getNumStates() << " states."  << std::endl;
    }
    if(context.getUseSymmetries()) {
      std::cout << "Symmetries reduced phonon band structure from "
          << bandStructure.getNumStates() << " to "
          << bandStructure.irrStateIterator().size() << " states." << std::endl;
    }
    std::cout << "Done computing phonon band structure.\n" << std::endl;
  }

  // load the 3phonon coupling
  auto coupling3Ph = IFC3Parser::parse(context, crystal);

  // build/initialize the scattering matrix and the smearing
  PhScatteringMatrix scatteringMatrix(context, statisticsSweep, bandStructure,
                                      bandStructure, &coupling3Ph, &phononH0);
  scatteringMatrix.setup();

  // solve the BTE at the relaxation time approximation level
  // we always do this, as it's the cheapest solver and is required to know
  // the diagonal for the exact method.

  if (mpi->mpiHead()) {
    std::cout << "\n"
              << std::string(80, '-') << "\n\n"
              << "Solving BTE within the relaxation time approximation."
              << std::endl;
  }

  // compute the phonon populations in the relaxation time approximation.
  // Note: this is the total phonon population n (n != f(1+f) Delta n)

  int dimensionality = context.getDimensionality();
  BulkTDrift drift(statisticsSweep, bandStructure, 3);
  VectorBTE phononRelTimes = scatteringMatrix.getSingleModeTimes();
  VectorBTE popRTA = drift * phononRelTimes;

  // output relaxation times
  scatteringMatrix.outputToJSON("rta_ph_relaxation_times.json");

  // compute the thermal conductivity
  PhononThermalConductivity phTCond(context, statisticsSweep, crystal,
                                    bandStructure);
  phTCond.calcFromPopulation(popRTA);
  phTCond.print();
  phTCond.outputToJSON("rta_phonon_thermal_cond.json");

  // compute the Wigner thermal conductivity
  WignerPhononThermalConductivity phTCondWigner(
      context, statisticsSweep, crystal, bandStructure, phononRelTimes);
  phTCondWigner.calcFromPopulation(popRTA);
  phTCondWigner.print();
  phTCondWigner.outputToJSON("wigner_phonon_thermal_cond.json");

  // compute the thermal conductivity
  PhononViscosity phViscosity(context, statisticsSweep, crystal, bandStructure);
  phViscosity.calcRTA(phononRelTimes);
  phViscosity.print();
  phViscosity.outputToJSON("rta_phonon_viscosity.json");

  // compute the specific heat
  SpecificHeat specificHeat(context, statisticsSweep, crystal, bandStructure);
  specificHeat.calc();
  specificHeat.print();
  specificHeat.outputToJSON("specific_heat.json");

  if (mpi->mpiHead()) {
    std::cout << "\n" << std::string(80, '-') << "\n" << std::endl;
  }

  // if requested, we solve the BTE exactly

  std::vector<std::string> solverBTE = context.getSolverBTE();

  bool doIterative = false;
  bool doVariational = false;
  bool doRelaxons = false;
  for (const auto &s : solverBTE) {
    if (s.compare("iterative") == 0)
      doIterative = true;
    if (s.compare("variational") == 0)
      doVariational = true;
    if (s.compare("relaxons") == 0)
      doRelaxons = true;
  }

  // here we do validation of the input, to check for consistency
  if (doRelaxons && !context.getScatteringMatrixInMemory()) {
    Error("Relaxons require matrix kept in memory");
  }
  if (doRelaxons && context.getUseSymmetries()) {
    Error("Relaxon solver only works without symmetries");
    // Note: this is a problem of the theory I suppose
    // because the scattering matrix may not be anymore symmetric
    // or may require some modifications to make it work
    // that we didn't have yet thought of
  }
  if (doVariational && context.getUseSymmetries()) {
    Error("Variational solver only works without symmetries");
    // Note: this is a problem of the theory I suppose
    // because the scattering matrix may not be anymore symmetric
    // or may require some modifications to make it work
    // that we didn't yet think of
  }
  if (context.getScatteringMatrixInMemory() &&
      statisticsSweep.getNumCalculations() != 1) {
    Error("If scattering matrix is kept in memory, only one "
          "temperature/chemical potential is allowed in a run");
  }

  if (context.getScatteringMatrixInMemory() && !context.getUseSymmetries()) {
    if (doVariational || doRelaxons || doIterative) {
      // reinforce the condition that the scattering matrix is symmetric
      // A -> ( A^T + A ) / 2
      // this helps removing negative eigenvalues which may appear due to noise
      scatteringMatrix.symmetrize();
      // it may not be necessary, so it's commented out
    }
  }

  if (doIterative) {

    if (mpi->mpiHead()) {
      std::cout << "Starting Omini Sparavigna BTE solver\n" << std::endl;
    }

    // initialize the (old) thermal conductivity
    PhononThermalConductivity phTCondOld = phTCond;

    VectorBTE fNext(statisticsSweep, bandStructure, 3);
    VectorBTE sMatrixDiagonal = scatteringMatrix.diagonal();

    // from n, we get f, such that n = bose(bose+1)f
    VectorBTE fRTA = popRTA;
    fRTA.population2Canonical();
    VectorBTE fOld = fRTA;

    auto threshold = context.getConvergenceThresholdBTE();

    for (int iter = 0; iter < context.getMaxIterationsBTE(); iter++) {

      fNext = scatteringMatrix.offDiagonalDot(fOld) / sMatrixDiagonal;
      fNext = fRTA - fNext;

      phTCond.calcFromCanonicalPopulation(fNext);
      phTCond.print(iter);

      // this exit condition must be improved
      // different temperatures might converge differently
      Eigen::Tensor<double,3> newCond = phTCond.getThermalConductivity();
      Eigen::Tensor<double,3> oldCond = phTCondOld.getThermalConductivity();
      double diff = findMaxRelativeDifference(newCond, oldCond);
      if (diff < threshold) {
        break;
      } else {
        phTCondOld = phTCond;
        fOld = fNext;
      }

      if (iter == context.getMaxIterationsBTE() - 1) {
        Error("Reached max BTE iterations without convergence");
      }
    }
    phTCond.print();
    phTCond.outputToJSON("omini_phonon_thermal_cond.json");

    if (mpi->mpiHead()) {
      std::cout << "Finished Omini Sparavigna BTE solver\n\n";
      std::cout << std::string(80, '-') << "\n" << std::endl;
    }
  }

  if (doVariational) {
    if (mpi->mpiHead()) {
      std::cout << "Starting variational BTE solver\n" << std::endl;
    }

    // note: each iteration should take approximately twice as long as
    // the iterative method above (in the way it's written here.

    // initialize the (old) thermal conductivity
    PhononThermalConductivity phTCondOld = phTCond;

    // load the conjugate gradient rescaling factor
    // VectorBTE preconditioning2 = scatteringMatrix.diagonal();
    // VectorBTE preconditioning = preconditioning2.sqrt();

    // initialize b
    VectorBTE b = drift; // / preconditioning;

    // initialize first population guess (the RTA solution)
    VectorBTE f = popRTA;
    f.population2Canonical();

    //  f = f * preconditioning;
    //  b = b / preconditioning;

    // residual
    VectorBTE w0 = scatteringMatrix.dot(f);
    //  VectorBTE w0 = scatteringMatrix.offDiagonalDot(f) + f;
    VectorBTE r = b - w0;

    // search direction
    VectorBTE d = b - w0;

    double threshold = context.getConvergenceThresholdBTE();

    for (int iter = 0; iter < context.getMaxIterationsBTE(); iter++) {
      // execute CG step, as in
      // https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf

      // compute w = E^-1 A E^-1 d
      VectorBTE w = scatteringMatrix.dot(d);
      //  VectorBTE w = scatteringMatrix.offDiagonalDot(d) + d;

      // amount of descent along the search direction
      // size of alpha: (numCalculations,3)
      Eigen::MatrixXd alpha = (r.dot(r)).array() / d.dot(w).array();
      // note: this is to avoid a 0/0 division for low dimensional problems
      for (int i=dimensionality; i<3; i++) {
        for (int iCalc=0; iCalc<statisticsSweep.getNumCalculations(); iCalc++) {
          alpha(iCalc,i) = 0.;
        }
      }

      // new guess of population
      VectorBTE fNew = d * alpha + f;

      // new estimate of residual
      VectorBTE tmp = w * alpha;
      VectorBTE rNew = r - tmp;

      // amount of correction for the search direction
      Eigen::MatrixXd beta = (rNew.dot(rNew)).array() / (r.dot(r)).array();
      // note: this is to avoid a 0/0 division for low dimensional problems
      for (int i=dimensionality; i<3; i++) {
        for (int iCalc=0; iCalc<statisticsSweep.getNumCalculations(); iCalc++) {
          beta(iCalc,i) = 0.;
        }
      }

      // new search direction
      VectorBTE dNew = d * beta + rNew;

      // compute thermal conductivity
      auto aF = scatteringMatrix.dot(fNew);
      //  auto aF = scatteringMatrix.offDiagonalDot(fNew) + fNew;
      phTCond.calcVariational(aF, f, b);

      phTCond.print(iter);

      // decide whether to exit or run the next iteration
      Eigen::Tensor<double,3> newCond = phTCond.getThermalConductivity();
      Eigen::Tensor<double,3> oldCond = phTCondOld.getThermalConductivity();
      double diff = findMaxRelativeDifference(newCond, oldCond);
      if (diff < threshold) {
        break;
      } else {
        phTCondOld = phTCond;
        f = fNew;
        r = rNew;
        d = dNew;
      }

      if (iter == context.getMaxIterationsBTE() - 1) {
        Error("Reached max BTE iterations without convergence");
      }
    }

    // nice formatting of the thermal conductivity at the last step
    phTCond.print();
    phTCond.outputToJSON("variational_phonon_thermal_cond.json");

    if (mpi->mpiHead()) {
      std::cout << "Finished variational BTE solver\n\n";
      std::cout << std::string(80, '-') << "\n" << std::endl;
    }
  }

  if (doRelaxons) {
    if (mpi->mpiHead()) {
      std::cout << "Starting relaxons BTE solver" << std::endl;
    }

    scatteringMatrix.a2Omega();

    auto tup2 = scatteringMatrix.diagonalize();
    auto eigenvalues = std::get<0>(tup2);
    auto eigenvectors = std::get<1>(tup2);
    // EV such that Omega = V D V^-1
    // eigenvectors(phonon index, eigenvalue index)

    phTCond.calcFromRelaxons(context, statisticsSweep, eigenvectors,
                             scatteringMatrix, eigenvalues);
    phTCond.print();
    phTCond.outputToJSON("relaxons_phonon_thermal_cond.json");
    // output relaxation times
    scatteringMatrix.relaxonsToJSON("relaxons_relaxation_times.json",
                                    eigenvalues);

    if (!context.getUseSymmetries()) {
      Vector0 boseEigenvector(statisticsSweep, bandStructure, specificHeat);
      phViscosity.calcFromRelaxons(boseEigenvector, eigenvalues,
                                   scatteringMatrix, eigenvectors);
      phViscosity.print();
      phViscosity.outputToJSON("relaxons_phonon_viscosity.json");
    }

    if (mpi->mpiHead()) {
      std::cout << "Finished relaxons BTE solver\n\n";
      std::cout << std::string(80, '-') << "\n" << std::endl;
    }
  }
  mpi->barrier();
}

void PhononTransportApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getPhFC2FileName(), "PhFC2FileName");
  throwErrorIfUnset(context.getQMesh(), "qMesh");
  throwWarningIfUnset(context.getSumRuleFC2(), "sumRuleFC2");
  throwErrorIfUnset(context.getPhFC3FileName(), "PhFC3FileName");
  throwErrorIfUnset(context.getTemperatures(), "temperatures");
  throwErrorIfUnset(context.getSmearingMethod(), "smearingMethod");
  if (context.getSmearingMethod() == DeltaFunction::gaussian) {
    throwErrorIfUnset(context.getSmearingWidth(), "smearingWidth");
  }
}
