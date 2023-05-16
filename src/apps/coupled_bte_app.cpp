#include "coupled_bte_app.h"
#include "bandstructure.h"
#include "context.h"
#include "drift.h"
#include "exceptions.h"
#include "ifc3_parser.h"
#include "observable.h"
#include "parser.h"
//#include "phel_scattering_matrix.h"
//#include "ph_scattering_matrix.h"
#include "c_scattering_matrix.h"
//#include "phonon_thermal_cond.h"
//#include "phonon_viscosity.h"
#include "points.h"
#include "specific_heat.h"
//#include "wigner_phonon_thermal_cond.h"
#include <iomanip>

void CoupledTransportApp::run(Context &context) {

  // there are four major possible contributions to this application
  // electron-phonon, phonon-phonon, phonon-electron, and electron-electron
  // Here we set these booleans based on if their relevant files are set
  // If elph files are available, we automatically calculate both phel and elph

  //bool useElElInteraction = !context.getElElFileName().empty(); // not implemented
  bool useElPhInteraction = !context.getElphFileName().empty();
  bool usePhPhInteraction = !context.getPhFC3FileName().empty();

  bool useCRTA = std::isnan(context.getConstantRelaxationTime());
  if(useCRTA) { // this isn't a possibility here
    Error("CRTA is not an option for the coupled BTE app!");
  }

  if(!usePhPhInteraction && !useElPhInteraction) {
    Error("To run the coupled BTE app supply a ph-ph or el-ph file!");
  }

  // TODO add an error message if anything but the population window is used
  // TODO only 1 num calculation at a time

  // Set up phonon bandstructure information ---------------------------------------------
  // Read the necessary input files
  auto tup = Parser::parsePhHarmonic(context);
  auto crystal = std::get<0>(tup);
  auto phononH0 = std::get<1>(tup);

  auto t1 = Parser::parseElHarmonicWannier(context, &crystal);
  auto crystalEl = std::get<0>(t1);
  auto electronH0 = std::get<1>(t1);

  // first we make compute the band structure on the fine grid
  Points fullPoints(crystal, context.getQMesh());

  if (mpi->mpiHead()) {
    std::cout << "\nComputing phonon band structure." << std::endl;
  }
  auto tup1 = ActiveBandStructure::builder(context, phononH0, fullPoints);
  auto phBandStructure = std::get<0>(tup1);
  auto statisticsSweep = std::get<1>(tup1);

  // stop the code if someone tries to run it with more than one value of (mu, T)
  if(statisticsSweep.getNumChemicalPotentials() != 1) { 
      Error("Can't run coupled BTE solve with more than one chemical potential or temperature "
        "at a time, as this would take up far too much memory at once!"); 
  }

  // print some info about state number reduction
  if (mpi->mpiHead()) {
    if(phBandStructure.hasWindow() != 0) {
        std::cout << "Window selection reduced phonon band structure from "
                << fullPoints.getNumPoints() * phononH0.getNumBands() << " to "
                << phBandStructure.getNumStates() << " states."  << std::endl;
    }
    if(context.getUseSymmetries()) {
      std::cout << "Symmetries reduced phonon band structure from "
          << phBandStructure.getNumStates() << " to "
          << phBandStructure.irrStateIterator().size() << " states." << std::endl;
    }
    std::cout << "Done computing phonon band structure.\n" << std::endl;
  }

  // Set up electron bandstructure information ---------------------------------------------
  if (mpi->mpiHead()) {
    std::cout << "\nComputing electronic band structure.\n"
              << std::endl;
  }

  // TODO we have to resolve having these two "fullPoints" meshes, as
  // we definitely need to allow for one phonon and one el one, since they converge
  // differently 
  // construct electronic band structure
  // TODO additionally, are the two statistics sweeps separate?

  //Points kPoints(crystal, context.getKMesh());
  auto t3 = ActiveBandStructure::builder(context, electronH0, fullPoints);
  auto elBandStructure = std::get<0>(t3);
  //auto statisticsSweep = std::get<1>(t3);

  // don't proceed if we use more than one doping concentration --
  // phph scattering only has 1 mu value, therefore the linewidths won't add to it correctly
  int numMu = statisticsSweep.getNumChemicalPotentials();
  if (numMu != 1) {
      Error("Can currenly only add ph-el scattering one doping "
        "concentration at the time. Let us know if you want to have multiple mu values as a feature.");
  }

    // print some info about how window and symmetries have reduced el bands
    if (mpi->mpiHead()) {
      if(elBandStructure.hasWindow() != 0) {
          std::cout << "Window selection reduced electronic band structure from "
                  << fullPoints.getNumPoints() * electronH0.getNumBands() << " to "
                  << elBandStructure.getNumStates() << " states."  << std::endl;
      }
      if(context.getUseSymmetries()) {
        std::cout << "Symmetries reduced electronic band structure from "
          << elBandStructure.getNumStates() << " to "
          << elBandStructure.irrStateIterator().size() << " states." << std::endl;
      }
      std::cout << "Done computing electronic band structure.\n" << std::endl;
    }

  // read in phonon-phonon coupling -------------------------------------
  // this also checks that the crystal is the same one read in for 3ph
  Interaction3Ph coupling3Ph = IFC3Parser::parse(context, crystal);

  // read in electron-phonon coupling ---------------------------------
  // this also checks that the crystal is the same one read in for 3ph
  InteractionElPhWan couplingElPh = 
                         InteractionElPhWan::parse(context, crystal, &phononH0);

  // Construct the full C matrix
  // the dimensions of this matrix are (numElStates + numPhStates, numElStates + numPhStates) 
  // We definitely need to supply both the phonon and electron bandstructures to this object. 
  // The phel scattering probably needs it's own separate electronic bandstructure object, 
  // as it requires much denser sampling to converge and we need to allow for that 

  // TODO create this object
  // TODO feels like we should be passing these thigns by reference? 
  CoupledScatteringMatrix scatteringMatrix(context, statisticsSweep, 
                                        elBandStructure, phBandStructure, 
                                        &coupling3Ph, &couplingElPh, 
                                        &phononH0, &electronH0);  
  scatteringMatrix.setup();   // adds in all the scattering rates
  
  // alternatively, we may want to add contributions one at a time to this matrix, 
  // especially if we're struggling in memory, so that we can only store either
  // the elph matrix elements or ph-ph matrix elements at once
  // This shouldn't be too hard -- instead of calling builder, we should be able to structure
  // things so it's possible to add contribution separately

  // load ph-el interaction and add in these linewidths
  if (mpi->mpiHead()) {
    std::cout << "\nStarting phonon-electron scattering calculation." << std::endl;
  }
  // helper function to return ph-el linewidths
  // TODO would it be better to make another friend function that just 
  // calculates these linewidths? 
  //VectorBTE phElLinewidths = getPhononElectronLinewidth(context, crystal,
  //                                                      phBandStructure, phononH0);

  // if we're using both phel and phph times, we should output
  // each independent linewidth set. PhEl is output above.
  scatteringMatrix.outputToJSON("rta_phph_relaxation_times.json");

  // add in the phel linewidths -- use getLinewidths to remove pop factors
  // TODO this needs to be changed so that we can add it to only one quadrant
  //VectorBTE totalRates = scatteringMatrix.getLinewidths();
  //totalRates = totalRates + phElLinewidths;
  //scatteringMatrix.setLinewidths(totalRates);


  // BTE Solvers ---------------------------------------------


  // solve the BTE at the relaxation time approximation level
  // we always do this, as it's the cheapest solver and is required to know
  // the diagonal for the exact method.
  if (mpi->mpiHead()) {
    std::cout << "\n" << std::string(80, '-') << "\n\n"
              << "Solving BTE within the relaxation time approximation."
              << std::endl;
  }

  // compute the phonon populations in the relaxation time approximation.
  // Note: this is the total phonon population n (n != f(1+f) Delta n)

  // TODO these feel like they need to be swaped into an option that can take both kinds of bandstructures? 
/*
  BulkTDrift drift(statisticsSweep, phBandStructure, 3);
  VectorBTE phononRelTimes = scatteringMatrix.getSingleModeTimes();
  VectorBTE popRTA = drift * phononRelTimes;

  // output relaxation times
  scatteringMatrix.outputToJSON("rta_ph_relaxation_times.json");

  // compute the thermal conductivity
  PhononThermalConductivity phTCond(context, statisticsSweep, crystal, phBandStructure);
  phTCond.calcFromPopulation(popRTA);
  phTCond.print();
  phTCond.outputToJSON("rta_phonon_thermal_cond.json");

  // compute the Wigner thermal conductivity
  WignerPhononThermalConductivity phTCondWigner(context, statisticsSweep, crystal, phBandStructure, phononRelTimes);
  phTCondWigner.calcFromPopulation(popRTA);
  phTCondWigner.print();
  phTCondWigner.outputToJSON("wigner_phonon_thermal_cond.json");

  // compute the thermal conductivity
  PhononViscosity phViscosity(context, statisticsSweep, crystal, phBandStructure);
  phViscosity.calcRTA(phononRelTimes);
  phViscosity.print();
  phViscosity.outputToJSON("rta_phonon_viscosity.json");

  // compute the specific heat
  SpecificHeat specificHeat(context, statisticsSweep, crystal, phBandStructure);
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
    if (s.compare("iterative") == 0) {   doIterative = true; }
    if (s.compare("variational") == 0) {  doVariational = true; }
    if (s.compare("relaxons") == 0) {    doRelaxons = true; }
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
      if ( context.getSymmetrizeMatrix() ) {
        // reinforce the condition that the scattering matrix is symmetric
        // A -> ( A^T + A ) / 2
        // this helps removing negative eigenvalues which may appear due to noise
        scatteringMatrix.symmetrize();
        // it may not be necessary, so it's commented out
      }
    }
  }

  if (doIterative) {

    if (mpi->mpiHead()) {
      std::cout << "Starting Omini Sparavigna BTE solver\n" << std::endl;
    }

    // initialize the (old) thermal conductivity
    PhononThermalConductivity phTCondOld = phTCond;

    VectorBTE fNext(statisticsSweep, phBandStructure, 3);
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

      // size of alpha: (numCalculations,3)
//      Eigen::MatrixXd alpha = (r.dot(r)).array() / d.dot(w).array();
//      // note: this is to avoid a 0/0 division for low dimensional problems
//      for (int i=dimensionality; i<3; i++) {
//        for (int iCalc=0; iCalc<statisticsSweep.getNumCalculations(); iCalc++) {
//          alpha(iCalc,i) = 0.;
//        }
//      }

      // amount of descent along the search direction
      int numCalculations = statisticsSweep.getNumCalculations();
      Eigen::MatrixXd alpha = Eigen::MatrixXd::Zero(numCalculations,3);
      {
        Eigen::MatrixXd num = r.dot(r);
        Eigen::MatrixXd den = d.dot(w);
        for (int iCalc=0; iCalc<numCalculations; iCalc++) {
          for (int i : {0,1,2}) {
            if (den(iCalc,i) != 0.) {
              alpha(iCalc, i) = num(iCalc, i) / den(iCalc, i);
            }
          }
        }
      }

      // new guess of population
      VectorBTE fNew = d * alpha + f;

      // new estimate of residual
      VectorBTE tmp = w * alpha;
      VectorBTE rNew = r - tmp;

      // amount of correction for the search direction
      // Eigen::MatrixXd beta = (rNew.dot(rNew)).array() / (r.dot(r)).array();
      Eigen::MatrixXd beta = Eigen::MatrixXd::Zero(numCalculations,3);
      {
        Eigen::MatrixXd num = rNew.dot(rNew);
        Eigen::MatrixXd den = r.dot(r);
        for (int iCalc=0; iCalc<numCalculations; iCalc++) {
          for (int i : {0,1,2}) {
            if (den(iCalc,i) != 0.) {
              beta(iCalc, i) = num(iCalc, i) / den(iCalc, i);
            }
          }
        }
      }

//      // note: this is to avoid a 0/0 division for low dimensional problems
//      for (int i=dimensionality; i<3; i++) {
//        for (int iCalc=0; iCalc<statisticsSweep.getNumCalculations(); iCalc++) {
//          beta(iCalc,i) = 0.;
//        }
//      }

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
      Vector0 boseEigenvector(statisticsSweep, phBandStructure, specificHeat);
      phViscosity.calcFromRelaxons(boseEigenvector, eigenvalues,
                                   scatteringMatrix, eigenvectors);
      phViscosity.print();
      phViscosity.outputToJSON("relaxons_phonon_viscosity.json");
    }

    if (mpi->mpiHead()) {
      std::cout << "Finished relaxons BTE solver\n\n";
      std::cout << std::string(80, '-') << "\n" << std::endl;
    }
  }*/
  mpi->barrier();
}
/*
// helper function to generate phEl rates
VectorBTE CoupledTransportApp::getPhononElectronLinewidth(Context& context, Crystal& crystal,
                                                         ActiveBandStructure &phBandStructure,
                                                         PhononH0& phononH0) {

    // load electron band structure
    auto t1 = Parser::parseElHarmonicWannier(context, &crystal);
    auto crystalEl = std::get<0>(t1);
    auto electronH0 = std::get<1>(t1);

    // check that the crystal in the elph calculation is the
    // same as the one in the phph calculation
    if (crystal.getDirectUnitCell() != crystalEl.getDirectUnitCell()) {
      Error("Phonon-electrons scattering requested, but crystals used for ph-ph and \n"
        "ph-el scattering are not the same!");
    }

    // load the elph coupling
    // Note: this file contains the number of electrons
    // which is needed to understand where to place the fermi level
    auto couplingElPh = InteractionElPhWan::parse(context, crystal, &phononH0);

    // compute the band structure on the fine grid
    if (mpi->mpiHead()) {
      std::cout << "\nComputing electronic band structure for ph-el scattering.\n"
                << std::endl;
    }

    // update the window so that it's only a very narrow
    // scale slice around mu, 1.25* the max phonon energy
    double maxPhEnergy = phBandStructure.getMaxEnergy();
    auto inputWindowType = context.getWindowType();
    context.setWindowType("muCenteredEnergy");
    if(mpi->mpiHead()) {
      std::cout << "Of the active phonon modes, the maximum energy state is " <<
          maxPhEnergy*energyRyToEv*1e3 << " meV." <<
          "\nSelecting states within +/- 1.25 x " << maxPhEnergy*energyRyToEv*1e3 << " meV"
          << " of max/min electronic mu values." << std::endl;
    }
    Eigen::Vector2d range = {-1.25*maxPhEnergy,1.25*maxPhEnergy};
    context.setWindowEnergyLimit(range);

    // construct electronic band structure
    Points fullPoints(crystal, context.getKMesh());
    auto t3 = ActiveBandStructure::builder(context, electronH0, fullPoints);
    auto elBandStructure = std::get<0>(t3);
    auto statisticsSweep = std::get<1>(t3);

    // don't proceed if we use more than one doping concentration --
    // phph scattering only has 1 mu value, therefore the linewidths won't add to it correctly
    int numMu = statisticsSweep.getNumChemicalPotentials();
    if (numMu != 1) {
        Error("Can currenly only add ph-el scattering one doping "
          "concentration at the time. Let us know if you want to have multiple mu values as a feature.");
    }

    // print some info about how window and symmetries have reduced el bands
    if (mpi->mpiHead()) {
      if(elBandStructure.hasWindow() != 0) {
          std::cout << "Window selection reduced electronic band structure from "
                  << fullPoints.getNumPoints() * electronH0.getNumBands() << " to "
                  << elBandStructure.getNumStates() << " states."  << std::endl;
      }
      if(context.getUseSymmetries()) {
        std::cout << "Symmetries reduced electronic band structure from "
          << elBandStructure.getNumStates() << " to "
          << elBandStructure.irrStateIterator().size() << " states." << std::endl;
      }
      std::cout << "Done computing electronic band structure.\n" << std::endl;
    }

    // build/initialize the scattering matrix and the smearing
    // it shouldn't take up too much memory, as it's just
    // the diagonal -- be careful to put elBandStruct first
    PhElScatteringMatrix phelScatteringMatrix(context, statisticsSweep,
                                              elBandStructure, phBandStructure,
                                              couplingElPh, electronH0);

    if(int(elBandStructure.irrStateIterator().size()) < mpi->getSize()) {
      Error("Cannot calculate PhEl scattering matrix with fewer states than\n"
        "number of MPI processes. Reduce the number of MPI processes,\n"
        "perhaps use more OMP threads instead.");
    }
    phelScatteringMatrix.setup();
    phelScatteringMatrix.outputToJSON("rta_phel_relaxation_times.json");

    // important to use getLinewidths here instead of diagonal() -- they
    // do not return the same thing (diagonal has population factors)
    VectorBTE phononElectronRates = phelScatteringMatrix.getLinewidths();
    return phononElectronRates;
}*/

void CoupledTransportApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getPhFC2FileName(), "PhFC2FileName");
  throwErrorIfUnset(context.getQMesh(), "qMesh");
  throwWarningIfUnset(context.getSumRuleFC2(), "sumRuleFC2");
  throwErrorIfUnset(context.getPhFC3FileName(), "PhFC3FileName");
  throwErrorIfUnset(context.getTemperatures(), "temperatures");
  throwErrorIfUnset(context.getSmearingMethod(), "smearingMethod");
  if (context.getSmearingMethod() == DeltaFunction::gaussian) {
    throwErrorIfUnset(context.getSmearingWidth(), "smearingWidth");
  }
  if (!context.getElphFileName().empty()) {
    throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");
    throwErrorIfUnset(context.getKMesh(), "kMesh");
    if (context.getDopings().size() == 0 &&
        context.getChemicalPotentials().size() == 0) {
      Error("Either chemical potentials or dopings must be set");
    }
  }
}
