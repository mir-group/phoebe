#include "electron_wannier_transport_app.h"
#include "bandstructure.h"
#include "context.h"
#include "drift.h"
#include "el_scattering.h"
#include "electron_viscosity.h"
#include "exceptions.h"
#include "observable.h"
#include "onsager.h"
#include "parser.h"
#include "wigner_electron.h"
#include "io.h"
#include <nlohmann/json.hpp>

void symmetrizeLinewidths(Context &context, BaseBandStructure &bandStructure,
                          ElScatteringMatrix &scatteringMatrix, StatisticsSweep &statisticsSweep) {
  int numCalculations = statisticsSweep.getNumCalculations();
  auto linewidths = scatteringMatrix.getLinewidths();

  // get points on the irr mesh reduced for mag field symmetries
  auto bfieldPoints = bandStructure.getPoints();
  auto bfieldCrystal = bfieldPoints.getCrystal();

  auto directCell = bfieldCrystal.getDirectUnitCell();
  auto atomicPositions = bfieldCrystal.getAtomicPositions();
  auto atomicSpecies = bfieldCrystal.getAtomicSpecies();
  auto speciesNames = bfieldCrystal.getSpeciesNames();
  auto speciesMasses = bfieldCrystal.getSpeciesMasses();
  // set up a new crystal and points mesh using symmetries
  // without a magnetic field
  Crystal noFieldCrystal(context, directCell, atomicPositions,
                         atomicSpecies, speciesNames, speciesMasses);

  // Points noFieldPoints(noFieldCrystal, context.getKMesh());
  Points noFieldPoints = bfieldPoints;
  noFieldPoints.swapCrystal(noFieldCrystal);

  // Here we want to find the symmetries that rotate both k-points and
  // the group velocities
  {
    std::vector<Eigen::MatrixXd> allVelocities;
    std::vector<Eigen::VectorXd> allEnergies;
    for (int ik = 0; ik < bandStructure.getNumPoints(); ik++) {
      WavevectorIndex ikIdx(ik);
      Eigen::MatrixXd v = bandStructure.getGroupVelocities(ikIdx);
      allVelocities.push_back(v);
      allEnergies.push_back(bandStructure.getEnergies(ikIdx));
    }
    noFieldPoints.setIrreduciblePoints(&allVelocities, &allEnergies);
  }

  // for each irr point, go over all equivalent points and average
  // the linewidths, then reset them
  for (int ikIrr : noFieldPoints.irrPointsIterator()) {

    // get all the kpoints which reduce to this kpoint
    auto reducibleList = noFieldPoints.getReducibleStarFromIrreducible(ikIrr);
    int nkRed = reducibleList.size();
    int numBands;
    {
      WavevectorIndex ikIdx(ikIrr);
      numBands = bandStructure.getNumBands(ikIdx);
    }

    for (int ib = 0; ib < numBands; ++ib) {

      Eigen::VectorXd avgLinewidths(numCalculations);
      avgLinewidths.setZero();
      // TODO: might be possible to switch this to vector and reduce
      //std::vector<double> avgLinewidths(numCalculations);
      //#pragma omp declare reduction(vec_double_plus : std::vector<double> : \
      //      std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
      //      initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

      // loop over all points to be averaged
      //#pragma omp parallel for reduction(vec_double_plus:avgLinewidths)
      for (auto ikRed : reducibleList) {

        // find the bte state index of this (band, kpoint) state
        int is = bandStructure.getIndex(WavevectorIndex(ikRed), BandIndex(ib));
        StateIndex isIdx(is);
        BteIndex iBteIdx = bandStructure.stateToBte(isIdx);
        int iBte = iBteIdx.get();

        for (int iCalc = 0; iCalc < numCalculations; ++iCalc) {
          avgLinewidths(iCalc) += linewidths(iCalc, 0, iBte) / double(nkRed);
        }
      }
      // copy back into all the reducible point linewidths
      for (auto ikRed : reducibleList) {

        // find the bte state index of this (band, kpoint) state
        int is = bandStructure.getIndex(WavevectorIndex(ikRed), BandIndex(ib));
        StateIndex isIdx(is);

        BteIndex iBteIdx = bandStructure.stateToBte(isIdx);
        int iBte = iBteIdx.get();

        for (int iCalc = 0; iCalc < numCalculations; ++iCalc) {
          linewidths(iCalc, 0, iBte) = avgLinewidths(iCalc);
        }
      }
    }
  }
  scatteringMatrix.setLinewidths(linewidths);
}

void unfoldLinewidths(Context &context, ElScatteringMatrix &oldMatrix,
                      ActiveBandStructure &bandStructure,
                      StatisticsSweep &statisticsSweep,
                      HarmonicHamiltonian &electronH0, Points& points,
                      Eigen::Vector3d bfield) {

  for (const std::string &s : context.getSolverBTE()) {
    if (s.compare("iterative") == 0 || s.compare("variational") == 0 || s.compare("relaxons") == 0) {
      Error("Cannot unfold linewidths with non-RTA solver.");
    }
  }

  // unfortunately for indexing, we need a copy of the old and new bandstructures
  ActiveBandStructure oldBandStructure = bandStructure;

  // create a new band structure object with the right symmetries
  points.magneticSymmetries(bfield);
  auto tup = ActiveBandStructure::builder(context, electronH0, points);
  bandStructure = std::get<0>(tup);

  // now the bands should have the right symmetries, and we can attempt the linewidths
  // we take an empty scattering matrix object.
  int numCalculations = statisticsSweep.getNumCalculations();

  VectorBTE newLinewidths(statisticsSweep, bandStructure, 1);
  VectorBTE oldLinewidths = oldMatrix.getLinewidths();

  // for each irr point of the new band structure, we will get the
  // wavevector of that point, look it up in the old list, and
  // use that index to set the new linewidth
  //
  // this is more straightforward because at the end, we only care to
  // save new irr points anyway
  for (int ikIrr : bandStructure.getPoints().irrPointsIterator()) {

    // this should work because all sym eq kpoints should have the same nbands
    // as a result of enforce function in active band structure
    int numBands;
    {
      WavevectorIndex ikIdx(ikIrr);
      numBands = oldBandStructure.getNumBands(ikIdx);
    }

    // get the actual wavevector and look it up manually
    WavevectorIndex temp = WavevectorIndex(ikIrr);
    auto kCoords = bandStructure.getWavevector(temp);
    kCoords = bandStructure.getPoints().cartesianToCrystal(kCoords);
    int ikOld = oldBandStructure.getPointIndex(kCoords);

    // for each band, we replace all the equivalent kpoints
    for (int ib = 0; ib < numBands; ++ib) {

      int iBteOld;
      {
        int is = oldBandStructure.getIndex(WavevectorIndex(ikOld), BandIndex(ib));
        StateIndex isIdx(is);
        BteIndex iBteIdx = oldBandStructure.stateToBte(isIdx);
        iBteOld = iBteIdx.get();
      }

      // find the bte state index of this (band, kpoint) state
      // must convert the reducible index (which should be the same for both band structures)
      // to the irr representation of the new band structure
      int is = bandStructure.getIndex(WavevectorIndex(ikIrr), BandIndex(ib));
      StateIndex isIdx(is);
      // check that indices map to exactly the same energies
      BteIndex iBteIdx = bandStructure.stateToBte(isIdx);
      int iBte = iBteIdx.get();

      for (int iCalc = 0; iCalc < numCalculations; ++iCalc) {
        newLinewidths(iCalc, 0, iBte) = oldLinewidths(iCalc, 0, iBteOld);
      }
    }
  }
  bool supressError = true; // normally, it's bad to set the matrix
                            // with different # of states, so we thrown an
                            // error. Here we have a special case.
  oldMatrix.setLinewidths(newLinewidths, supressError);
}

void ElectronWannierTransportApp::run(Context &context) {

  auto t2 = Parser::parsePhHarmonic(context);
  auto crystal = std::get<0>(t2);
  auto phononH0 = std::get<1>(t2);

  // if there's a magnetic field, reduce symmetry.
  // we intentionally do this after the electronic structure
  //auto bfield = context.getBField();
  //if(bfield.norm() != 0) {
  //  crystal.magneticSymmetries(context);
  //}

  auto t1 = Parser::parseElHarmonicWannier(context, &crystal);
  auto crystalEl = std::get<0>(t1);
  auto electronH0 = std::get<1>(t1);

  // load the elph coupling
  // Note: this file contains the number of electrons
  // which is needed to understand where to place the fermi level
  InteractionElPhWan couplingElPh(crystal);
  if (std::isnan(context.getConstantRelaxationTime())) {
    couplingElPh = InteractionElPhWan::parse(context, crystal, &phononH0);
  }

  // compute the band structure on the fine grid
  if (mpi->mpiHead()) {
    std::cout << "\nComputing electronic band structure." << std::endl;
  }

  Points fullPoints(crystal, context.getKMesh());
  auto t3 = ActiveBandStructure::builder(context, electronH0, fullPoints);
  auto bandStructure = std::get<0>(t3);
  auto statisticsSweep = std::get<1>(t3);

  // print some info about how window and symmetries have reduced things
  if (mpi->mpiHead()) {
    if (bandStructure.hasWindow() != 0) {
      std::cout << "Window selection reduced electronic band structure from "
                << fullPoints.getNumPoints() * electronH0.getNumBands() << " to "
                << bandStructure.getNumStates() << " states." << std::endl;
    }
    if (context.getUseSymmetries()) {
      std::cout << "Symmetries reduced electronic band structure from "
                << bandStructure.getNumStates() << " to "
                << bandStructure.irrStateIterator().size() << " states." << std::endl;
    }
    std::cout << "Done computing electronic band structure.\n"
              << std::endl;
  }

  // build/initialize the scattering matrix and the smearing
  ElScatteringMatrix scatteringMatrix(context, statisticsSweep, bandStructure,
                                      bandStructure, phononH0, &couplingElPh);
  scatteringMatrix.setup();

  // symmetrize the linewidths on the diagonal
  //symmetrizeLinewidths(context, bandStructure, scatteringMatrix, statisticsSweep);

  // TODO we might want to throw a warning if someone runs with a
  // bfield that adds a contribution on the order of the diagonal
  // scattering matrix elements, or something like that -- if we have effective mass, we could also estimate omega_c for comparison

  // Add magnetotransport term to scattering matrix if found in input file
  std::vector<Eigen::Vector3d> magneticFields = context.getBField();

// TODO MPI or OMP Parallel this
  // if there are multiple bfields, we run everything for each of them
  for(int b = 0; b < int(magneticFields.size()); b++) {

    Eigen::Vector3d magneticField = magneticFields[b];

    // crystal inside points is a reference, so we need to be careful
    // and make an all new copy
    auto directCell = crystal.getDirectUnitCell();
    auto atomicPositions = crystal.getAtomicPositions();
    auto atomicSpecies = crystal.getAtomicSpecies();
    auto speciesNames = crystal.getSpeciesNames();
    auto speciesMasses = crystal.getSpeciesMasses();
    // set up a new crystal and points mesh using symmetries
    Crystal tempCrystal(context, directCell, atomicPositions,
                           atomicSpecies, speciesNames, speciesMasses);
    Points tempPoints(tempCrystal, context.getKMesh());

    // make a new copy of the scattering matrix
    auto tup = ActiveBandStructure::builder(context, electronH0, tempPoints);
    ActiveBandStructure tempBS = std::get<0>(tup);
    auto tempStatisticsSweep = std::get<1>(tup);

    // make a new scattering matrix with the new bandstructure
    ElScatteringMatrix tempScatteringMatrix(scatteringMatrix, tempBS, tempBS);

    // TODO only do this if it's RTA. Also, turn on mag syms earlier if it's iterative, relaxons, etc
    // unfold the symmetries
    if(context.getUseSymmetries()) {
      unfoldLinewidths(context, tempScatteringMatrix, tempBS, statisticsSweep, electronH0, tempPoints, magneticField);
    }

    // add -e (vxB) . \nabla f correction to the scattering matrix
    tempScatteringMatrix.addMagneticTerm(magneticField);

    if (mpi->mpiHead()) {
      std::cout << "Solving BTE within the relaxation time approximation "
            "for B = " << magneticField.transpose() * 1./teslaToAu << "." << std::endl;
    }

    // compute the populations in the relaxation time approximation.
    BulkEDrift driftE(tempStatisticsSweep, tempBS, 3);
    BulkTDrift driftT(tempStatisticsSweep, tempBS, 3);
    VectorBTE relaxationTimes = tempScatteringMatrix.getSingleModeTimes();
    VectorBTE nERTA = -driftE * relaxationTimes;
    VectorBTE nTRTA = -driftT * relaxationTimes;

    // compute the electrical conductivity -- TODO check crystal syms and make sure it's not wrongly syming?
    OnsagerCoefficients transportCoefficients(tempStatisticsSweep, tempCrystal, tempBS, context);
    transportCoefficients.calcFromPopulation(nERTA, nTRTA);
    transportCoefficients.print();
    transportCoefficients.outputToJSON("B" + std::to_string(b) + "_rta_onsager_coefficients.json");

    if(mpi->mpiHead()) {
      std::ifstream i("B" +std::to_string(b) + "_rta_onsager_coefficients.json");
      nlohmann::json j; i >> j;
      std::vector<double> vec(3);
      vec[0] = magneticField(0)/teslaToAu; vec[1] = magneticField(1)/teslaToAu; vec[2] = magneticField(2)/teslaToAu;
      j["magneticField"] = vec;
      j["magneticFieldUnit"] = "Tesla";
      std::ofstream o("B" +std::to_string(b) + "_rta_onsager_coefficients.json");
      o << std::setw(5) << j << std::endl;
      o.close();
    }
   tempScatteringMatrix.outputToJSON("B" +std::to_string(b) + "_rta_el_relaxation_times.json");
  }

  scatteringMatrix.outputToJSON("rta_el_relaxation_times.json");

  // solve the BTE at the relaxation time approximation level
  // we always do this, as it's the cheapest solver and is required to know
  // the diagonal for the exact method.

  if (mpi->mpiHead()) {
    std::cout << "\n"
              << std::string(80, '-') << "\n\n";
    std::cout << "Solving BTE within the relaxation time approximation.\n";
  }

  // compute the populations in the relaxation time approximation.
  // TODO fix this for electrons --
  // Note: this is the total phonon population n (n != f(1+f) Delta n)

  BulkEDrift driftE(statisticsSweep, bandStructure, 3);
  BulkTDrift driftT(statisticsSweep, bandStructure, 3);
  VectorBTE relaxationTimes = scatteringMatrix.getSingleModeTimes();
  VectorBTE nERTA = -driftE * relaxationTimes;
  VectorBTE nTRTA = -driftT * relaxationTimes;

  // compute the electrical conductivity
  OnsagerCoefficients transportCoefficients(statisticsSweep, crystal,
                                            bandStructure, context);

  transportCoefficients.calcFromPopulation(nERTA, nTRTA);
  transportCoefficients.print();
  transportCoefficients.outputToJSON("rta_onsager_coefficients.json");

  // compute the Wigner transport coefficients
  WignerElCoefficients wignerCoefficients(
      statisticsSweep, crystal, bandStructure, context, relaxationTimes);
  wignerCoefficients.calcFromPopulation(nERTA, nTRTA);
  wignerCoefficients.print();
  wignerCoefficients.outputToJSON("rta_wigner_coefficients.json");

  // compute the electron viscosity
  ElectronViscosity elViscosity(context, statisticsSweep, crystal,
                                bandStructure);
  elViscosity.calcRTA(relaxationTimes);
  elViscosity.print();
  elViscosity.outputToJSON("rta_electron_viscosity.json");

  // compute the specific heat
  SpecificHeat specificHeat(context, statisticsSweep, crystal, bandStructure);
  specificHeat.calc();
  specificHeat.print();
  specificHeat.outputToJSON("el_specific_heat.json");

  if (mpi->mpiHead()) {
    std::cout << "\n"
              << std::string(80, '-') << "\n"
              << std::endl;
  }

  if (!std::isnan(context.getConstantRelaxationTime())) {
    return;// if we used the constant RTA, we can't solve the BTE exactly
  }

  //---------------------------------------------------------------------------
  // start section on exact BTE solvers

  std::vector<std::string> solverBTE = context.getSolverBTE();

  bool doIterative = false;
  bool doVariational = false;
  bool doRelaxons = false;
  for (const std::string &s : solverBTE) {
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
  if (context.getScatteringMatrixInMemory() && statisticsSweep.getNumCalculations() != 1) {
    Error("If scattering matrix is kept in memory, only one "
          "temperature/chemical potential is allowed in a run");
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

  if (context.getScatteringMatrixInMemory() && !context.getUseSymmetries()) {
    if (doVariational || doRelaxons || doIterative) {
      if (context.getSymmetrizeMatrix()) {
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
      std::cout << "Starting Omini Sparavigna BTE solver\n"
                << std::endl;
    }

    Eigen::Tensor<double, 3> elCond =
        transportCoefficients.getElectricalConductivity();
    Eigen::Tensor<double, 3> thCond =
        transportCoefficients.getThermalConductivity();
    Eigen::Tensor<double, 3> elCondOld = elCond;
    Eigen::Tensor<double, 3> thCondOld = thCond;

    VectorBTE nENext(statisticsSweep, bandStructure, 3);
    VectorBTE nTNext(statisticsSweep, bandStructure, 3);

    VectorBTE lineWidths = scatteringMatrix.getLinewidths();
    VectorBTE nEOld = nERTA;
    VectorBTE nTOld = nTRTA;

    auto threshold = context.getConvergenceThresholdBTE();

    for (int iter = 0; iter < context.getMaxIterationsBTE(); iter++) {

      std::vector<VectorBTE> nIn;
      nIn.push_back(nEOld);
      nIn.push_back(nTOld);
      auto nOut = scatteringMatrix.offDiagonalDot(nIn);
      nENext = nOut[0] / lineWidths;
      nTNext = nOut[1] / lineWidths;
      nENext = nERTA - nENext;
      nTNext = nTRTA - nTNext;

      transportCoefficients.calcFromPopulation(nENext, nTNext);
      transportCoefficients.print(iter);
      elCond = transportCoefficients.getElectricalConductivity();
      thCond = transportCoefficients.getThermalConductivity();

      // this exit condition must be improved
      // different temperatures might converge differently

      double dE = findMaxRelativeDifference(elCond, elCondOld);
      double dT = findMaxRelativeDifference(thCond, thCondOld);
      if ((dE < threshold) && (dT < threshold)) {
        break;
      } else {
        elCondOld = elCond;
        thCondOld = thCond;
        nEOld = nENext;
        nTOld = nTNext;
      }

      if (iter == context.getMaxIterationsBTE() - 1) {
        Error("Reached max BTE iterations without convergence");
      }
    }
    transportCoefficients.print();
    transportCoefficients.outputToJSON("omini_onsager_coefficients.json");

    if (mpi->mpiHead()) {
      std::cout << "Finished Omini-Sparavigna BTE solver\n\n";
      std::cout << std::string(80, '-') << "\n"
                << std::endl;
    }
  }

  if (doVariational) {
    runVariationalMethod(context, crystal, statisticsSweep, bandStructure,
                         scatteringMatrix);
  }

  if (doRelaxons) {
    if (mpi->mpiHead()) {
      std::cout << "Starting relaxons BTE solver" << std::endl;
    }
    // scatteringMatrix.a2Omega(); Important!! Must use the matrix A, non-scaled
    // this because the scaling used for phonons here would cause instability
    // as it could contains factors like 1/0
    auto tup = scatteringMatrix.diagonalize();
    Eigen::VectorXd eigenvalues = std::get<0>(tup);
    ParallelMatrix<double> eigenvectors = std::get<1>(tup);
    // EV such that Omega = V D V^-1

    transportCoefficients.calcFromRelaxons(eigenvalues, eigenvectors,
                                           scatteringMatrix);
    transportCoefficients.print();
    transportCoefficients.outputToJSON("relaxons_onsager_coefficients.json");
    scatteringMatrix.relaxonsToJSON("exact_relaxation_times.json", eigenvalues);

    if (!context.getUseSymmetries()) {
      Vector0 fermiEigenvector(statisticsSweep, bandStructure, specificHeat);
      elViscosity.calcFromRelaxons(eigenvalues, eigenvectors);
      elViscosity.print();
      elViscosity.outputToJSON("relaxons_electron_viscosity.json");
    }

    if (mpi->mpiHead()) {
      std::cout << "Finished relaxons BTE solver\n\n";
      std::cout << std::string(80, '-') << "\n"
                << std::endl;
    }
  }
}

void ElectronWannierTransportApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");
  throwErrorIfUnset(context.getKMesh(), "kMesh");
  throwErrorIfUnset(context.getTemperatures(), "temperatures");

  if (std::isnan(context.getConstantRelaxationTime())) {// non constant tau
    throwErrorIfUnset(context.getElphFileName(), "elphFileName");
    throwErrorIfUnset(context.getSmearingMethod(), "smearingMethod");
    if (context.getSmearingMethod() == DeltaFunction::gaussian) {
      throwErrorIfUnset(context.getSmearingWidth(), "smearingWidth");
    }
  } else {
    if (std::isnan(context.getNumOccupiedStates()) && std::isnan(context.getFermiLevel())) {
      Error("For constant tau calculations, you must provide either the number "
            "of occupied Kohn-Sham states in the valence band or the Fermi "
            "level at T=0K");
    }
  }

  if (context.getDopings().size() == 0 && context.getChemicalPotentials().size() == 0) {
    Error("Either chemical potentials or dopings must be set");
  }
}

void ElectronWannierTransportApp::runVariationalMethod(
    Context &context, Crystal &crystal, StatisticsSweep &statisticsSweep,
    ActiveBandStructure &bandStructure, ElScatteringMatrix &scatteringMatrix) {

  // here we implement a conjugate gradient solution to Az = b

  if (mpi->mpiHead()) {
    std::cout << "Starting variational BTE solver\n";
    std::cout << std::endl;
  }

  OnsagerCoefficients transportCoefficients(statisticsSweep, crystal,
                                            bandStructure, context);

  // initialize the transport coefficients
  // to check how these change during iterations
  Eigen::Tensor<double, 3> elCond =
      transportCoefficients.getElectricalConductivity();
  Eigen::Tensor<double, 3> thCond =
      transportCoefficients.getThermalConductivity();
  auto elCondOld = elCond.setConstant(1.);
  auto thCondOld = thCond.setConstant(1.);

  int numCalculations = statisticsSweep.getNumCalculations();

  // initialize b
  // Note: we solve Az=b not with z as electron population, but as canonical
  // population, i.e. with z being f = f0 + f0 (1-f0) z
  // where f0 is the Fermi-Dirac distribution
  // so, here we compute b without the factor n(1-n)
  // note also, we cannot do divisions by n(1-n) because it could be zero
  VectorBTE bE(statisticsSweep, bandStructure, 3);
  VectorBTE bT(statisticsSweep, bandStructure, 3);
  VectorBTE preconditioning(statisticsSweep, bandStructure, 3);
  Particle particle = bandStructure.getParticle();
  auto parallelIrrStates = bandStructure.parallelIrrStateIterator();
  size_t numParallelIrrStates = parallelIrrStates.size();
#pragma omp parallel for
  for (size_t iss = 0; iss < numParallelIrrStates; iss++) {
    int is = parallelIrrStates[iss];
    StateIndex isIdx(is);
    double energy = bandStructure.getEnergy(isIdx);
    Eigen::Vector3d vel = bandStructure.getGroupVelocity(isIdx);
    int iBte = bandStructure.stateToBte(isIdx).get();
    for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      auto chemPot = calcStat.chemicalPotential;
      auto temp = calcStat.temperature;
      for (int i : {0, 1, 2}) {
        bE(iCalc, i, iBte) = vel(i) / temp;
        bT(iCalc, i, iBte) = -vel(i) * (energy - chemPot) / temp / temp;
        preconditioning(iCalc, i, iBte) = 1.;
      }
    }
  }
  mpi->allReduceSum(&bE.data);
  mpi->allReduceSum(&bT.data);
  mpi->allReduceSum(&preconditioning.data);

  // populations, initial guess
  VectorBTE zNewE = bE;
  VectorBTE zNewT = bT;
  VectorBTE zE = zNewE;
  VectorBTE zT = zNewT;

  std::vector<VectorBTE> inVec;
  inVec.push_back(zE);
  inVec.push_back(zT);
  std::vector<VectorBTE> w0s = scatteringMatrix.dot(inVec);
  VectorBTE w0E = w0s[0];
  VectorBTE w0T = w0s[1];

  // residual
  VectorBTE rE = bE - w0E;
  VectorBTE rT = bT - w0T;

  // search direction
  VectorBTE dE = rE;
  VectorBTE dT = rT;

  double threshold = context.getConvergenceThresholdBTE();

  for (int iter = 0; iter < context.getMaxIterationsBTE(); iter++) {
    // execute CG step, as in
    // https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf

    // compute w = A d
    std::vector<VectorBTE> inW;
    inW.push_back(dE);
    inW.push_back(dT);
    std::vector<VectorBTE> outW = scatteringMatrix.dot(inW);
    VectorBTE wE = outW[0];
    VectorBTE wT = outW[1];

    // amount of descent along the search direction
    //    Eigen::MatrixXd alphaE = (rE.dot(rE)).array() / (dE.dot(wE)).array();
    //    Eigen::MatrixXd alphaT = (rT.dot(rT)).array() / (dT.dot(wT)).array();
    int numCalculations = statisticsSweep.getNumCalculations();
    Eigen::MatrixXd alphaE = Eigen::MatrixXd::Zero(numCalculations, 3);
    Eigen::MatrixXd alphaT = Eigen::MatrixXd::Zero(numCalculations, 3);
    {
      Eigen::MatrixXd numE = rE.dot(rE);
      Eigen::MatrixXd numT = rT.dot(rT);
      Eigen::MatrixXd denE = dE.dot(wE);
      Eigen::MatrixXd denT = dT.dot(wT);
      for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
        for (int i : {0, 1, 2}) {
          if (denE(iCalc, i) != 0.) {
            alphaE(iCalc, i) = numE(iCalc, i) / denE(iCalc, i);
          }
          if (denT(iCalc, i) != 0.) {
            alphaT(iCalc, i) = numT(iCalc, i) / denT(iCalc, i);
          }
        }
      }
    }

    // new guess of population
    zNewE = dE * alphaE + zE;
    zNewT = dT * alphaT + zT;

    // new estimate of residual
    VectorBTE tmpE = wE * alphaE;
    VectorBTE tmpT = wT * alphaT;
    VectorBTE rNewE = rE - tmpE;
    VectorBTE rNewT = rT - tmpT;

    // amount of correction for the search direction
    //    Eigen::MatrixXd betaE = (rNewE.dot(rNewE)).array() / (rE.dot(rE)).array();
    //    Eigen::MatrixXd betaT = (rNewT.dot(rNewT)).array() / (rT.dot(rT)).array();
    Eigen::MatrixXd betaE = Eigen::MatrixXd::Zero(numCalculations, 3);
    Eigen::MatrixXd betaT = Eigen::MatrixXd::Zero(numCalculations, 3);
    {
      Eigen::MatrixXd numE = rNewE.dot(rNewE);
      Eigen::MatrixXd numT = rNewT.dot(rNewT);
      Eigen::MatrixXd denE = rE.dot(rE);
      Eigen::MatrixXd denT = rT.dot(rT);
      for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
        for (int i : {0, 1, 2}) {
          if (denE(iCalc, i) != 0.) {
            betaE(iCalc, i) = numE(iCalc, i) / denE(iCalc, i);
          }
          if (denT(iCalc, i) != 0.) {
            betaT(iCalc, i) = numT(iCalc, i) / denT(iCalc, i);
          }
        }
      }
    }

    // new search direction
    VectorBTE dNewE = dE * betaE + rNewE;
    VectorBTE dNewT = dT * betaT + rNewT;

    // now we update the guess of transport coefficients
    // first though, we add the factors n(1-n) that we removed from the CG
    // we need to do this because we didn't remove this factor from sMatrix
    VectorBTE z2E = zNewE;
    VectorBTE z2T = zNewT;
    auto irrStates = bandStructure.irrStateIterator();
    size_t numIrrStates = irrStates.size();
#pragma omp parallel for
    for (size_t iss = 0; iss < numIrrStates; iss++) {
      int is = irrStates[iss];
      StateIndex isIdx(is);
      double energy = bandStructure.getEnergy(isIdx);
      int iBte = bandStructure.stateToBte(isIdx).get();
      for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
        auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
        auto chemPot = calcStat.chemicalPotential;
        auto temp = calcStat.temperature;
        double pop = sqrt(particle.getPopPopPm1(energy, temp, chemPot));
        for (int i : {0, 1, 2}) {
          z2E(iCalc, i, iBte) *= pop;
          z2T(iCalc, i, iBte) *= pop;
        }
      }
    }
    std::vector<VectorBTE> inF;
    inF.push_back(z2E);
    inF.push_back(z2T);
    auto outF = scatteringMatrix.dot(inF);
    VectorBTE az2E = outF[0];
    VectorBTE az2T = outF[1];
    transportCoefficients.calcVariational(az2E, az2T, z2E, z2T, zNewE, zNewT,
                                          preconditioning);
    transportCoefficients.print(iter);
    elCond = transportCoefficients.getElectricalConductivity();
    thCond = transportCoefficients.getThermalConductivity();

    // decide whether to exit or run the next iteration
    double deltaE = findMaxRelativeDifference(elCond, elCondOld);
    double deltaT = findMaxRelativeDifference(thCond, thCondOld);
    if ((deltaE < threshold) && (deltaT < threshold)) {
      // recompute our final guess for transport coefficients
      transportCoefficients.calcFromCanonicalPopulation(zNewE, zNewT);
      break;
    } else {
      elCondOld = elCond;
      thCondOld = thCond;
      zE = zNewE;
      zT = zNewT;
      rE = rNewE;
      rT = rNewT;
      dE = dNewE;
      dT = dNewT;
    }

    if (iter == context.getMaxIterationsBTE() - 1) {
      Error("Reached max BTE iterations without convergence");
    }
  }

  // nice formatting of the transport properties at the last step
  transportCoefficients.print();
  transportCoefficients.outputToJSON("variational_onsager_coefficients.json");

  if (mpi->mpiHead()) {
    std::cout << "Finished variational BTE solver\n\n";
    std::cout << std::string(80, '-') << "\n"
              << std::endl;
  }
}
