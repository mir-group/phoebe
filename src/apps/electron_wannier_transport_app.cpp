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


void symmetrizeStuff(Context &context, BaseBandStructure &bandStructure,
                     VectorBTE &lifetimes, StatisticsSweep &statisticsSweep) {

  int numCalculations = statisticsSweep.getNumCalculations();

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
  Crystal nofieldCrystal(context, directCell, atomicPositions,
        atomicSpecies, speciesNames, speciesMasses);

  // Points nofieldPoints(nofieldCrystal, context.getKMesh());
  Points nofieldPoints = bfieldPoints;
  nofieldPoints.swapCrystal(nofieldCrystal);

  // Here we want to find the symmetries that rotate both k-points and
  // the group velocities
  {
    std::vector<Eigen::MatrixXd> allVelocities;
    for (int ik = 0; ik < bandStructure.getNumPoints(); ik++) {
      WavevectorIndex ikIdx(ik);
      Eigen::MatrixXd v = bandStructure.getGroupVelocities(ikIdx);
      allVelocities.push_back(v);
    }
    nofieldPoints.setIrreduciblePoints(&allVelocities);
  }
  // nofieldPoints.setIrreduciblePoints();

  // TODO why can we not loop over the irr points from the
  // mesh with the bfield? This mesh will include more irr points
  // than the noField mesh, as some of them will be duplicate under the symmetries
  // of the noField mesh. We could symmetrize those points...
  //
  // I think this would also work in the noField case
  //
  //Strategy -- tell me why it won't work so I understand:
  // could we loop over the irr points on the bfield mesh, and for each one...
  //   * look up the corresponding irr point from the nofield mesh, ik.
  //   * use noFieldPoints.getReducibleStarFromIrreducible(ik) to find
  //     all the other points indices which reduce to this point
  //   * look up these points on the bfield bandstructure, and average them
  //     Use the bfield bandstructure because it's already been computed.
  //
  //   * keep track of which irr points we've visited already, and break the
  //     loop if we try to do the same one twice
  //   * skip points/bands that are removed by active band structure
  //
  //   QUESTION: are the grid point labels for these two grids the same, even
  //   if only for the reducible points? I don't see any reason they should be different
  //   QUESTION: is the bandstructure only saved on the irr mesh, or on the whole reducible mesh?


  // loop over all irreducible points on the noField mesh
  for (int ik : nofieldPoints.irrPointsIterator()) {
    // get coordinates of point ik
    // TODO first we wrote getWavevector here, which is a bands class function.
    // However, because the nofield band structure doesn't exist, can we ask for this from the bfield band structure
    Eigen::Vector3d kCoords = nofieldPoints.getPointCoordinates(ik, Points::cartesianCoordinates);
    // reconstruct the star of wavevectors reducible to ik, by coordinates
    // TODO QUESTION -- what was the reason we didn't use
    // noFieldPoints.getReducibleStarFromIrreducible here? wouldn't that be simpler?
    // maybe not, still need to look up kpoint vectors
    std::vector<Eigen::Vector3d> equivalent_kCoords;
    auto rots = nofieldPoints.getRotationsStar(ik);
    for (const Eigen::Matrix3d &rot : rots) {
      Eigen::Vector3d kRed = rot * kCoords;
      equivalent_kCoords.push_back(kRed);
    }

    // see which wavevectors are in the irreducible set of bfieldPoints
    // save index if the point in bfieldPoints
    std::vector<int> tmp_indices;
    for (auto testVector : equivalent_kCoords) {
      auto testVectorCrystal = bfieldPoints.cartesianToCrystal(testVector);
      int bf_ik = bfieldPoints.isPointStored(testVectorCrystal);
      if (bf_ik == -1) {
        continue;
      }
      tmp_indices.push_back(bf_ik);
    }

    if (tmp_indices.empty()) { continue; }

    // check which indices are irreeducible in bfieldPoints
    std::vector<int> indices_to_be_symmetrized;
    auto irrPointsList = bfieldPoints.irrPointsIterator();
    for (int i : tmp_indices) {
      if (std::find(irrPointsList.begin(), irrPointsList.end(), i) != irrPointsList.end() ) {
        indices_to_be_symmetrized.push_back(i);
      }
    }

    //--------------------
    // symmetrize linewidths
    if ( indices_to_be_symmetrized.size() > 1 ) {
      int nb;
      {
        WavevectorIndex ikIdx(indices_to_be_symmetrized[0]);
        nb = int(bandStructure.getEnergies(ikIdx).size());
      }
      int nk = indices_to_be_symmetrized.size();

      Eigen::Tensor<double,3> avg_linewidths(numCalculations, nk, nb);
      avg_linewidths.setZero();

      for (int ib=0; ib<nb; ++ib) {
        for (int j=0; j<nk; ++j) {
          int ikk = indices_to_be_symmetrized[j];
          WavevectorIndex ikIdx(ikk);
          int is = bandStructure.getIndex(ikIdx, BandIndex(ib));
          StateIndex isIdx(is);
          BteIndex iBteIdx = bandStructure.stateToBte(isIdx);
          int iBte = iBteIdx.get();

          for (int iCalc=0; iCalc<numCalculations; ++iCalc) {
            avg_linewidths(iCalc, j, ib) += lifetimes(iCalc, 0, iBte) / double(nk);
          }
        }
      }
      // copy back
      for (int ib=0; ib<nb; ++ib) {
        for (int j=0; j<nk; ++j) {
          int ikk = indices_to_be_symmetrized[j];
          WavevectorIndex ikIdx(ikk);
          int is = bandStructure.getIndex(ikIdx, BandIndex(ib));
          StateIndex isIdx(is);
          BteIndex iBteIdx = bandStructure.stateToBte(isIdx);
          int iBte = iBteIdx.get();
          for (int iCalc=0; iCalc<numCalculations; ++iCalc) {
            lifetimes(iCalc, 0, iBte) = avg_linewidths(iCalc, j, ib);
          }
        }
      }

      //---------------
      // symmetrize energies and velocities
      Eigen::VectorXd energies(nb);
      Eigen::Tensor<std::complex<double>, 3> velocities(nb, nb, 3);
      energies.setZero();
      velocities.setZero();
      for (int j = 0; j < nk; ++j) {
        int ikk = indices_to_be_symmetrized[j];
        WavevectorIndex ikIdx(ikk);
        for (int ib=0; ib<nb; ++ib) {
          int is = bandStructure.getIndex(ikIdx, BandIndex(ib));
          StateIndex isIdx(is);
          energies(ib) += bandStructure.getEnergy(isIdx) / double(nk);
        }
        auto tmpVel = bandStructure.getVelocities(ikIdx);

        // v is a vector, so it must be rotated before performing the average
        Eigen::Matrix3d rot = nofieldPoints.getRotationFromReducibleIndex(ikk);
        rot = rot.inverse();

        for (int ib1=0; ib1<nb; ++ib1) {
          for (int ib2 = 0; ib2 < nb; ++ib2) {

            Eigen::VectorXcd tmp(3);
            for (int iCart : {0,1,2}) {
              tmp(iCart) = velocities(ib1, ib2, iCart);
            }
            tmp = rot * tmp;
            for (int iCart : {0,1,2}) {
              velocities(ib1, ib2, iCart) += tmp(iCart) / double(nk);
            }
          }
        }
      }

      for (int j = 0; j < nk; ++j) {
        int ikk = indices_to_be_symmetrized[j];
        Point point = bfieldPoints.getPoint(ikk);
        bandStructure.setEnergies(point, energies);

        Eigen::Matrix3d rot = nofieldPoints.getRotationFromReducibleIndex(ikk);

        int nb = velocities.dimensions()[0];
        Eigen::Tensor<std::complex<double>,3> tmpVel(nb,nb,3);
        tmpVel.setZero();
        for (int ib1=0; ib1<nb; ++ib1) {
          for (int ib2 = 0; ib2 < nb; ++ib2) {
            for (int iCart : {0,1,2}) {
              for (int jCart : {0, 1, 2}) {
                tmpVel(ib1,ib2,iCart) += rot(iCart,jCart) * velocities(ib1, ib2, jCart);
              }
            }
          }
        }
        bandStructure.setVelocities(point, tmpVel);
      }
    }
  }
}

void ElectronWannierTransportApp::run(Context &context) {

  auto t2 = Parser::parsePhHarmonic(context);
  auto crystal = std::get<0>(t2);
  auto phononH0 = std::get<1>(t2);

  // if there's a magnetic field, reduce symmetry
  auto bfield = context.getBField();
  if(bfield.norm() != 0) {
    crystal.magneticSymmetries(context);
  }

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
    if(bandStructure.hasWindow() != 0) {
        std::cout << "Window selection reduced electronic band structure from "
                << fullPoints.getNumPoints()*electronH0.getNumBands() << " to "
                << bandStructure.getNumStates() << " states."  << std::endl;
    }
    if(context.getUseSymmetries()) {
      std::cout << "Symmetries reduced electronic band structure from "
          << bandStructure.getNumStates() << " to "
          << bandStructure.irrStateIterator().size() << " states." << std::endl;
    }
    std::cout << "Done computing electronic band structure.\n" << std::endl;
  }

  // Old code for using all the band structure
  //  bool withVelocities = true;
  //  bool withEigenvectors = true;
  //  FullBandStructure bandStructure = electronH0.populate(
  //      fullPoints, withVelocities, withEigenvectors);
  //  // set the chemical potentials to zero, load temperatures
  //  StatisticsSweep statisticsSweep(context, &bandStructure);

  // build/initialize the scattering matrix and the smearing
  ElScatteringMatrix scatteringMatrix(context, statisticsSweep, bandStructure,
                                      bandStructure, phononH0, &couplingElPh);
  scatteringMatrix.setup();

  // TODO there's an internalDiagonal object that has the lifetmes -- we can get this info
  // with .diagonal() or getLinewidths() to get a vector BTE object,
  // but if we want to put them back in the matrix after symmetrizing we will need to
  // add some kind of set function
  auto lifetimes = scatteringMatrix.diagonal();
  symmetrizeStuff(context, bandStructure, lifetimes, statisticsSweep);

  // Add magnetotransport term to scattering matrix if found in input file
  auto magneticField = context.getBField();
  if(magneticField.squaredNorm() != 0) {

    // TODO we might want to throw a warning if someone runs with a
    // bfield that adds a contribution on the order of the diagonal
    // scattering matrix elements, or something like that

    // add -e (vxB) . \nabla f correction to the scattering matrix
    scatteringMatrix.addMagneticTerm(magneticField);
  }

  scatteringMatrix.outputToJSON("rta_el_relaxation_times.json");

  // solve the BTE at the relaxation time approximation level
  // we always do this, as it's the cheapest solver and is required to know
  // the diagonal for the exact method.

  if (mpi->mpiHead()) {
    std::cout << "\n" << std::string(80, '-') << "\n\n";
    std::cout << "Solving BTE within the relaxation time approximation.\n";
  }

  // compute the phonon populations in the relaxation time approximation.
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
    std::cout << "\n" << std::string(80, '-') << "\n" << std::endl;
  }

  if (!std::isnan(context.getConstantRelaxationTime())) {
    return; // if we used the constant RTA, we can't solve the BTE exactly
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
  if (context.getScatteringMatrixInMemory() &&
      statisticsSweep.getNumCalculations() != 1) {
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
      std::cout << std::string(80, '-') << "\n" << std::endl;
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
      std::cout << std::string(80, '-') << "\n" << std::endl;
    }
  }
}

void ElectronWannierTransportApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");
  throwErrorIfUnset(context.getKMesh(), "kMesh");
  throwErrorIfUnset(context.getTemperatures(), "temperatures");

  if (std::isnan(context.getConstantRelaxationTime())) { // non constant tau
    throwErrorIfUnset(context.getElphFileName(), "elphFileName");
    throwErrorIfUnset(context.getSmearingMethod(), "smearingMethod");
    if (context.getSmearingMethod() == DeltaFunction::gaussian) {
      throwErrorIfUnset(context.getSmearingWidth(), "smearingWidth");
    }
  } else {
    if (std::isnan(context.getNumOccupiedStates()) &&
        std::isnan(context.getFermiLevel())) {
      Error("For constant tau calculations, you must provide either the number "
            "of occupied Kohn-Sham states in the valence band or the Fermi "
            "level at T=0K");
    }
  }

  if (context.getDopings().size() == 0 &&
      context.getChemicalPotentials().size() == 0) {
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
  for (size_t iss=0; iss<numParallelIrrStates; iss++) {
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
    Eigen::MatrixXd alphaE = Eigen::MatrixXd::Zero(numCalculations,3);
    Eigen::MatrixXd alphaT = Eigen::MatrixXd::Zero(numCalculations,3);
    {
      Eigen::MatrixXd numE = rE.dot(rE);
      Eigen::MatrixXd numT = rT.dot(rT);
      Eigen::MatrixXd denE = dE.dot(wE);
      Eigen::MatrixXd denT = dT.dot(wT);
      for (int iCalc=0; iCalc<numCalculations; iCalc++) {
        for (int i : {0,1,2}) {
          if (denE(iCalc,i) != 0.) {
            alphaE(iCalc, i) = numE(iCalc, i) / denE(iCalc, i);
          }
          if (denT(iCalc,i) != 0.) {
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
    Eigen::MatrixXd betaE = Eigen::MatrixXd::Zero(numCalculations,3);
    Eigen::MatrixXd betaT = Eigen::MatrixXd::Zero(numCalculations,3);
    {
      Eigen::MatrixXd numE = rNewE.dot(rNewE);
      Eigen::MatrixXd numT = rNewT.dot(rNewT);
      Eigen::MatrixXd denE = rE.dot(rE);
      Eigen::MatrixXd denT = rT.dot(rT);
      for (int iCalc=0; iCalc<numCalculations; iCalc++) {
        for (int i : {0,1,2}) {
          if (denE(iCalc,i) != 0.) {
            betaE(iCalc, i) = numE(iCalc, i) / denE(iCalc, i);
          }
          if (denT(iCalc,i) != 0.) {
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
    for (size_t iss=0; iss<numIrrStates; iss++) {
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
    std::cout << std::string(80, '-') << "\n" << std::endl;
  }
}
