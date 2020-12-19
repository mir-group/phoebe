#include "electron_wannier_transport_app.h"
#include "bandstructure.h"
#include "context.h"
#include "drift.h"
#include "el_scattering.h"
#include "electron_viscosity.h"
#include "exceptions.h"
#include "io.h"
#include "observable.h"
#include "onsager.h"
#include "qe_input_parser.h"
#include "wigner_electron.h"

void ElectronWannierTransportApp::run(Context &context) {

  auto t2 = QEParser::parsePhHarmonic(context);
  auto crystal = std::get<0>(t2);
  auto phononH0 = std::get<1>(t2);

  auto t1 = QEParser::parseElHarmonicWannier(context, &crystal);
  auto crystalEl = std::get<0>(t1);
  auto electronH0 = std::get<1>(t1);

  // load the 3phonon coupling
  // Note: this file contains the number of electrons
  // which is needed to understand where to place the fermi level
  auto couplingElPh = InteractionElPhWan::parse(context, crystal, &phononH0);

  // compute the band structure on the fine grid
  FullPoints fullPoints(crystal, context.getKMesh());
  auto t3 = ActiveBandStructure::builder(context, electronH0, fullPoints);
  auto bandStructure = std::get<0>(t3);
  auto statisticsSweep = std::get<1>(t3);

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

  auto dimensionality = context.getDimensionality();
  BulkEDrift driftE(statisticsSweep, bandStructure, dimensionality);
  BulkTDrift driftT(statisticsSweep, bandStructure, dimensionality);
  VectorBTE relaxationTimes = scatteringMatrix.getSingleModeTimes();
  VectorBTE nERTA = driftE * relaxationTimes;
  VectorBTE nTRTA = driftT * relaxationTimes;

  // compute the electrical conductivity
  OnsagerCoefficients transportCoeffs(statisticsSweep, crystal, bandStructure,
                                      context);
  transportCoeffs.calcFromPopulation(nERTA, nTRTA);
  transportCoeffs.print();
  transportCoeffs.outputToJSON("rta_onsager_coefficients.json");

  // compute the Wigner transport coefficients
  WignerElCoefficients wignerCoeffs(statisticsSweep, crystal, bandStructure,
                                    context, relaxationTimes);
  wignerCoeffs.calcFromPopulation(nERTA, nTRTA);
  wignerCoeffs.print();
  wignerCoeffs.outputToJSON("rta_wigner_coefficients.json");

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

  //---------------------------------------------------------------------------

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
    Error e("Relaxons require matrix kept in memory");
  }
  if (context.getScatteringMatrixInMemory() &&
      statisticsSweep.getNumCalcs() != 1) {
    Error e("If scattering matrix is kept in memory, only one "
            "temperature/chemical potential is allowed in a run");
  }

  mpi->barrier();

  if (doIterative) {

    if (mpi->mpiHead()) {
      std::cout << "Starting Omini Sparavigna BTE solver\n" << std::endl;
    }

    // initialize the (old) transport coefficients
    OnsagerCoefficients transportCoeffsOld = transportCoeffs;

    Eigen::Tensor<double, 3> elCond =
        transportCoeffs.getElectricalConductivity();
    Eigen::Tensor<double, 3> thCond = transportCoeffs.getThermalConductivity();
    Eigen::Tensor<double, 3> elCondOld = elCond;
    Eigen::Tensor<double, 3> thCondOld = thCond;

    VectorBTE fENext(statisticsSweep, bandStructure, dimensionality);
    VectorBTE fTNext(statisticsSweep, bandStructure, dimensionality);
    VectorBTE sMatrixDiagonal = scatteringMatrix.diagonal();

    // from n, we get f, such that n = bose(bose+1)f
    VectorBTE lineWidths = scatteringMatrix.getLinewidths();
    VectorBTE fERTA = driftE / lineWidths;
    VectorBTE fTRTA = driftT / lineWidths;
    VectorBTE fEOld = fERTA;
    VectorBTE fTOld = fTRTA;

    auto threshold = context.getConvergenceThresholdBTE();

    for (int iter = 0; iter < context.getMaxIterationsBTE(); iter++) {

      std::vector<VectorBTE> fIn;
      fIn.push_back(fEOld);
      fIn.push_back(fTOld);
      auto fOut = scatteringMatrix.offDiagonalDot(fIn);
      fENext = fOut[0] / sMatrixDiagonal;
      fTNext = fOut[1] / sMatrixDiagonal;
      fENext = fERTA - fENext;
      fTNext = fTRTA - fTNext;

      transportCoeffs.calcFromCanonicalPopulation(fENext, fTNext);
      transportCoeffs.print(iter);
      elCond = transportCoeffs.getElectricalConductivity();
      thCond = transportCoeffs.getThermalConductivity();

      // this exit condition must be improved
      // different temperatures might converge differently

      Eigen::Tensor<double, 3> diffE = ((elCond - elCondOld) / elCondOld).abs();
      Eigen::Tensor<double, 3> diffT = ((thCond - thCondOld) / thCondOld).abs();
      double dE = 10000.;
      double dT = 10000.;
      for (int i = 0; i < diffE.dimension(0); i++) {
        for (int j = 0; j < diffE.dimension(0); j++) {
          for (int k = 0; k < diffE.dimension(0); k++) {
            if (diffE(i, j, k) < dE)
              dE = diffE(i, j, k);
            if (diffT(i, j, k) < dT)
              dT = diffT(i, j, k);
          }
        }
      }
      if ((dE < threshold) && (dT < threshold)) {
        break;
      } else {
        elCondOld = elCond;
        thCondOld = thCond;
        fEOld = fENext;
        fTOld = fTNext;
      }

      if (iter == context.getMaxIterationsBTE() - 1) {
        Error e("Reached max BTE iterations without convergence");
      }
    }
    transportCoeffs.print();
    transportCoeffs.outputToJSON("omini_onsager_coefficients.json");

    if (mpi->mpiHead()) {
      std::cout << "Finished Omini-Sparavigna BTE solver\n\n";
      std::cout << std::string(80, '-') << "\n" << std::endl;
    }
  }

  if (doVariational) {
    if (mpi->mpiHead()) {
      std::cout << "Starting variational BTE solver\n";
      std::cout << std::endl;
    }

    // initialize the (old) transport coefficients
    OnsagerCoefficients transportCoeffsOld = transportCoeffs;
    Eigen::Tensor<double, 3> elCond =
        transportCoeffs.getElectricalConductivity();
    Eigen::Tensor<double, 3> thCond = transportCoeffs.getThermalConductivity();
    auto elCondOld = elCond;
    auto thCondOld = thCond;

    // set the initial guess to the RTA solution
    VectorBTE lineWidths = scatteringMatrix.getLinewidths();
    VectorBTE fENew = driftE / lineWidths;
    VectorBTE fTNew = driftT / lineWidths;

    VectorBTE preconditioning = lineWidths.sqrt();

    // CG rescaling
    fENew = fENew * preconditioning;
    fTNew = fTNew * preconditioning;

    // save the population of the previous step
    VectorBTE fEOld = fENew;
    VectorBTE fTOld = fTNew;

    // do the conjugate gradient method for transport coefficients
    std::vector<VectorBTE> fIn;
    fIn.push_back(fENew);
    fIn.push_back(fTNew);
    auto gOld = scatteringMatrix.offDiagonalDot(fIn);
    VectorBTE gEOld = gOld[0] / lineWidths; // CG scaling
    VectorBTE gTOld = gOld[1] / lineWidths; // CG scaling
    gEOld = gEOld - fEOld;
    gTOld = gTOld - fTOld;
    VectorBTE hEOld = -gEOld;
    VectorBTE hTOld = -gTOld;

    std::vector<VectorBTE> tIn;
    tIn.push_back(hEOld);
    tIn.push_back(hTOld);
    auto tOut = scatteringMatrix.dot(tIn);
    VectorBTE tEOld = tOut[0];
    VectorBTE tTOld = tOut[1];
    tEOld = tEOld / lineWidths; // CG scaling
    tTOld = tTOld / lineWidths; // CG scaling

    double threshold = context.getConvergenceThresholdBTE();

    for (int iter = 0; iter < context.getMaxIterationsBTE(); iter++) {
      // execute CG step, as in

      Eigen::MatrixXd alphaE =
          (gEOld.dot(hEOld)).array() / (hEOld.dot(tEOld)).array();
      Eigen::MatrixXd alphaT =
          (gTOld.dot(hTOld)).array() / (hTOld.dot(tTOld)).array();

      fENew = hEOld * alphaE;
      fTNew = hTOld * alphaT;
      fENew = fEOld - fENew;
      fTNew = fTOld - fTNew;

      VectorBTE gENew = tEOld * alphaE;
      VectorBTE gTNew = tTOld * alphaT;
      gENew = gEOld - gENew;
      gTNew = gTOld - gTNew;

      Eigen::MatrixXd betaE =
          (gENew.dot(gENew)).array() / (gEOld.dot(gEOld)).array();
      Eigen::MatrixXd betaT =
          (gTNew.dot(gTNew)).array() / (gTOld.dot(gTOld)).array();
      VectorBTE hENew = hEOld * betaE;
      VectorBTE hTNew = hTOld * betaT;
      hENew = -gENew + hENew;
      hTNew = -gTNew + hTNew;

      std::vector<VectorBTE> inVectors;
      inVectors.push_back(fENew);
      inVectors.push_back(hENew); // note: at next step hNew is hOld -> gives tOld
      inVectors.push_back(fTNew);
      inVectors.push_back(hTNew); // note: at next step hNew is hOld -> gives tOld
      auto outVectors = scatteringMatrix.dot(inVectors);
      tEOld = outVectors[1];
      tTOld = outVectors[3];
      tEOld = tEOld / lineWidths; // CG scaling
      tTOld = tTOld / lineWidths; // CG scaling

      transportCoeffs.calcVariational(outVectors[0], outVectors[2], fENew,
                                      fTNew, preconditioning);
      transportCoeffs.print(iter);
      elCond = transportCoeffs.getElectricalConductivity();
      thCond = transportCoeffs.getThermalConductivity();

      // decide whether to exit or run the next iteration
      Eigen::Tensor<double, 3> diffE = ((elCond - elCondOld) / elCondOld).abs();
      Eigen::Tensor<double, 3> diffT = ((thCond - thCondOld) / thCondOld).abs();
      double dE = 10000.;
      double dT = 10000.;
      for (int i = 0; i < diffE.dimension(0); i++) {
        for (int j = 0; j < diffE.dimension(0); j++) {
          for (int k = 0; k < diffE.dimension(0); k++) {
            if (diffE(i, j, k) < dE)
              dE = diffE(i, j, k);
            if (diffT(i, j, k) < dT)
              dT = diffT(i, j, k);
          }
        }
      }
      if ((dE < threshold) && (dT < threshold)) {
        // this because calcVariational computes LTT and LEE
        fENew = fENew / preconditioning;
        fTNew = fTNew / preconditioning;
        transportCoeffs.calcFromCanonicalPopulation(fENew, fTNew);
        break;
      } else {
        elCondOld = elCond;
        thCondOld = thCond;
        fEOld = fENew;
        fTOld = fTNew;
        gEOld = gENew;
        gTOld = gTNew;
        hEOld = hENew;
        hTOld = hTNew;
      }

      if (iter == context.getMaxIterationsBTE() - 1) {
        Error e("Reached max BTE iterations without convergence");
      }
    }

    // nice formatting of the transport properties at the last step
    transportCoeffs.print();
    transportCoeffs.outputToJSON("variational_onsager_coefficients.json");

    if (mpi->mpiHead()) {
      std::cout << "Finished variational BTE solver\n\n";
      std::cout << std::string(80, '-') << "\n" << std::endl;
    }
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

    transportCoeffs.calcFromRelaxons(eigenvalues, eigenvectors,
                                     scatteringMatrix);
    transportCoeffs.print();
    transportCoeffs.outputToJSON("relaxons_onsager_coefficients.json");

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
  throwErrorIfUnset(context.getEpwFileName(), "epwFileName");
  throwErrorIfUnset(context.getTemperatures(), "temperatures");
  throwErrorIfUnset(context.getSmearingMethod(), "smearingMethod");
  if (context.getSmearingMethod() == DeltaFunction::gaussian) {
    throwErrorIfUnset(context.getSmearingWidth(), "smearingWidth");
  }

  if (context.getDopings().size() == 0 &&
      context.getChemicalPotentials().size() == 0) {
    Error e("Either chemical potentials or dopings must be set");
  }
}
