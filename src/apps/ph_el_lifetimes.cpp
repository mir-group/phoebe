#include "ph_el_lifetimes.h"
#include "bandstructure.h"
#include "context.h"
#include "phel_scattering.h"
#include "exceptions.h"
#include "parser.h"

void PhElLifetimesApp::run(Context &context) {

  // load harmonic hamiltonians
  auto t2 = Parser::parsePhHarmonic(context);
  auto crystal = std::get<0>(t2);
  auto phononH0 = std::get<1>(t2);

  auto t1 = Parser::parseElHarmonicWannier(context, &crystal);
  auto crystalEl = std::get<0>(t1);
  auto electronH0 = std::get<1>(t1);

  // load the elph coupling
  // Note: this file contains the number of electrons
  // which is needed to understand where to place the fermi level
  auto couplingElPh = InteractionElPhWan::parse(context, crystal, &phononH0);

  // Compute the window filtered phonon band structure ---------------------------
  if (mpi->mpiHead()) {
    std::cout << "\nComputing phonon band structure." << std::endl;
  }

  bool withVelocities = true;
  bool withEigenvectors = true;
  Points qPoints(crystal, context.getQMesh());
  auto t3 = ActiveBandStructure::builder(context, phononH0, qPoints);
  auto phBandStructure = std::get<0>(t3);

  // print some info about how window and symmetries have reduced things
  if (mpi->mpiHead()) {
    if(phBandStructure.hasWindow() != 0) {
        std::cout << "Window selection reduced phonon band structure from "
          << qPoints.getNumPoints() * phononH0.getNumBands() << " to "
          << phBandStructure.getNumStates() << " states."  << std::endl;
    }
    if(context.getUseSymmetries()) {
      std::cout << "Symmetries reduced phonon band structure from "
          << phBandStructure.getNumStates() << " to "
          << phBandStructure.irrStateIterator().size() << " states.\n" << std::endl;
    }
    std::cout << "Done computing phonon band structure.\n" << std::endl;
  }

  // compute the el band structure on the fine grid -----------------------------
  if (mpi->mpiHead()) {
    std::cout << "Computing electronic band structure.\n" << std::endl;
  }

  // manually setting the window to 1.25 the maximum phonon
  double maxPhEnergy = phBandStructure.getMaxEnergy();
  auto inputWindowType = context.getWindowType();
  context.setWindowType("muCenteredEnergy");
  if(mpi->mpiHead()) {
    std::cout << "Of the active phonon modes, the maximum energy state is " <<
        maxPhEnergy * energyRyToEv * 1e3 << " meV." <<
        "\nSelecting states within +/- 1.25 x " << maxPhEnergy*energyRyToEv*1e3 << " meV"
        << " of max/min electronic mu values." << std::endl;
  }
  Eigen::Vector2d range = {-1.25*maxPhEnergy, 1.25*maxPhEnergy};
  context.setWindowEnergyLimit(range);

  Points kPoints(crystal, context.getKMesh());
  auto t4 = ActiveBandStructure::builder(context, electronH0, kPoints);
  auto elBandStructure = std::get<0>(t4);
  auto statisticsSweep = std::get<1>(t4);

  // print some info about how window and symmetries have reduced things
  if (mpi->mpiHead()) {
    if(elBandStructure.hasWindow() != 0) {
        std::cout << "Window selection reduced electronic band structure from "
          << kPoints.getNumPoints() * electronH0.getNumBands() << " to "
          << elBandStructure.getNumStates() << " states."  << std::endl;
    }
    if(context.getUseSymmetries()) {
      std::cout << "Symmetries reduced electronic band structure from "
          << elBandStructure.getNumStates() << " to "
          << elBandStructure.irrStateIterator().size() << " states.\n" << std::endl;
    }
    std::cout << "Done computing electronic band structure.\n" << std::endl;
  }

  // build/initialize the scattering matrix and the smearing
  PhElScatteringMatrix scatteringMatrix(context, statisticsSweep,
                                        elBandStructure, phBandStructure,
                                        couplingElPh, electronH0);
  scatteringMatrix.setup();
  scatteringMatrix.outputToJSON("rta_phel_relaxation_times.json");

  // solve the BTE at the relaxation time approximation level
  // we always do this, as it's the cheapest solver and is required to know
  // the diagonal for the exact method.

  if (mpi->mpiHead()) {
    std::cout << "\n" << std::string(80, '-') << "\n\n";
    std::cout << "Phonon-electron lifetimes computed.\n";
  }
}

void PhElLifetimesApp::checkRequirements(Context &context) {

  throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");
  throwErrorIfUnset(context.getPhFC2FileName(), "phFC2FileName");
  throwErrorIfUnset(context.getQMesh(), "qMesh");
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
