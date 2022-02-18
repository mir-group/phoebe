#include "ph_el_lifetimes.h"
#include "bandstructure.h"
#include "context.h"
#include "phel_scattering.h"
#include "exceptions.h"
#include "parser.h"

void PhElLifetimesApp::run(Context &context) {

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

  // compute the band structure on the fine grid
  if (mpi->mpiHead()) {
    std::cout << "\nComputing electronic band structure." << std::endl;
  }
  Points fullPoints(crystal, context.getQMesh());
  auto t3 = ActiveBandStructure::builder(context, electronH0, fullPoints);
  auto elBandStructure = std::get<0>(t3);
  auto statisticsSweep = std::get<1>(t3);

  // print some info about how window and symmetries have reduced things
  if (mpi->mpiHead()) {
    if(elBandStructure.hasWindow() != 0) {
        std::cout << "Window selection reduced electronic band structure from "
                << fullPoints.getNumPoints()*electronH0.getNumBands() << " to "
                << elBandStructure.getNumStates() << " states."  << std::endl;
    }
    if(context.getUseSymmetries()) {
      std::cout << "Symmetries reduced electronic band structure from "
          << elBandStructure.getNumStates() << " to "
          << elBandStructure.irrStateIterator().size() << " states." << std::endl;
    }
    std::cout << "Done computing electronic band structure.\n" << std::endl;
  }

  // Compute the full phonon band structure
  bool withVelocities = true;
  bool withEigenvectors = true;
  FullBandStructure phBandStructure = phononH0.populate(
      fullPoints, withVelocities, withEigenvectors);
  // set the chemical potentials to zero, load temperatures

  // build/initialize the scattering matrix and the smearing
  PhElScatteringMatrix scatteringMatrix(context, statisticsSweep,
                                        elBandStructure, phBandStructure,
                                        couplingElPh);
  scatteringMatrix.setup();
  scatteringMatrix.outputToJSON("rta_ph_el_relaxation_times.json");

  // solve the BTE at the relaxation time approximation level
  // we always do this, as it's the cheapest solver and is required to know
  // the diagonal for the exact method.

  if (mpi->mpiHead()) {
    std::cout << "\n" << std::string(80, '-') << "\n\n";
    std::cout << "Phonon-electron lifetimes computed.\n";
  }
}

void PhElLifetimesApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getPhFC2FileName(), "phFc2FileName");
  throwErrorIfUnset(context.getQMesh(), "qMesh");
  throwErrorIfUnset(context.getTemperatures(), "temperatures");
  throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");

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
