#include "lifetimes_app.h"
#include "bandstructure.h"
#include "context.h"
#include "el_scattering.h"
#include "ph_scattering.h"
#include "exceptions.h"
#include "io.h"
#include "path_points.h"
#include "qe_input_parser.h"
#include "ifc3_parser.h"

void ElectronLifetimesApp::run(Context &context) {
  context.setScatteringMatrixInMemory(false);

  auto t2 = QEParser::parsePhHarmonic(context);
  auto crystal = std::get<0>(t2);
  auto phononH0 = std::get<1>(t2);

  auto t1 = QEParser::parseElHarmonicWannier(context, &crystal);
  auto crystalEl = std::get<0>(t1);
  auto electronH0 = std::get<1>(t1);

  // load the el-ph coupling
  // Note: this file contains the number of electrons
  // which is needed to understand where to place the fermi level
  auto couplingElPh = InteractionElPhWan::parse(context, crystal, &phononH0);

  // set k and q point meshes and paths
  PathPoints pathKPoints(crystal, context.getPathExtrema(),
                         context.getDeltaPath());
  FullPoints fullKPoints(crystal, context.getKMesh());
  FullPoints fullQPoints(crystal, context.getQMesh());

  //----------------------------------------------------------------------------

  bool withVelocities = true;
  bool withEigenvectors = true;
  FullBandStructure fullElBandStructure =
      electronH0.populate(fullKPoints, withVelocities, withEigenvectors);
  FullBandStructure pathElBandStructure =
      electronH0.populate(pathKPoints, withVelocities, withEigenvectors);

  StatisticsSweep statisticsSweep(context, &fullElBandStructure);

  //----------------------------------------------------------------------------

  // build/initialize the scattering matrix and the smearing
  ElScatteringMatrix scatteringMatrix(context, statisticsSweep,
                                      fullElBandStructure, pathElBandStructure,
                                      phononH0, &couplingElPh);
  scatteringMatrix.setup();
  VectorBTE relaxationTimes = scatteringMatrix.getSingleModeTimes();

  mpi->barrier();
}

void PhononLifetimesApp::run(Context &context) {
  context.setScatteringMatrixInMemory(false);

  auto t2 = QEParser::parsePhHarmonic(context);
  auto crystal = std::get<0>(t2);
  auto phononH0 = std::get<1>(t2);

  // load the 3phonon coupling
  auto coupling3Ph = IFC3Parser::parse(context, crystal);

  // set k and q point meshes and paths
  PathPoints pathPoints(crystal, context.getPathExtrema(),
                        context.getDeltaPath());
  FullPoints fullPoints(crystal, context.getQMesh());

  //----------------------------------------------------------------------------

  bool withVelocities = true;
  bool withEigenvectors = true;
  FullBandStructure fullBandStructure =
      phononH0.populate(fullPoints, withVelocities, withEigenvectors);
  FullBandStructure pathBandStructure =
      phononH0.populate(pathPoints, withVelocities, withEigenvectors);

  StatisticsSweep statisticsSweep(context);

  //----------------------------------------------------------------------------

  // build/initialize the scattering matrix and the smearing
  PhScatteringMatrix scatteringMatrix(context, statisticsSweep,
                                      fullBandStructure, pathBandStructure,
                                      &coupling3Ph, &phononH0);
  scatteringMatrix.setup();
  VectorBTE relaxationTimes = scatteringMatrix.getSingleModeTimes();

  mpi->barrier();
}

void ElectronLifetimesApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");
  throwErrorIfUnset(context.getPhD2FileName(), "phD2FileName");
  throwErrorIfUnset(context.getPathExtrema(), "points path extrema");
  throwErrorIfUnset(context.getEpwFileName(), "epwFileName");
  throwErrorIfUnset(context.getKMesh(), "kMesh");
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

void PhononLifetimesApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getPhD2FileName(), "phD2FileName");
  throwWarningIfUnset(context.getSumRuleD2(), "sumRuleD2");
  throwErrorIfUnset(context.getPhD3FileName(), "phD3FileName");
  throwErrorIfUnset(context.getPathExtrema(), "points path extrema");
  throwErrorIfUnset(context.getKMesh(), "kMesh");
  throwErrorIfUnset(context.getTemperatures(), "temperatures");
  throwErrorIfUnset(context.getSmearingMethod(), "smearingMethod");
  if (context.getSmearingMethod() == DeltaFunction::gaussian) {
    throwErrorIfUnset(context.getSmearingWidth(), "smearingWidth");
  }
}
