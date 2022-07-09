#include "active_bandstructure.h"
#include "points.h"
#include "el_scattering.h"
#include "qe_input_parser.h"
#include <fstream>
#include <gtest/gtest.h>
#include "electron_wannier_transport_app.h"
#include "elph_qe_to_phoebe_app.h"
#include "parser.h"


/* TODO describe the test here
 */

TEST(MagneticUnfoldingTest, Test1) {

  // setup input file
  Context context;
  context.setPhFC2FileName("../test/data/silicon.fc");
  context.setElectronH0Name("../test/data/si_tb.dat");
  context.setWannier90Prefix("../test/data/si");
  context.setQuantumEspressoPrefix("../test/data/silicon");
  Eigen::VectorXd x3(1);
  x3(0) = 300. / temperatureAuToSi;
  context.setTemperatures(x3);
  context.setWindowType("magnetotransport");
  context.setScatteringMatrixInMemory(false);
  context.setUseSymmetries(true);
  Eigen::VectorXd dopings(1);
  dopings(0) = -1e14;
  context.setDopings(dopings);

  context.setSmearingMethod(0);
  context.setSmearingWidth(0.1/energyRyToEv);
  context.setFixedCouplingConstant(1.e-6);

  std::vector<Eigen::Vector3d> bfields;
  Eigen::Vector3d bfield = {1.0,0.,0.};
  bfields.push_back(bfield);
  context.setBField(bfields);

  auto t2 = Parser::parsePhHarmonic(context);
  auto crystal = std::get<0>(t2);
  auto phononH0 = std::get<1>(t2);

  auto t1 = Parser::parseElHarmonicWannier(context, &crystal);
  auto crystalEl = std::get<0>(t1);
  auto electronH0 = std::get<1>(t1);

  Eigen::Vector3i kMesh = {5,5,5};
  context.setKMesh(kMesh);

  context.setElphFileName("../test/data/silicon.phoebe.elph.hdf5");

  // load the elph coupling
  // Note: this file contains the number of electrons
  // which is needed to understand where to place the fermi level
  InteractionElPhWan couplingElPh = InteractionElPhWan::parse(context, crystal, &phononH0);
  Crystal fullCrystal = crystal;
  Crystal magCrystal = crystal;
  magCrystal.magneticSymmetries(bfield);

  Points fullSymPoints(fullCrystal, context.getKMesh());
  Points magSymPoints(magCrystal, context.getKMesh());

  // generate mag field linewidths
  auto t3 = ActiveBandStructure::builder(context, electronH0, magSymPoints);
  ActiveBandStructure magSymsBandStructure = std::get<0>(t3);
  auto magStatisticsSweep = std::get<1>(t3);

  // now the one fullSyms without magnetic symmetries
  auto t4 = ActiveBandStructure::builder(context, electronH0, fullSymPoints);
  ActiveBandStructure fullSymsBandStructure = std::get<0>(t4);
  auto fullStatisticsSweep = std::get<1>(t4);

  ASSERT_LT(magSymsBandStructure.getPoints().getCrystal().getNumSymmetries(),
            fullSymsBandStructure.getPoints().getCrystal().getNumSymmetries());

  // build/initialize the scattering matrix and the smearing
  ElScatteringMatrix fullScatteringMatrix(context, fullStatisticsSweep, fullSymsBandStructure,
                                          fullSymsBandStructure, phononH0, &couplingElPh);
  fullScatteringMatrix.setup();

  // unfold the linewidths from the calculation done with full symmetries to the set of states
  // with only the magnetic symmetries, to check that this is equal to the calculation done
  // with magnetic symmetries from the beginning
  Crystal crystalTemp = fullSymsBandStructure.getPoints().getCrystal();
  Points points(crystalTemp, context.getKMesh());

  unfoldLinewidths(context, fullScatteringMatrix, fullSymsBandStructure, fullStatisticsSweep,
                electronH0, points, bfield);

  VectorBTE fullSymsLinewidths = fullScatteringMatrix.getLinewidths();
  // after this call, fullSymsBandStructure is not anymore representative
  // of the system with the full symmetries without magnetic field.
  // now it should coincide with magSymsBandStructure

  ElScatteringMatrix magScatteringMatrix(context, fullStatisticsSweep, magSymsBandStructure,
                                         magSymsBandStructure, phononH0, &couplingElPh);
  magScatteringMatrix.setup();
  VectorBTE magSymLinewidths = magScatteringMatrix.getLinewidths();

  // now we compare the linewidths
  int numStatesMagSym = magSymLinewidths.getNumStates();
  int numStatesFullSyms = fullSymsLinewidths.getNumStates();
  ASSERT_EQ(numStatesMagSym, numStatesFullSyms);

  auto irrFullSyms = fullSymsBandStructure.irrStateIterator();
  auto irrMagSyms = magSymsBandStructure.irrStateIterator();
  ASSERT_EQ(irrFullSyms.size(), irrMagSyms.size());

  // check that the two band structures are exactly the same
  for (int is=0; is<numStatesMagSym; is++) {
    StateIndex isIdx(is);
    double en1 = fullSymsBandStructure.getEnergy(isIdx);
    double en2 = magSymsBandStructure.getEnergy(isIdx);
    ASSERT_EQ(en1, en2);
  }

  /// check the that the velocities make sense
  for (int is=0; is<numStatesMagSym; is++) {
    StateIndex isIdx(is);
    Eigen::Vector3d v2 = fullSymsBandStructure.getGroupVelocity(isIdx);
    Eigen::Vector3d v1 = magSymsBandStructure.getGroupVelocity(isIdx);
    double diff = (v1 - v2).norm();
    ASSERT_NEAR(diff, 0., 1e-8);

    auto t = magSymsBandStructure.getIndex(isIdx);
    WavevectorIndex ikIdx = std::get<0>(t);
    BandIndex ibIdx = std::get<1>(t);

    Eigen::Vector3d k1 = magSymsBandStructure.getWavevector(ikIdx);
    auto t2 = magSymsBandStructure.getRotationToIrreducible(k1, Points::cartesianCoordinates);
    int ikIrr = std::get<0>(t2);
    Eigen::Matrix3d rot = std::get<1>(t2);

    t2 = magSymsBandStructure.getRotationToIrreducible(k1, Points::cartesianCoordinates);
    ikIrr = std::get<0>(t2);
    rot = std::get<1>(t2);

    WavevectorIndex ikIrrIdx(ikIrr);
    Eigen::Vector3d kIrr = magSymsBandStructure.getWavevector(ikIrrIdx);
    auto v1sIrr = magSymsBandStructure.getGroupVelocities(ikIrrIdx);

    Eigen::Vector3d v1IrrReconstructed = rot * v1;
    Eigen::Vector3d k1IrrReconstructed = rot * k1;
    Eigen::Vector3d v1sIrrVector = v1sIrr.row(ibIdx.get());

    k1IrrReconstructed = magSymsBandStructure.getPoints().foldToBz(k1IrrReconstructed, Points::cartesianCoordinates);
    kIrr = magSymsBandStructure.getPoints().foldToBz(kIrr, Points::cartesianCoordinates);

    ASSERT_NEAR( (k1IrrReconstructed-kIrr).norm() , 0., 1e-4);
    ASSERT_NEAR( (v1IrrReconstructed-v1sIrrVector).norm() , 0., 1e-4);
}
  auto irrKPtsFullSyms = fullSymsBandStructure.irrPointsIterator();
  auto irrKPtsMagSyms = magSymsBandStructure.irrPointsIterator();
  ASSERT_EQ(irrKPtsFullSyms.size(), irrKPtsMagSyms.size());


  double maxLinewidth = magSymLinewidths.data.maxCoeff();
  double cutoff = maxLinewidth * 0.001;

  double lineWidthError = 0;
  // loop over irr points for each band structure
  for(int ik = 0; ik < int(magSymsBandStructure.irrPointsIterator().size()); ik++) {

    // the index of the irreducible point for each band structure in the MP list of k-points
    int irrKUnfold = irrKPtsFullSyms[ik];
    int irrKMagSyms = irrKPtsMagSyms[ik];
    ASSERT_EQ(irrKUnfold, irrKMagSyms);

    WavevectorIndex ikidx = WavevectorIndex(irrKUnfold);
    for(int ib = 0; ib < magSymsBandStructure.getNumBands(ikidx); ib++) {

      auto isIdxUnfold = StateIndex(fullSymsBandStructure.getIndex(
                WavevectorIndex(irrKUnfold),BandIndex(ib)));
      auto isIdxMagSym = StateIndex(magSymsBandStructure.getIndex(
                WavevectorIndex(irrKMagSyms),BandIndex(ib)));

      // check that the indices are the same
      ASSERT_EQ(isIdxUnfold.get(),isIdxMagSym.get());

      auto ibteUnfold = fullSymsBandStructure.stateToBte(isIdxUnfold).get();
      auto ibteMagSym = magSymsBandStructure.stateToBte(isIdxMagSym).get();
      ASSERT_EQ(ibteUnfold, ibteMagSym);

      // check the energies
      double enUnfold = fullSymsBandStructure.getEnergy(isIdxUnfold);
      double enMagSyms= magSymsBandStructure.getEnergy(isIdxMagSym);
      ASSERT_EQ(enUnfold, enMagSyms);

      // check the linewidths by relative error
      double linewidthUnfold = fullSymsLinewidths(0,0,ibteUnfold);
      double linewidthMagSyms = magSymLinewidths(0,0,ibteMagSym);
      double thisDiff = pow(linewidthUnfold-linewidthMagSyms,2);
      lineWidthError += thisDiff;
      // check ... anything else?
      if ( thisDiff > cutoff ) {
        thisDiff /= linewidthMagSyms;
        ASSERT_NEAR(thisDiff, 0, 1e-6);
      }
    }
  }
  ASSERT_NEAR(lineWidthError,0, 0.001);
}

TEST(MagneticUnfoldingTest, FullBandStructure) {
  // the test is similar to the above, but we don't discard active states

  // setup input file
  Context context;
  context.setPhFC2FileName("../test/data/silicon.fc");
  context.setElectronH0Name("../test/data/si_tb.dat");
  context.setWannier90Prefix("../test/data/si");
  context.setQuantumEspressoPrefix("../test/data/silicon");
  Eigen::VectorXd x3(1);
  x3(0) = 300. / temperatureAuToSi;
  context.setTemperatures(x3);
  context.setScatteringMatrixInMemory(false);
  context.setUseSymmetries(true);
  Eigen::VectorXd dopings(1);
  dopings(0) = -1e14;
  context.setDopings(dopings);
  context.setSmearingMethod(0);
  context.setSmearingWidth(0.1/energyRyToEv);
  context.setFixedCouplingConstant(1.e-4);

  std::vector<Eigen::Vector3d> bfields;
  Eigen::Vector3d bfield = {1.0,0.,0.};
  bfields.push_back(bfield);
  context.setBField(bfields);

  Eigen::Vector3i kMesh = {3,3,3};
  context.setKMesh(kMesh);

  context.setElphFileName("../test/data/silicon.phoebe.elph.hdf5");

  auto t2 = Parser::parsePhHarmonic(context);
  auto crystal = std::get<0>(t2);
  auto phononH0 = std::get<1>(t2);

  auto t1 = Parser::parseElHarmonicWannier(context, &crystal);
  auto crystalEl = std::get<0>(t1);
  auto electronH0 = std::get<1>(t1);

  // load the elph coupling
  // Note: this file contains the number of electrons
  // which is needed to understand where to place the fermi level
  InteractionElPhWan couplingElPh = InteractionElPhWan::parse(context, crystal, &phononH0);
  Crystal fullCrystal = crystal;
  Crystal magCrystal = crystal;
  magCrystal.magneticSymmetries(bfield);

  Points fullSymPoints(fullCrystal, context.getKMesh());
  Points magSymPoints(magCrystal, context.getKMesh());

  // generate mag field linewidths
  auto t3 = ActiveBandStructure::builder(context, electronH0, magSymPoints);
  ActiveBandStructure magSymsBandStructure = std::get<0>(t3);
  auto magStatisticsSweep = std::get<1>(t3);

  // now the one fullSyms without magnetic symmetries
  auto t4 = ActiveBandStructure::builder(context, electronH0, fullSymPoints);
  ActiveBandStructure fullSymsBandStructure = std::get<0>(t4);
  auto fullStatisticsSweep = std::get<1>(t4);

  ASSERT_LT(magSymsBandStructure.getPoints().getCrystal().getNumSymmetries(),
            fullSymsBandStructure.getPoints().getCrystal().getNumSymmetries());

  // build/initialize the scattering matrix and the smearing
  ElScatteringMatrix fullScatteringMatrix(context, fullStatisticsSweep, fullSymsBandStructure,
                                          fullSymsBandStructure, phononH0, &couplingElPh);
  fullScatteringMatrix.setup();

  Crystal crystalTemp = fullSymsBandStructure.getPoints().getCrystal();
  Points points(crystalTemp, context.getKMesh());

  // with only the magnetic symmetries, to check that this is equal to the calculation done
  // with magnetic symmetries from the beginning
  unfoldLinewidths(context, fullScatteringMatrix, fullSymsBandStructure, fullStatisticsSweep,
                 electronH0, points, bfield);
  VectorBTE fullSymsLinewidths = fullScatteringMatrix.getLinewidths();
  // after this call, fullSymsBandStructure is not anymore representative
  // of the system with the full symmetries without magnetic field.
  // now it should coincide with magSymBandStructure

  ElScatteringMatrix magScatteringMatrix(context, fullStatisticsSweep, magSymsBandStructure,
                                         magSymsBandStructure, phononH0, &couplingElPh);
  magScatteringMatrix.setup();
  VectorBTE magSymsLinewidths = magScatteringMatrix.getLinewidths();

  // now we compare the linewidths
  int numStatesMagSyms = magSymsLinewidths.getNumStates();
  int numStatesFullSyms = fullSymsLinewidths.getNumStates();
  ASSERT_EQ(numStatesMagSyms, numStatesFullSyms);

  auto irrFullSyms = fullSymsBandStructure.irrStateIterator();
  auto irrMagSyms = magSymsBandStructure.irrStateIterator();
  ASSERT_EQ(irrFullSyms.size(), irrMagSyms.size());

  // check that the two band structures are exactly the same
  for (int is=0; is< numStatesMagSyms; is++) {
    StateIndex isIdx(is);
    double en1 = fullSymsBandStructure.getEnergy(isIdx);
    double en2 = magSymsBandStructure.getEnergy(isIdx);
    ASSERT_EQ(en1, en2);
  }

  auto irrKPtsFullSyms = fullSymsBandStructure.irrPointsIterator();
  auto irrKPtsMagSyms = magSymsBandStructure.irrPointsIterator();
  ASSERT_EQ(irrKPtsFullSyms.size(), irrKPtsMagSyms.size());

  double maxLinewidth = magSymsLinewidths.data.maxCoeff();
  double cutoff = maxLinewidth * 0.001;

  // here we assert the equality between state and bte indices
  for(int ik = 0; ik < int(magSymsBandStructure.irrPointsIterator().size()); ik++) {
    // the index of the irreducible point for each band structure in the MP list of k-points
    int irrKUnfold = irrKPtsFullSyms[ik];
    int irrKMagSyms = irrKPtsMagSyms[ik];
    ASSERT_EQ(irrKUnfold, irrKMagSyms);
    WavevectorIndex ikidx = WavevectorIndex(irrKUnfold);
    for (int ib = 0; ib < magSymsBandStructure.getNumBands(ikidx); ib++) {
      auto isIdxUnfold = StateIndex(fullSymsBandStructure.getIndex(WavevectorIndex(irrKUnfold), BandIndex(ib)));
      auto isIdxMagSym = StateIndex(magSymsBandStructure.getIndex(WavevectorIndex(irrKMagSyms), BandIndex(ib)));
      // check that the indices are the same
      ASSERT_EQ(isIdxUnfold.get(), isIdxMagSym.get());
      auto ibteUnfold = fullSymsBandStructure.stateToBte(isIdxUnfold).get();
      auto ibteMagSym = magSymsBandStructure.stateToBte(isIdxMagSym).get();
      ASSERT_EQ(ibteUnfold, ibteMagSym);
    }
  }

  double lineWidthError = 0;

  // loop over irr points for each band structure
  for(int ik = 0; ik < int(magSymsBandStructure.irrPointsIterator().size()); ik++) {

    // the index of the irreducible point for each band structure in the MP list of k-points
    int irrKUnfold = irrKPtsFullSyms[ik];
    int irrKMagSyms = irrKPtsMagSyms[ik];
    ASSERT_EQ(irrKUnfold, irrKMagSyms);

    {
      WavevectorIndex ikIdx(ik);
      Eigen::VectorXd ens = fullSymsBandStructure.getEnergies(ikIdx);
      std::cout << ik << " " << ens.transpose() << "\n";
    }

    WavevectorIndex ikidx = WavevectorIndex(irrKUnfold);
    for(int ib = 0; ib < magSymsBandStructure.getNumBands(ikidx); ib++) {
      auto isIdxUnfold = StateIndex(fullSymsBandStructure.getIndex(
                        WavevectorIndex(irrKUnfold),BandIndex(ib)));
      auto isIdxMagSym = StateIndex(magSymsBandStructure.getIndex(
                        WavevectorIndex(irrKMagSyms),BandIndex(ib)));

      // check that the indices are the same
      ASSERT_EQ(isIdxUnfold.get(),isIdxMagSym.get());

      auto ibteUnfold = fullSymsBandStructure.stateToBte(isIdxUnfold).get();
      auto ibteMagSym = magSymsBandStructure.stateToBte(isIdxMagSym).get();
      ASSERT_EQ(ibteUnfold, ibteMagSym);

      // check the energies
      double enUnfold = fullSymsBandStructure.getEnergy(isIdxUnfold);
      double enMagSyms= magSymsBandStructure.getEnergy(isIdxMagSym);
      ASSERT_EQ(enUnfold, enMagSyms);

      // check the linewidths by relative error
      double linewidthUnfold = fullSymsLinewidths(0,0,ibteUnfold);
      double linewidthMagSyms = magSymsLinewidths(0,0,ibteMagSym);

      double thisDiff = pow(linewidthUnfold - linewidthMagSyms, 2);
      lineWidthError += thisDiff;
      // check linewidth error
      if ( thisDiff > cutoff ) {
        thisDiff /= linewidthMagSyms;
        ASSERT_NEAR(thisDiff, 0, 1e-3);
      }
    }
  }
  ASSERT_NEAR(lineWidthError, 0, 0.001);
}
