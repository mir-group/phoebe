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
  context.setFixedCouplingConstant(1.e-4);

  Eigen::Vector3d bfield = {1.0,0.,0.};
  context.setBField(bfield);

  auto t2 = Parser::parsePhHarmonic(context);
  auto crystal = std::get<0>(t2);
  auto phononH0 = std::get<1>(t2);

  auto t1 = Parser::parseElHarmonicWannier(context, &crystal);
  auto crystalEl = std::get<0>(t1);
  auto electronH0 = std::get<1>(t1);

  Eigen::Vector3i kMesh = {3,3,3};

  int numModes = 3 * crystal.getNumAtoms();
  context.setKMesh(kMesh);

  context.setElphFileName("../test/data/silicon.phoebe.elph.hdf5");

  // load the elph coupling
  // Note: this file contains the number of electrons
  // which is needed to understand where to place the fermi level
  InteractionElPhWan couplingElPh = InteractionElPhWan::parse(context, crystal, &phononH0);
  Crystal fullCrystal = crystal;
  Crystal magCrystal = crystal;
  magCrystal.magneticSymmetries(context);

  Points fullSymPoints(fullCrystal, context.getKMesh());
  Points magSymPoints(magCrystal, context.getKMesh());

  // generate mag field linewidths
  auto t3 = ActiveBandStructure::builder(context, electronH0, magSymPoints);
  ActiveBandStructure magSymBandStructure = std::get<0>(t3);
  auto magStatisticsSweep = std::get<1>(t3);

  // now the one fullSyms without magnetic symmetries
  auto t4 = ActiveBandStructure::builder(context, electronH0, fullSymPoints);
  ActiveBandStructure fullSymsBandStructure = std::get<0>(t4);
  auto fullStatisticsSweep = std::get<1>(t4);
/*
  // make band structure with mag syms
  ActiveBandStructure magSymBandStructure = fullSymsBandStructure;
  // can we make a new copy of the points class so we don't mess up the old one?
  Points magSymPoints(fullCrystal, context.getKMesh());
  magSymPoints.magneticSymmetries(context);
  //magSymBandStructure.getPoints().magneticSymmetries(context);
  magSymBandStructure.swapPoints(magSymPoints);
  magSymBandStructure.rebuildSymmetries();
*/

  ASSERT_LT(magSymBandStructure.getPoints().getCrystal().getNumSymmetries(),
            fullSymsBandStructure.getPoints().getCrystal().getNumSymmetries());

  // build/initialize the scattering matrix and the smearing
  ElScatteringMatrix fullScatteringMatrix(context, fullStatisticsSweep, fullSymsBandStructure,
                                          fullSymsBandStructure, phononH0, &couplingElPh);
  fullScatteringMatrix.setup();

  // unfold the linewidths from the calculation done with full symmetries to the set of states
  // with only the magnetic symmetries, to check that this is equal to the calculation done
  // with magnetic symmetries from the beginning
  unfoldLinewidths(context, fullScatteringMatrix, fullSymsBandStructure, fullStatisticsSweep,electronH0);
  VectorBTE fullSymsLinewidths = fullScatteringMatrix.getLinewidths();
  // after this call, fullSymsBandStructure is not anymore representative
  // of the system with the full symmetries without magnetic field.
  // now it should coincide with magSymBandStructure

  ElScatteringMatrix magScatteringMatrix(context, fullStatisticsSweep, magSymBandStructure,
                                          magSymBandStructure, phononH0, &couplingElPh);
  magScatteringMatrix.setup();
  VectorBTE magSymLinewidths = magScatteringMatrix.getLinewidths();

  // NOTE: the number of active points may be different, due to some numerical
  // noise

  // now we compare the linewidths
  int numStatesMagSym = magSymLinewidths.getNumStates();
  int numStatesFullSyms = fullSymsLinewidths.getNumStates();
  ASSERT_EQ(numStatesMagSym, numStatesFullSyms);

  auto irrFullSyms = fullSymsBandStructure.irrStateIterator();
  auto irrMagSyms = magSymBandStructure.irrStateIterator();
  ASSERT_EQ(irrFullSyms.size(), irrMagSyms.size());

  // check that the two band structures are exactly the same
  for (int is=0; is<numStatesMagSym; is++) {
    StateIndex isIdx(is);
    double en1 = fullSymsBandStructure.getEnergy(isIdx);
    double en2 = magSymBandStructure.getEnergy(isIdx);
    ASSERT_EQ(en1, en2);
  }

  auto irrKPtsFullSyms = fullSymsBandStructure.irrPointsIterator();
  auto irrKPtsMagSyms = magSymBandStructure.irrPointsIterator();
  ASSERT_EQ(irrKPtsFullSyms.size(), irrKPtsMagSyms.size());


  double maxLinewidth = magSymLinewidths.data.maxCoeff(); //std::max_element(magSymLinewidths.n(), magSymLinewidths.data.end());
  double cutoff = maxLinewidth * 0.001;

  double lineWidthError = 0;
  // loop over irr points for each band structure
  for(int ik = 0; ik < int(magSymBandStructure.irrPointsIterator().size()); ik++) {

    // the index of the irreducible point for each band structure in the MP list of k-points
    int irrKUnfold = irrKPtsFullSyms[ik];
    int irrKMagSyms = irrKPtsMagSyms[ik];
    ASSERT_EQ(irrKUnfold, irrKMagSyms);

    WavevectorIndex ikidx = WavevectorIndex(irrKUnfold);
    for(int ib = 0; ib < magSymBandStructure.getNumBands(ikidx); ib++) {

      auto isIdxUnfold = StateIndex(fullSymsBandStructure.getIndex(WavevectorIndex(irrKUnfold),BandIndex(ib)));
      auto isIdxMagSym = StateIndex(magSymBandStructure.getIndex(WavevectorIndex(irrKMagSyms),BandIndex(ib)));

      //std::cout << "irrKUnfold irrKMagSyms " << irrKUnfold << " " << irrKMagSyms << std::endl;
      //std::cout << "stateUnfold stateMagSym " << isIdxUnfold.get() << " " << isIdxMagSym.get() << std::endl;
      // check that the indices are the same
      ASSERT_EQ(isIdxUnfold.get(),isIdxMagSym.get());

      auto ibteUnfold = fullSymsBandStructure.stateToBte(isIdxUnfold).get();
      auto ibteMagSym = magSymBandStructure.stateToBte(isIdxMagSym).get();
      ASSERT_EQ(ibteUnfold, ibteMagSym);
      // check the ibte states -- removed as this is not necessarily true
      //ASSERT_EQ(ibteUnfold, ibteMagSym);
      /*if(ibteUnfold != ibteMagSym)
      {  std::cout << "ibteUnfold ibteMagSym " << ibteUnfold << " " << ibteMagSym << std::endl;
        auto t1 = fullSymsBandStructure.getIndex(isIdxUnfold);
        int kunfold = std::get<0>(t1).get();
        int bunfold = std::get<1>(t1).get();
        auto t2 = magSymBandStructure.getIndex(isIdxMagSym);
        int kmag = std::get<0>(t2).get();
        int bmag = std::get<1>(t2).get();
        std::cout << "get Index " << kunfold << " " << bunfold << ", " << kmag << " " << bmag << std::endl;
      }*/
      // check the energies
      double enUnfold = fullSymsBandStructure.getEnergy(isIdxUnfold);
      double enMagSyms= magSymBandStructure.getEnergy(isIdxMagSym);
      //std::cout << "Energies "<< enUnfold << " " << enMagSyms << std::endl;
      //double relativeErr = (enUnfold - enMagSyms)/(enUnfold) * 100;
      ASSERT_EQ(enUnfold, enMagSyms);

      // check the linewidths by relative error
      double linewidthUnfold = fullSymsLinewidths(0,0,ibteUnfold);
      double linewidthMagSyms = magSymLinewidths(0,0,ibteMagSym);
      //std::cout << linewidthUnfold << " " << linewidthMagSyms << std::endl;
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
  //std::cout << "linewidthError" << lineWidthError << std::endl;

// temporary print statements for looking at irr points from each case
/*
  for(auto is : magSymBandStructure.irrPointsIterator()) {
    std::cout << is << " ";
  }
  std::cout << std::endl;
  std::cout << "now fullSyms" << std::endl;
  for(auto is : fullSymsBandStructure.irrPointsIterator()) {
    std::cout << is << " ";
  }
  std::cout << std::endl;
  std::cout << "now fullSyms from smatrix" << std::endl;
  for(auto is : fullScatteringMatrix.getBandStructure().irrPointsIterator()) {
    std::cout << is << " " ;
  }
  std::cout << std::endl;
*/

}

TEST(MagneticUnfoldingTest, FullBS) {
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

  Eigen::Vector3d bfield = {1.0,0.,0.};
  context.setBField(bfield);

  Eigen::Vector3i kMesh = {3,3,3};
  context.setKMesh(kMesh);

  context.setElphFileName("../test/data/silicon.phoebe.elph.hdf5");

  auto t2 = Parser::parsePhHarmonic(context);
  auto crystal = std::get<0>(t2);
  auto phononH0 = std::get<1>(t2);

  auto t1 = Parser::parseElHarmonicWannier(context, &crystal);
  auto crystalEl = std::get<0>(t1);
  auto electronH0 = std::get<1>(t1);

  int numModes = 3 * crystal.getNumAtoms();

  // load the elph coupling
  // Note: this file contains the number of electrons
  // which is needed to understand where to place the fermi level
  InteractionElPhWan couplingElPh = InteractionElPhWan::parse(context, crystal, &phononH0);
  Crystal fullCrystal = crystal;
  Crystal magCrystal = crystal;
  magCrystal.magneticSymmetries(context);

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
  unfoldLinewidths(context, fullScatteringMatrix, fullSymsBandStructure, fullStatisticsSweep,electronH0);
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

  double maxLinewidth = magSymsLinewidths.data.maxCoeff(); //std::max_element(magSymsLinewidths.n(), magSymsLinewidths.data.end());
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

//  {
//    for (int ik = 0; ik < int(magSymsBandStructure.irrPointsIterator().size()); ik++) {
//      // the index of the irreducible point for each band structure in the MP list of k-points
//      int irrKUnfold = irrKPtsFullSyms[ik];
//      int irrKMagSyms = irrKPtsMagSyms[ik];
//      ASSERT_EQ(irrKUnfold, irrKMagSyms);
//
//      WavevectorIndex ikidx = WavevectorIndex(irrKUnfold);
//      int ib = 0;
//      auto isIdxUnfold = StateIndex(fullSymsBandStructure.getIndex(WavevectorIndex(irrKUnfold), BandIndex(ib)));
//      auto isIdxMagSym = StateIndex(magSymsBandStructure.getIndex(WavevectorIndex(irrKMagSyms), BandIndex(ib)));
//      auto ibteUnfold = fullSymsBandStructure.stateToBte(isIdxUnfold).get();
//      auto ibteMagSym = magSymsBandStructure.stateToBte(isIdxMagSym).get();
//
//      // check the linewidths by relative error
//      double linewidthUnfold = fullSymsLinewidths(0, 0, ibteUnfold);
//      double linewidthMagSyms = magSymsLinewidths(0, 0, ibteMagSym);
//      std::cout << ik << " " << linewidthUnfold << " " << linewidthMagSyms << "\n";
//    }
//    for (int ik : magSymsBandStructure.irrPointsIterator()) {
//      WavevectorIndex ikIdx = WavevectorIndex(ik);
//      int ib = 0;
//      auto isIdx = StateIndex(magSymsBandStructure.getIndex(WavevectorIndex(ikIdx), BandIndex(ib)));
//      auto ibte = magSymsBandStructure.stateToBte(isIdx).get();
//      std::cout << ik << " " << magSymsLinewidths(0, 0, ibte) << " - m\n";
//    }
//    for (int ik : fullSymsBandStructure.irrPointsIterator()) {
//      WavevectorIndex ikIdx = WavevectorIndex(ik);
//      int ib = 0;
//      auto isIdx = StateIndex(fullSymsBandStructure.getIndex(WavevectorIndex(ikIdx), BandIndex(ib)));
//      auto ibte = fullSymsBandStructure.stateToBte(isIdx).get();
//      std::cout << ik << " " << fullSymsLinewidths(0, 0, ibte) << " - m\n";
//    }
//  }

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

      auto isIdxUnfold = StateIndex(fullSymsBandStructure.getIndex(WavevectorIndex(irrKUnfold),BandIndex(ib)));
      auto isIdxMagSym = StateIndex(magSymsBandStructure.getIndex(WavevectorIndex(irrKMagSyms),BandIndex(ib)));

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

      std::cout << ik << " " << ib << " "
                << linewidthUnfold << " " << linewidthMagSyms << "\n";

      //std::cout << linewidthUnfold << " " << linewidthMagSyms << std::endl;
      double thisDiff = pow(linewidthUnfold - linewidthMagSyms, 2);
      lineWidthError += thisDiff;
      // check ... anything else?
      if ( thisDiff > cutoff ) {
        thisDiff /= linewidthMagSyms;
        ASSERT_NEAR(thisDiff, 0, 1e-3);
      }
    }
  }

  {
//    double m = fullSymsLinewidths.data.mean();
//    std::cout << m << "!\n";
    ASSERT_NEAR(lineWidthError, 0, 0.001);
  }
}
