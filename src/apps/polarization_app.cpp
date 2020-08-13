#include "polarization_app.h"

#include <fstream>
#include <string>

#include "constants.h"
#include "exceptions.h"
#include "periodic_table.h"
#include "qe_input_parser.h"
#include "statistics_sweep.h"

// Compute the electronic polarization using the Berry connection
void ElectronPolarizationApp::run(Context &context) {
  std::cout << "Starting electron polarization calculation" << std::endl;

  // Read the necessary input files
  auto tup = QEParser::parseElHarmonicWannier(context);
  auto crystal = std::get<0>(tup);
  auto h0 = std::get<1>(tup);

  // first we make compute the band structure on the fine grid
  FullPoints points(crystal, context.getKMesh());
  bool withVelocities = false;
  bool withEigenvectors = false;
  FullBandStructure bandStructure =
      h0.populate(points, withVelocities, withEigenvectors);

  // now we build the Berry connection

  Eigen::Tensor<double, 3> berryConnection(points.getNumPoints(),
                                           h0.getNumBands(), 3);
  berryConnection.setZero();

  for (long ik = 0; ik < points.getNumPoints(); ik++) {
    auto point = points.getPoint(ik);
    std::vector<Eigen::MatrixXcd> thisBerryConnection =
        h0.getBerryConnection(point);
    for (long ib = 0; ib < h0.getNumBands(); ib++) {
      for (long i = 0; i < 3; i++) {
        auto x = thisBerryConnection[i](ib, ib);
        berryConnection(ik, ib, i) = x.real();
      }
    }
  }

  Particle particle = h0.getParticle();

  // before moving on, we need to fix the chemical potential
  StatisticsSweep statisticsSweep(context, &bandStructure);
  auto numCalcs = statisticsSweep.getNumCalcs();

  // now we can compute the polarization

  Eigen::MatrixXd polarization(numCalcs, 3);
  polarization.setZero();

  for (long is = 0; is < bandStructure.getNumStates(); is++) {
    double energy = bandStructure.getEnergy(is);
    auto t = bandStructure.getIndex(is);
    int ik = std::get<0>(t).get();
    int ib = std::get<0>(t).get();

    for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
      auto sc = statisticsSweep.getCalcStatistics(iCalc);
      double temp = sc.temperature;
      double chemPot = sc.chemicalPotential;
      double population = particle.getPopulation(energy, temp, chemPot);
      for (long i = 0; i < 3; i++) {
        polarization(iCalc, i) -= population * berryConnection(ik, ib, i);
      }
    }
  }

  double volume = crystal.getVolumeUnitCell();
  polarization.array() /= points.getNumPoints() * volume;

  // now we add the ionic term

  PeriodicTable periodicTable;
  Eigen::MatrixXd atomicPositions = crystal.getAtomicPositions();
  std::vector<std::string> atomicNames = crystal.getAtomicNames();
  auto numAtoms = crystal.getNumAtoms();
  for (long ia = 0; ia < numAtoms; ia++) {
    Eigen::Vector3d position = atomicPositions.row(ia);
    auto charge = periodicTable.getIonicCharge(atomicNames[ia]);
    for (long i = 0; i < 3; i++) {
      for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
        polarization(iCalc, i) += charge * position(i);
      }
    }
  }
  polarization.array() /= volume;

  // Save results to file
  std::ofstream outfile("./polarization.dat");
  outfile << "# Electrical polarization density: "
             "chemical potential (eV), doping (cm^-3), temperature (K)"
             "polarization[x,y,z] (a.u.)\n";
  for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
    auto sc = statisticsSweep.getCalcStatistics(iCalc);
    auto temp = sc.temperature;
    auto chemPot = sc.chemicalPotential;
    auto doping = sc.doping;
    outfile << chemPot * energyRyToEv << "\t" << doping;
    outfile << "\t" << temp * temperatureAuToSi;
    for (long i = 0; i < 3; i++) {
      outfile << "\t" << polarization(iCalc, i);
    }
    outfile << "\n";
  }

  std::cout << "Electron polarization computed" << std::endl;
}
