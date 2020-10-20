#include "elph_plot_app.h"
#include "bandstructure.h"
#include "constants.h"
#include "context.h"
#include "drift.h"
#include "exceptions.h"
#include "io.h"
#include "observable.h"
#include "particle.h"
#include "el_scattering.h"
#include "onsager.h"
#include "qe_input_parser.h"
#include "specific_heat.h"
#include "path_points.h"

void ElPhCouplingPlotApp::run(Context &context) {

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

  std::vector<std::pair<Eigen::Vector3d,Eigen::Vector3d>> pointsPath;

  // create a list of (k,q) pairs, where k(q) is on a path and q(k) is fixed
  if ( context.getG2PlotStyle() == "fixedQ" ) {

    PathPoints kPoints(crystal, context.getPathExtrema(), context.getDeltaPath());
    for (long ik = 0; ik < kPoints.getNumPoints(); ik++) {
      auto thisK = kPoints.getPointCoords(ik, Points::cartesianCoords);
      std::pair<Eigen::Vector3d, Eigen::Vector3d> thisPair;
      thisPair.first = thisK;
      thisPair.second = context.getG2PlotFixedPoint();
      pointsPath.push_back(thisPair);
    }

  } else {

    PathPoints qPoints(crystal, context.getPathExtrema(), context.getDeltaPath());
    for (long iq = 0; iq < qPoints.getNumPoints(); iq++) {
      auto thisQ = qPoints.getPointCoords(iq, Points::cartesianCoords);
      std::pair<Eigen::Vector3d, Eigen::Vector3d> thisPair;
      thisPair.first = context.getG2PlotFixedPoint();
      thisPair.second = thisQ;
      pointsPath.push_back(thisPair);
    }

  }

  std::pair<int,int> g2PlotEl1Bands = context.getG2PlotEl1Bands();
  std::pair<int,int> g2PlotEl2Bands = context.getG2PlotEl2Bands();
  std::pair<int,int> g2PlotPhBands = context.getG2PlotPhBands();

  // Compute the coupling

  std::vector<double> allGs;
  for (auto thisPair : pointsPath) {
    Eigen::Vector3d k1C = thisPair.first;
    Eigen::Vector3d q3C = thisPair.second;
    Eigen::Vector3d k2C = k1C + q3C;

    // I need to get the eigenvectors at these three wavevectors

    auto t1 = electronH0.diagonalizeFromCoords(k1C);
    auto eigvec1 = std::get<1>(t1);

    auto t2 = electronH0.diagonalizeFromCoords(k2C);
    auto eigvec2 = std::get<1>(t2);
    std::vector<Eigen::MatrixXcd> eigvecs2;
    eigvecs2.push_back(eigvec2);
    std::vector<Eigen::Vector3d> k2Cs;
    k2Cs.push_back(k2C);

    auto t3 = phononH0.diagonalizeFromCoords(q3C);
    auto eigvec3 = std::get<1>(t3);
    std::vector<Eigen::MatrixXcd> eigvecs3;
    eigvecs3.push_back(eigvec3);
    std::vector<Eigen::Vector3d> q3Cs;
    q3Cs.push_back(q3C);

    couplingElPh.calcCouplingSquared(eigvec1, eigvecs2, eigvecs3, k1C, k2Cs,
                                     q3Cs);
    auto coupling = couplingElPh.getCouplingSquared(0);

    double sum1 = 0.;
    for (int ib1 = g2PlotEl1Bands.first; ib1 <= g2PlotEl1Bands.second; ib1++) {
      for (int ib2 = g2PlotEl2Bands.first; ib2 <= g2PlotEl2Bands.second;
           ib2++) {
        for (int ib3 = g2PlotPhBands.first; ib3 <= g2PlotPhBands.second;
             ib3++) {
          sum1 += coupling(ib1, ib2, ib3);
        }
      }
    }
    allGs.push_back(sum1);
  }

  if (mpi->mpiHead()) {
    std::ofstream outfile("./g2_coupling_plot.dat");
    int i = 0;
    for ( auto x : allGs ) {
      outfile << i << "\t" << x << "\n";
      i += 1;
    }
  }


  if ( mpi->mpiHead()) {
    std::cout << "\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << "\n";
  }

  mpi->barrier();
}

void ElPhCouplingPlotApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");
  throwErrorIfUnset(context.getPhD2FileName(), "phD2FileName");
  throwErrorIfUnset(context.getPathExtrema(), "points path extrema");
  throwErrorIfUnset(context.getEpwFileName(), "epwFileName");
}
