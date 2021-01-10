#include "elph_plot_app.h"
#include "bandstructure.h"
#include "context.h"
#include "el_scattering.h"
#include "exceptions.h"
#include "io.h"
#include "points.h"
#include "qe_input_parser.h"

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

  std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> pointsPath;

  // create a list of (k,q) pairs, where k(q) is on a path and q(k) is fixed
  if (context.getG2PlotStyle() == "fixedQ") {

    Points kPoints(crystal, context.getPathExtrema(),
                       context.getDeltaPath());
    for (int ik = 0; ik < kPoints.getNumPoints(); ik++) {
      auto thisK =
          kPoints.getPointCoordinates(ik, Points::cartesianCoordinates);
      std::pair<Eigen::Vector3d, Eigen::Vector3d> thisPair;
      thisPair.first = thisK;
      thisPair.second = context.getG2PlotFixedPoint();
      pointsPath.push_back(thisPair);
    }

  } else {

    Points qPoints(crystal, context.getPathExtrema(),
                       context.getDeltaPath());
    for (int iq = 0; iq < qPoints.getNumPoints(); iq++) {
      auto thisQ =
          qPoints.getPointCoordinates(iq, Points::cartesianCoordinates);
      std::pair<Eigen::Vector3d, Eigen::Vector3d> thisPair;
      thisPair.first = context.getG2PlotFixedPoint();
      thisPair.second = thisQ;
      pointsPath.push_back(thisPair);
    }
  }

  std::pair<int, int> g2PlotEl1Bands = context.getG2PlotEl1Bands();
  std::pair<int, int> g2PlotEl2Bands = context.getG2PlotEl2Bands();
  std::pair<int, int> g2PlotPhBands = context.getG2PlotPhBands();

  // Compute the coupling

  std::vector<double> allGs;
  for (const auto &thisPair : pointsPath) {
    Eigen::Vector3d k1C = thisPair.first;
    Eigen::Vector3d q3C = thisPair.second;
    Eigen::Vector3d k2C = k1C + q3C;

    // I need to get the eigenvectors at these three wavevectors

    auto t3 = electronH0.diagonalizeFromCoordinates(k1C);
    auto eigenVector1 = std::get<1>(t3);

    auto t4 = electronH0.diagonalizeFromCoordinates(k2C);
    auto eigenVector2 = std::get<1>(t4);
    std::vector<Eigen::MatrixXcd> eigenVectors2;
    eigenVectors2.push_back(eigenVector2);
    std::vector<Eigen::Vector3d> k2Cs;
    k2Cs.push_back(k2C);

    auto t5 = phononH0.diagonalizeFromCoordinates(q3C);
    auto eigenVector3 = std::get<1>(t5);
    std::vector<Eigen::MatrixXcd> eigenVectors3;
    eigenVectors3.push_back(eigenVector3);
    std::vector<Eigen::Vector3d> q3Cs;
    q3Cs.push_back(q3C);

    couplingElPh.calcCouplingSquared(eigenVector1, eigenVectors2, eigenVectors3,
                                     k1C, k2Cs, q3Cs);
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
    for (auto x : allGs) {
      outfile << i << "\t" << x << "\n";
      i += 1;
    }
  }

  if (mpi->mpiHead()) {
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
