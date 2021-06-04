#include "elph_qe_to_phoebe_app.h"
#include "bandstructure.h"
#include "eigen.h"
#include "interaction_elph.h"
#include "io.h"
#include "qe_input_parser.h"
#include <iomanip>
#include <sstream>
#include <string>
#include <exception>

#ifdef HDF5_AVAIL
#include <highfive/H5Easy.hpp>
#endif

void ElPhQeToPhoebeApp::epaPostProcessing(Context &context, Eigen::MatrixXd &elEnergies,
                       Points &kPoints, Points &qPoints,
                       const int &numElectrons, const int &numSpin,
                       const int &numModes, const int &numIrrQPoints,
                       const int &numQEBands, const Eigen::MatrixXd &energies,
                       const Eigen::MatrixXd &kGridFull) {

  if (mpi->mpiHead()) {
    std::cout << "Starting EPA post-processing\n" << std::endl;
  }

  auto t2 = QEParser::parseElHarmonicFourier(context);
  auto electronH0 = std::get<1>(t2);
  int numBands = int(electronH0.getNumBands());

  // read coupling from file
  auto t5 =
      readGFromQEFile(context, numModes, numBands, numBands, kPoints, qPoints,
                      kGridFull, numIrrQPoints, numQEBands, energies);
  auto gFull = std::get<0>(t5);          // (nBands, nBands, nModes, numK, numQ)
  auto phEigenvectors = std::get<1>(t5); // (numModes, numModes, numQPoints)
  auto phEnergies = std::get<2>(t5);     // (numModes, numQPoints)

  assert(numBands == gFull.dimension(0));
  assert(numModes == gFull.dimension(2));

  // input
  double smearing = context.getEpaSmearingEnergy();
  double smearing2 = 2. * smearing * smearing;

  // prepare energy bins
  double minEnergy = context.getEpaMinEnergy();
  double maxEnergy = context.getEpaMaxEnergy();
  if ( maxEnergy < minEnergy ) {
    Error("Problems in setting the EPA energy ranges");
  }

  double deltaEnergy = context.getEpaDeltaEnergy();
  int numEpaEnergies = context.getEpaNumBins();
  if ( std::isnan(deltaEnergy)) {
    context.getEpaNumBins();
    deltaEnergy = (maxEnergy-minEnergy)/numEpaEnergies;
  } else {
    numEpaEnergies = int( (maxEnergy - minEnergy) / deltaEnergy ) + 1;
  }

  Eigen::VectorXd epaEnergies(numEpaEnergies);
#pragma omp parallel for default(none) shared(numEpaEnergies, epaEnergies, deltaEnergy, minEnergy)
  for (int i = 0; i < numEpaEnergies; i++) {
    epaEnergies[i] = i * deltaEnergy + minEnergy;
  }

  if (mpi->mpiHead()) {
    std::cout << "Building EPA with " << numEpaEnergies << " energy bins.";
  }

  int numKPoints = gFull.dimension(3);
  int numQPoints = gFull.dimension(4);

  Eigen::Tensor<double, 3> gaussian(numEpaEnergies, numBands, numKPoints);
#pragma omp parallel for collapse(3) default(none) shared(numBands,numKPoints, numEpaEnergies, elEnergies, epaEnergies, smearing2, gaussian)
  for (int ib1 = 0; ib1 < numBands; ib1++) {
    for (int ik = 0; ik < numKPoints; ik++) {
      for (int i = 0; i < numEpaEnergies; i++) {
        double arg = pow(elEnergies(ib1, ik) - epaEnergies(i), 2) / smearing2;
        gaussian(i, ib1, ik) = exp(-arg);
      }
    }
  }

  Eigen::Tensor<double,5> g2Full(numBands,numBands,numModes,numKPoints,numQPoints);
  for (int iq = 0; iq < numQPoints; iq++) {
    for (int ik = 0; ik < numKPoints; ik++) {
      for (int nu = 0; nu < numModes; nu++) {
        for (int ib2 = 0; ib2 < numBands; ib2++) {
          for (int ib1 = 0; ib1 < numBands; ib1++) {
            g2Full(ib1, ib2, nu, ik, iq) =
                std::norm(gFull(ib1, ib2, nu, ik, iq));
          }
        }
      }
    }
  }
  gFull.resize(0,0,0,0,0);

  Eigen::Tensor<double, 3> g2Epa(numModes, numEpaEnergies, numEpaEnergies);
  g2Epa.setZero();

  LoopPrint loopPrint("Computing coupling EPA", "q-points", numQPoints);
  for (int iq : mpi->divideWorkIter(numQPoints)) {
    loopPrint.update();
    Eigen::Vector3d q =
        qPoints.getPointCoordinates(iq, Points::cartesianCoordinates);
    for (int ik = 0; ik < numKPoints; ik++) {
      Eigen::Vector3d k =
          kPoints.getPointCoordinates(ik, Points::cartesianCoordinates);

      // Coordinates and index of k+q point
      Eigen::Vector3d kq = k + q;
      Eigen::Vector3d kqCrystal = kPoints.cartesianToCrystal(kq);
      int ikq = int(kPoints.getIndex(kqCrystal));

 #pragma omp parallel for collapse(3) default(none) shared(numEpaEnergies, numModes, numBands, gaussian, g2Full, phEnergies, g2Epa, ik, ikq, iq)
      for (int j = 0; j < numEpaEnergies; j++) {
        for (int i = 0; i < numEpaEnergies; i++) {
          for (int nu = 0; nu < numModes; nu++) {

            for (int ib2 = 0; ib2 < numBands; ib2++) {
              for (int ib1 = 0; ib1 < numBands; ib1++) {

              double gaussianX = gaussian(i, ib1, ik) * gaussian(j, ib2, ikq);

                g2Epa(nu, i, j) += g2Full(ib1, ib2, nu, ik, iq) *
                    gaussianX / 2. / phEnergies(nu, iq);
                // /2omega, because there is a difference between the
                // coupling <k+q| dV_q |k> from quantum espresso
                // and the coupling g to be used for transport calculations
              }
            }
          }
        }
      }
    }
  }
  mpi->allReduceSum(&g2Epa);
  loopPrint.close();

  Eigen::VectorXd phAvgEnergies(numModes); // phEnergies(numModes, numQPoints);
  for (int nu = 0; nu < numModes; nu++) {
    phAvgEnergies(nu) = phEnergies.row(nu).sum() / phEnergies.cols();
  }

  if (mpi->mpiHead()) {
    std::cout << "\nStart writing el-ph coupling to file." << std::endl;
    std::string phoebePrefixQE = context.getQuantumEspressoPrefix();
    std::string outFileName = phoebePrefixQE + ".phoebe.epa.dat";
    std::ofstream outfile(outFileName);
    if (not outfile.is_open()) {
      Error("Output file couldn't be opened");
    }
    outfile << numElectrons << " " << numSpin << "\n";
    outfile << phAvgEnergies.size() << "\n";
    outfile << phAvgEnergies.transpose() << "\n";
    outfile << numEpaEnergies << "\n";
    outfile << epaEnergies.transpose() << "\n";
    for (auto i = 0; i < numModes; ++i) {
      for (auto j = 0; j < numEpaEnergies; ++j) {
        for (auto k = 0; k < numEpaEnergies; ++k) {
          outfile << g2Epa(i, j, k) << "\n";
        }
      }
    }
    std::cout << "Done writing el-ph coupling to file.\n" << std::endl;
  }
}
