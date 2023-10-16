#include "viscosity_io.h"
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>

void printViscosity(std::string& viscosityName, Eigen::Tensor<double, 5>& viscosityTensor,
                                StatisticsSweep& statisticsSweep, int& dimensionality) {

  if (!mpi->mpiHead()) return;

  std::string units;
  if (dimensionality == 1)      { units = "Pa s / m^2"; } // 3d
  else if (dimensionality == 2) { units = "Pa s / m";   } // 2d
  else                          { units = "Pa s";       } // 1d

  int numCalculations = statisticsSweep.getNumCalculations();

  std::cout << "\n";
  std::cout << viscosityName << " viscosity (" << units << ")\n";
  std::cout << "i, j, k, eta[i,j,k,0], eta[i,j,k,1], eta[i,j,k,2]\n";

  for (int iCalc = 0; iCalc < numCalculations; iCalc++) {

    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStat.temperature;

    std::cout << std::fixed;
    std::cout.precision(2);
    std::cout << "Temperature: " << temp * temperatureAuToSi << " (K)\n";
    std::cout.precision(5);
    std::cout << std::scientific;
    for (int i = 0; i < dimensionality; i++) {
      for (int j = 0; j < dimensionality; j++) {
        for (int k = 0; k < dimensionality; k++) {
          std::cout << i << " " << j << " " << k;
          for (int l = 0; l < dimensionality; l++) {
            std::cout << " " << std::setw(12) << std::right
                      << viscosityTensor(iCalc, i, j, k, l) * viscosityAuToSi;
          }
          std::cout << "\n";
        }
      }
    }
    std::cout << std::endl;
  }
}


void outputViscosityToJSON(const std::string& outFileName, const std::string& viscosityName,
                Eigen::Tensor<double, 5>& viscosityTensor, const bool& append,
                StatisticsSweep& statisticsSweep, int& dimensionality) {

  if (!mpi->mpiHead()) return;

  int numCalculations = statisticsSweep.getNumCalculations();

  std::string units;
  if (dimensionality == 1) {      units = "Pa s / m^2"; } // 3d
  else if (dimensionality == 2) { units = "Pa s / m";   } // 2d
  else {                          units = "Pa s";       } // 1d

  std::vector<double> temps;
  // this vector mess is of shape (iCalculations, iRows, iColumns, k, l)
  std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> viscosity;

  for (int iCalc = 0; iCalc < numCalculations; iCalc++) {

    // store temperatures
    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStat.temperature;
    temps.push_back(temp * temperatureAuToSi);

    // store viscosity
    std::vector<std::vector<std::vector<std::vector<double>>>> rows;
    for (int i = 0; i < dimensionality; i++) {
      std::vector<std::vector<std::vector<double>>> cols;
      for (int j = 0; j < dimensionality; j++) {
        std::vector<std::vector<double>> ijk;
        for (int k = 0; k < dimensionality; k++) {
          std::vector<double> ijkl;
          for (int l = 0; l < dimensionality; l++) {
            ijkl.push_back(viscosityTensor(iCalc, i, j, k, l) * viscosityAuToSi);
          }
          ijk.push_back(ijkl);
        }
        cols.push_back(ijk);
      }
      rows.push_back(cols);
    }
    viscosity.push_back(rows);
  }

  // output to json
  // TODO need to update this with low dim viscosities
  nlohmann::json output;
  if(append) { // we're going to add this to an existing one (as in the coupled case)
    std::ifstream f(outFileName);
    output = nlohmann::json::parse(f);
  }
  output["temperatures"] = temps;
  output[viscosityName] = viscosity;
  output["temperatureUnit"] = "K";
  output["viscosityUnit"] = units;
  std::ofstream o(outFileName);
  o << std::setw(3) << output << std::endl;
  o.close();
}


// TODO we need to fix the dimensionality to work for
// low dim materials in all the coefficients!
void genericOutputRealSpaceToJSON(ScatteringMatrix& scatteringMatrix,
                                BaseBandStructure& bandStructure,
                                StatisticsSweep& statisticsSweep,
                                Eigen::VectorXd& theta0,
                                Eigen::VectorXd& theta_e,
                                Eigen::MatrixXd& phi,
                                double& C, Eigen::Vector3d& A) {

  // write D to file before diagonalizing, as the scattering matrix
  // will be destroyed by scalapack

  bool isPhonon = bandStructure.getParticle().isPhonon();
  int dimensionality = bandStructure.getPoints().getCrystal().getDimensionality();

  auto calcStat = statisticsSweep.getCalcStatistics(0); // only one calc for relaxons
  double kBT = calcStat.temperature;

  Eigen::MatrixXd Du(3,3); Du.setZero();
  Eigen::MatrixXd Wji0(3,3); Wji0.setZero();
  // below used only for electrons
  Eigen::MatrixXd Wjie(3,3); Wjie.setZero();

  // sum over the alpha and v states that this process owns
  for (auto tup : scatteringMatrix.getAllLocalStates()) {

    auto is1 = std::get<0>(tup);
    auto is2 = std::get<1>(tup);
    for (auto i : {0, 1, 2}) {
      for (auto j : {0, 1, 2}) {
        Du(i,j) += phi(i,is1) * scatteringMatrix(is1,is2) * phi(j,is2) * energyRyToEv;
      }
    }
  }
  mpi->allReduceSum(&Du);

  for (int is : bandStructure.parallelStateIterator()) {
    auto isIdx = StateIndex(is);
    double en = bandStructure.getEnergy(isIdx);
    // discard acoustic phonon modes
    if (en < 0.001 / ryToCmm1) { continue; }
    auto v = bandStructure.getGroupVelocity(isIdx);
    for (auto j : {0, 1, 2}) {
      for (auto i : {0, 1, 2}) {
        // note: phi and theta here are elStates long, so we need to shift the state
        // index to account for the fact that we summed over the electronic part above
        // calculate qunatities for the real-space solve
        Wji0(j,i) += phi(i,is) * v(j) * theta0(is) * velocityRyToSi;
        Wjie(j,i) += phi(i,is) * v(j) * theta_e(is) * velocityRyToSi;
      }
    }
  }
  mpi->allReduceSum(&Wji0); mpi->allReduceSum(&Wjie);

  // NOTE we cannot use nested vectors from the start, as
  // vector<vector> is not necessarily contiguous and MPI
  // cannot all reduce on it
  std::vector<std::vector<double>> vecDu;
  std::vector<std::vector<double>> vecWji0;
  std::vector<std::vector<double>> vecWjie;
  for (auto i : {0, 1, 2}) {
    std::vector<double> temp1,temp2,temp3;
    for (auto j : {0, 1, 2}) {
      temp1.push_back(Du(i,j));
      temp2.push_back(Wji0(i,j));
      temp3.push_back(Wjie(i,j));
    }
    vecDu.push_back(temp1);
    vecWji0.push_back(temp2);
    vecWjie.push_back(temp3);
  }

  // convert Ai to SI, in units of picograms/(mu m^3)
  double Aconversion = electronMassSi /
                       std::pow(distanceBohrToMum,dimensionality) * // convert AU mass / V -> SI
                       2 *   // factor of two is a Ry->Ha conversion required here
                       1e15; // convert electronMassSi in kg to pico g

                       // Michele's version of this, gives thes same answer
                       // double altConv =  1./rydbergSi * // convert kBT
                       // std::pow(hBarSi/bohrRadiusSi,2) * // convert (hbar * q)^2
                       // 1./std::pow(bohrRadiusSi, dimensionality) * // convert 1/V
                       // 1e-3; //convert from kg->pg, 1/m^3 -> 1/mum^3; // converting to pico and mu

  std::string specificHeatUnits;
  std::string AiUnits;
  if (dimensionality == 1) {
    specificHeatUnits = "J / K / m";
    AiUnits = "pg/(mum)";
  } else if (dimensionality == 2) {
    specificHeatUnits = "J / K / m^2";
    AiUnits = "pg/(mum)^2";
  } else {
    specificHeatUnits = "J / K / m^3";
    AiUnits = "pg/(mum)^3";
  }

  // this extra kBoltzmannRy is required when we calculate specific heat ...
  // TODO need to keep track of this and figure out where it's coming from
  double specificHeatConversion = kBoltzmannSi / pow(bohrRadiusSi, 3) / kBoltzmannRy;
  auto particle = bandStructure.getParticle();

  if(mpi->mpiHead()) {
    // output to json
    std::string outFileName = "el_relaxons_real_space_coeffs.json";
    if(isPhonon) outFileName = "ph_relaxons_real_space_coeffs.json";
    nlohmann::json output;
    output["temperature"] = kBT * temperatureAuToSi;
    output["Wji0"] = vecWji0;
    if(!isPhonon) output["Wjie"] = vecWjie;
    output["Du"] = vecDu;
    output["temperatureUnit"] = "K";
    output["wUnit"] = "m/s";
    output["DuUnit"] = "eV";
    output["specificHeat"] = C * specificHeatConversion;
    output["specificHeatUnit"] = specificHeatUnits;
    output["particleType"] = particle.isPhonon() ? "phonon" : "electron";
    std::vector<double> Atemp;
    for(int i = 0; i < 3; i++) {
      Atemp.push_back(A(i) * Aconversion );
    }
    output["Ai"] = Atemp;
    output["AiUnit"] = AiUnits;
    std::ofstream o(outFileName);
    o << std::setw(3) << output << std::endl;
    o.close();
  }
}

