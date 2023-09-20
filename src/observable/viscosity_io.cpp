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
                Eigen::Tensor<double, 5>& viscosityTensor, const bool& isPhonon, const bool& append,
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
  nlohmann::json output;
  if(append) { // we're going to add this to an existing one (as in the coupled case)
    std::ifstream f(outFileName);
    output = nlohmann::json::parse(f);
  }
  output["temperatures"] = temps;
  output[viscosityName] = viscosity;
  output["temperatureUnit"] = "K";
  output["viscosityUnit"] = units;
  //output["particleType"] = "phonon";
  std::ofstream o(outFileName);
  o << std::setw(3) << output << std::endl;
  o.close();
}
