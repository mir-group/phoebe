#include "phonon_thermal_cond.h"

#include <time.h>
#include <nlohmann/json.hpp>
#include <fstream>
#include <iomanip>
#include "constants.h"
#include "mpiHelper.h"

PhononThermalConductivity::PhononThermalConductivity(
    StatisticsSweep &statisticsSweep_, Crystal &crystal_,
    BaseBandStructure &bandStructure_)
    : Observable(statisticsSweep_, crystal_), bandStructure(bandStructure_) {
  tensordxd =
      Eigen::Tensor<double, 3>(numCalcs, dimensionality, dimensionality);
  tensordxd.setZero();
}

// copy constructor
PhononThermalConductivity::PhononThermalConductivity(
    const PhononThermalConductivity &that)
    : Observable(that), bandStructure(that.bandStructure) {}

// copy assigmnent
PhononThermalConductivity &PhononThermalConductivity::operator=(
    const PhononThermalConductivity &that) {
  Observable::operator=(that);
  if (this != &that) {
    bandStructure = that.bandStructure;
  }
  return *this;
}

PhononThermalConductivity PhononThermalConductivity::operator-(
    const PhononThermalConductivity &that) {
  PhononThermalConductivity newObservable(statisticsSweep, crystal,
                                          bandStructure);
  baseOperatorMinus(newObservable, that);
  return newObservable;
}

void PhononThermalConductivity::calcFromCanonicalPopulation(VectorBTE &f) {
  VectorBTE n = f;
  n.canonical2Population();  // n = bose (bose+1) f
  calcFromPopulation(n);
}

void PhononThermalConductivity::calcFromPopulation(VectorBTE &n) {
  double norm = 1. / bandStructure.getNumPoints(true) /
                crystal.getVolumeUnitCell(dimensionality);


  auto excludeIndeces = n.excludeIndeces;

  tensordxd.setZero();

  #pragma omp parallel
  {
    // we do manually the reduction, to avoid custom type declaration
    // which is not always allowed by the compiler e.g. by clang

    // first omp parallel for on a private variable
    Eigen::Tensor<double,3> tensorPrivate(numCalcs,dimensionality,dimensionality);
    tensorPrivate.setZero();

    #pragma omp for nowait
    for ( long is : bandStructure.parallelStateIterator() ) {
      double en = bandStructure.getEnergy(is);
      Eigen::Vector3d vel = bandStructure.getGroupVelocity(is);

      // skip the acoustic phonons
      if (std::find(excludeIndeces.begin(), excludeIndeces.end(), is) !=
          excludeIndeces.end())
        continue;

      for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
        for (long i = 0; i < dimensionality; i++) {
          for (long j = 0; j < dimensionality; j++) {
            tensorPrivate(iCalc, i, j) += n(iCalc, i, is) * vel(j) * en * norm;
          }
        }
      }
    }

    // now we do the reduction thread by thread
    #pragma omp critical
    {
      for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
        for (long i = 0; i < dimensionality; i++) {
          for (long j = 0; j < dimensionality; j++) {
            tensordxd(iCalc, i, j) += tensorPrivate(iCalc, i, j);
          }
        }
      }
    }
  }
  // lastly, the states were distributed with MPI
  mpi->allReduceSum(&tensordxd);
}

void PhononThermalConductivity::calcVariational(VectorBTE &af, VectorBTE &f,
                                                VectorBTE &scalingCG) {
  double norm = 1. / bandStructure.getNumPoints(true) /
                crystal.getVolumeUnitCell(dimensionality);

  auto fUnscaled = f;
  fUnscaled = fUnscaled / scalingCG;

  calcFromCanonicalPopulation(fUnscaled);

  tensordxd *= tensordxd.constant(2.);

  Eigen::Tensor<double,3> tmpTensor = tensordxd.constant(0.); // retains shape
  #pragma omp parallel
  {
    Eigen::Tensor<double,3> tmpTensorPrivate = tensordxd.constant(0.); // retains shape
    #pragma omp for nowait
    for ( long is : bandStructure.parallelStateIterator() ) {
      for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
        for (long i = 0; i < dimensionality; i++) {
          for (long j = 0; j < dimensionality; j++) {
            tmpTensorPrivate(iCalc, i, j) +=
                f(iCalc, i, is) * af(iCalc, j, is) * norm;
          }
        }
      }
    }
    #pragma omp critical
    for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
      for (long i = 0; i < dimensionality; i++) {
        for (long j = 0; j < dimensionality; j++) {
          tmpTensor(iCalc, i, j) += tmpTensorPrivate(iCalc, i, j);
        }
      }
    }
  }
  mpi->allReduceSum(&tmpTensor);
  tensordxd -= tmpTensor;
}

void PhononThermalConductivity::calcFromRelaxons(SpecificHeat &specificHeat,
                                                 VectorBTE &relaxonV,
                                                 VectorBTE &relaxationTimes) {
  // we decide to skip relaxon states
  // 1) there is a relaxon with zero (or epsilon) eigenvalue -> infinite tau
  // 2) if we include (3) acoustic modes at gamma, we have 3 zero eigenvalues
  //    because we set some matrix rows/cols to zero
  long firstState = 1;
  firstState += relaxationTimes.excludeIndeces.size();

  tensordxd.setZero();

  #pragma omp parallel
  {
    Eigen::Tensor<double,3> tmpTensor = tensordxd.constant(0.);

    #pragma omp for nowait
    for (long is = firstState; is < relaxationTimes.numStates; is++) {
      for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
        auto tup = loc2Glob(iCalc);
        auto imu = std::get<0>(tup);
        auto it = std::get<1>(tup);
        double c = specificHeat.get(imu, it);

        if (relaxationTimes(iCalc, 0, is) <= 0.) continue;

        for (long i = 0; i < dimensionality; i++) {
          for (long j = 0; j < dimensionality; j++) {
            tmpTensor(iCalc, i, j) += c * relaxonV(iCalc, i, is) *
                                      relaxonV(iCalc, j, is) *
                                      relaxationTimes(iCalc, 0, is);
          }
        }
      }
    }
    #pragma omp critical
    for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
      for (long i = 0; i < dimensionality; i++) {
        for (long j = 0; j < dimensionality; j++) {
          tensordxd(iCalc, i, j) += tmpTensor(iCalc, i, j);
        }
      }
    }

  }
}

void PhononThermalConductivity::print() {
  if (!mpi->mpiHead()) return;

  std::string units;
  if (dimensionality == 1) {
    units = "W m / K";
  } else if (dimensionality == 2) {
    units = "W / K";
  } else {
    units = "W / m / K";
  }

  std::cout << "\n";
  std::cout << "Thermal Conductivity (" << units << ")\n";

  for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStat.temperature;

    std::cout << std::fixed;
    std::cout.precision(2);
    std::cout << "Temperature: " << temp * temperatureAuToSi << " (K)\n";
    std::cout.precision(5);
    for (long i = 0; i < dimensionality; i++) {
      std::cout << "  " << std::scientific;
      for (long j = 0; j < dimensionality; j++) {
        std::cout << " " << std::setw(13) << std::right;
        std::cout << tensordxd(iCalc, i, j) * thConductivityAuToSi;
      }
      std::cout << "\n";
    }
    std::cout << std::endl;
  }
}

void PhononThermalConductivity::outputToJSON(std::string outFileName) {

  if(mpi->mpiHead()) {
    std::string units;
    if (dimensionality == 1) {
      units = "W m / K";
    } else if (dimensionality == 2) {
      units = "W / K";
    } else {
      units = "W / m / K";
    }

    std::vector<double> temps;
    std::vector<std::vector<std::vector<double>>> conds;
    for (long iCalc = 0; iCalc < numCalcs; iCalc++) {

      // store temperatures
      auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
      double temp = calcStat.temperature;
      temps.push_back(temp*temperatureAuToSi);

      // store conductivity
      std::vector<std::vector<double>> rows;
      for (long i = 0; i < dimensionality; i++) {
        std::vector<double> cols;
        for (long j = 0; j < dimensionality; j++) {
          cols.push_back(tensordxd(iCalc, i, j) * thConductivityAuToSi);
        }
        rows.push_back(cols);
      }
      conds.push_back(rows);
    }

    // output to json
    nlohmann::json output;
    output["temperatures"] = temps;
    output["thermalConductivity"] = conds;
    output["temperatureUnit"] = "K";
    output["thermalConductivityUnits"] = units;
    std::ofstream o(outFileName);
    o << std::setw(3) << output << std::endl;
    o.close();
  }
}

void PhononThermalConductivity::print(const int &iter) {
  if (!mpi->mpiHead()) return;

  // get the time
  time_t currentTime;
  currentTime = time(NULL);
  // and format the time nicely
  char s[200];
  struct tm *p = localtime(&currentTime);
  strftime(s, 200, "%F, %T", p);

  std::cout << "Iteration: " << iter << " | " << s << "\n";
  for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
    auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
    double temp = calcStat.temperature;
    std::cout << std::fixed;
    std::cout.precision(2);
    std::cout << "T = " << temp * temperatureAuToSi << ", k = ";
    std::cout.precision(5);
    for (long i = 0; i < dimensionality; i++) {
      std::cout << std::scientific;
      std::cout << tensordxd(iCalc, i, i) * thConductivityAuToSi << " ";
    }
    std::cout << "\n";
  }
  std::cout << std::endl;
}

int PhononThermalConductivity::whichType() { return is2Tensor; }
