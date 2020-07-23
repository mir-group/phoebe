#include "vector_bte.h"
#include "constants.h"

// default constructor
VectorBTE::VectorBTE(StatisticsSweep &statisticsSweep_,
                     BaseBandStructure &bandStructure_,
                     const long &dimensionality_)
    : BaseVectorBTE(statisticsSweep_, bandStructure_.getNumStates(),
                    dimensionality_),
      bandStructure(bandStructure_) {

  if (bandStructure.getParticle().isPhonon()) {
    for (long is = 0; is < numStates; is++) {
      double en = bandStructure.getEnergy(is);
      if (en < 0.1 / ryToCmm1) { // cutoff at 0.1 cm^-1
        excludeIndeces.push_back(is);
      }
    }
  }
}

// copy constructor
VectorBTE::VectorBTE(const VectorBTE &that)
    : BaseVectorBTE(that), bandStructure(that.bandStructure) {}

// copy assignment
VectorBTE &VectorBTE::operator=(const VectorBTE &that) {
  BaseVectorBTE::operator=(that);
  if (this != &that) {
    bandStructure = that.bandStructure;
  }
  return *this;
}

// product operator overload
Eigen::VectorXd VectorBTE::dot(const VectorBTE &that) {
  if (that.numCalcs != numCalcs || that.numStates != numStates) {
    Error e("The 2 VectorBTE must be aligned for dot() to work.");
  }
  Eigen::VectorXd result(numCalcs);
  result.setZero();
  for (long is : bandStructure.parallelStateIterator()) {
    for (long i = 0; i < numCalcs; i++) {
      result(i) += this->data(i, is) * that.data(i, is);
    }
  }
  mpi->allReduceSum(&result);
  return result;
}

VectorBTE VectorBTE::baseOperator(VectorBTE &that, const int &operatorType) {
  VectorBTE newPopulation(statisticsSweep, bandStructure, dimensionality);

  if (dimensionality == that.dimensionality) {

    if (operatorType == operatorSums) {
      newPopulation.data << this->data.array() + that.data.array();
    } else if (operatorType == operatorDivs) {
      newPopulation.data << this->data.array() / that.data.array();
    } else if (operatorType == operatorProd) {
      newPopulation.data << this->data.array() * that.data.array();
    } else if (operatorType == operatorDiff) {
      newPopulation.data << this->data.array() - that.data.array();
    } else {
      Error e("Operator type for VectorBTE not recognized");
    }

  } else if (that.dimensionality == 1) {

    for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
      auto tup = loc2Glob(iCalc);
      auto imu = std::get<0>(tup);
      auto it = std::get<1>(tup);
      auto idim = std::get<2>(tup);
      auto i2 = that.glob2Loc(imu, it, DimIndex(0));

      if (operatorType == operatorSums) {
        newPopulation.data.row(iCalc) =
            this->data.row(iCalc).array() + that.data.row(i2).array();
      } else if (operatorType == operatorDivs) {
        newPopulation.data.row(iCalc) =
            this->data.row(iCalc).array() / that.data.row(i2).array();
      } else if (operatorType == operatorProd) {
        newPopulation.data.row(iCalc) =
            this->data.row(iCalc).array() * that.data.row(i2).array();
      } else if (operatorType == operatorDiff) {
        newPopulation.data.row(iCalc) =
            this->data.row(iCalc).array() - that.data.row(i2).array();
      } else {
        Error e("Operator type for VectorBTE not recognized");
      }
    }
  } else {
    Error e("VectorBTE can't handle dimensionality for this case");
  }
  for (auto is : excludeIndeces) {
    newPopulation.data.col(is).setZero();
  }
  return newPopulation;
}

// product operator overload
VectorBTE VectorBTE::operator*(VectorBTE &that) {
  return baseOperator(that, operatorProd);
}

// product operator overload
VectorBTE VectorBTE::operator*(const double &scalar) {
  VectorBTE newPopulation(statisticsSweep, bandStructure, dimensionality);
  for (long i = 0; i < numCalcs; i++) {
    newPopulation.data.row(i) = this->data.row(i) * scalar;
  }
  return newPopulation;
}

// product operator overload
VectorBTE VectorBTE::operator*(const Eigen::VectorXd &vector) {
  VectorBTE newPopulation(statisticsSweep, bandStructure, dimensionality);
  for (long i = 0; i < numCalcs; i++) {
    newPopulation.data.row(i) = this->data.row(i) * vector(i);
  }
  return newPopulation;
}

// product operator overload
VectorBTE VectorBTE::operator*(ParallelMatrix<double> &matrix) {

  if (numCalcs != dimensionality) {
    // you'd need to keep in memory a lot of matrices.
    Error e("We didn't implement VectorBTE * matrix for numCalcs > 1");
  }
  if (matrix.rows() != numStates) {
    Error e("VectorBTE and Matrix not aligned");
  }
  VectorBTE newPopulation(statisticsSweep, bandStructure, dimensionality);
  newPopulation.data.setZero();
  for (auto ijtup : matrix.getAllLocalStates()) {
    auto i = std::get<0>(ijtup);
    auto j = std::get<1>(ijtup);
    for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
      newPopulation.data(iCalc, j) += data(iCalc, i) * matrix(i, j); // -5e-12
    }
  }
  mpi->allReduceSum(&newPopulation.data);
  return newPopulation;
}

// sum operator overload
VectorBTE VectorBTE::operator+(VectorBTE &that) {
  return baseOperator(that, operatorSums);
}

// product operator overload
VectorBTE VectorBTE::operator-(VectorBTE &that) {
  return baseOperator(that, operatorDiff);
}

// difference operator overload
VectorBTE VectorBTE::operator-() {
  VectorBTE newPopulation(statisticsSweep, bandStructure, dimensionality);
  newPopulation.data = -this->data;
  return newPopulation;
}

// division operator overload
VectorBTE VectorBTE::operator/(VectorBTE &that) {
  return baseOperator(that, operatorDivs);
}

VectorBTE VectorBTE::sqrt() {
  VectorBTE newPopulation(statisticsSweep, bandStructure, dimensionality);
  newPopulation.data << this->data.array().sqrt();
  return newPopulation;
}

VectorBTE VectorBTE::reciprocal() {
  VectorBTE newPopulation(statisticsSweep, bandStructure, dimensionality);
  newPopulation.data << 1. / this->data.array();
  return newPopulation;
}

void VectorBTE::canonical2Population() {
  auto particle = bandStructure.getParticle();
  for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
    auto tup = loc2Glob(iCalc);
    auto imu = std::get<0>(tup);
    auto it = std::get<1>(tup);
    auto idim = std::get<2>(tup);
    auto calcStatistics = statisticsSweep.getCalcStatistics(it, imu);
    auto temp = calcStatistics.temperature;
    auto chemPot = calcStatistics.chemicalPotential;
    for (long is = 0; is < numStates; is++) {
      double en = bandStructure.getEnergy(is);
      double term = particle.getPopPopPm1(en, temp, chemPot);
      data(iCalc, is) *= term;
    }
  }
}

void VectorBTE::population2Canonical() {
  auto particle = bandStructure.getParticle();
  for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
    auto tup = loc2Glob(iCalc);
    auto imu = std::get<0>(tup);
    auto it = std::get<1>(tup);
    auto idim = std::get<2>(tup);
    auto calcStatistics = statisticsSweep.getCalcStatistics(it, imu);
    auto temp = calcStatistics.temperature;
    auto chemPot = calcStatistics.chemicalPotential;
    for (long is = 0; is < numStates; is++) {
      double en = bandStructure.getEnergy(is);
      double term = particle.getPopPopPm1(en, temp, chemPot);
      data(iCalc, is) /= term;
    }
  }
}
