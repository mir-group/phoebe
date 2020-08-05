#include "constants.h"
#include "vector_bte.h"

// default constructor
BaseVectorBTE::BaseVectorBTE(StatisticsSweep &statisticsSweep_,
                             const long &numStates_,
                             const long &dimensionality_)
    : statisticsSweep(statisticsSweep_) {
  if (dimensionality_ <= 0) {
    Error e("BaseVectorBTE doesn't accept <=0 dimensions");
  }
  if (numStates_ <= 0) {
    Error e("BaseVectorBTE doesn't accept <=0 number of states");
    numStates = 0;
  }

  dimensionality = dimensionality_;
  numStates = numStates_;

  numCalcs = statisticsSweep.getNumCalcs();
  numCalcs *= dimensionality;

  numChemPots = statisticsSweep.getNumChemicalPotentials();
  numTemps = statisticsSweep.getNumTemperatures();

  data = Eigen::MatrixXd::Zero(numCalcs, numStates);
}

// copy constructor
BaseVectorBTE::BaseVectorBTE(const BaseVectorBTE &that)
    : statisticsSweep(that.statisticsSweep) {
  numCalcs = that.numCalcs;
  numStates = that.numStates;
  numChemPots = that.numChemPots;
  numTemps = that.numTemps;
  dimensionality = that.dimensionality;
  data = that.data;
  excludeIndeces = that.excludeIndeces;
}

// copy assignment
BaseVectorBTE &BaseVectorBTE::operator=(const BaseVectorBTE &that) {
  if (this != &that) {
    statisticsSweep = that.statisticsSweep;
    numCalcs = that.numCalcs;
    numStates = that.numStates;
    numChemPots = that.numChemPots;
    numTemps = that.numTemps;
    dimensionality = that.dimensionality;
    data = that.data;
    excludeIndeces = that.excludeIndeces;
  }
  return *this;
}

// product operator overload
Eigen::VectorXd BaseVectorBTE::dot(const BaseVectorBTE &that) {
  if (that.numCalcs != numCalcs || that.numStates != numStates) {
    Error e("The 2 BaseVectorBTE must be aligned for dot() to work.");
  }
  Eigen::VectorXd result(numCalcs);
  result.setZero();
  for (long is : mpi->divideWorkIter(numStates)) {
    for (long i = 0; i < numCalcs; i++) {
      result(i) += this->data(i, is) * that.data(i, is);
    }
  }
  mpi->allReduceSum(&result);
  return result;
}

BaseVectorBTE BaseVectorBTE::baseOperator(BaseVectorBTE &that,
                                          const int &operatorType) {
  BaseVectorBTE newPopulation = *this;

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
      Error e("Operator type for BaseVectorBTE not recognized");
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
        Error e("Operator type for BaseVectorBTE not recognized");
      }
    }
  } else {
    Error e("BaseVectorBTE can't handle dimensionality for this case");
  }
  for (auto is : excludeIndeces) {
    newPopulation.data.col(is).setZero();
  }
  return newPopulation;
}

// product operator overload
BaseVectorBTE BaseVectorBTE::operator*(BaseVectorBTE &that) {
  return baseOperator(that, operatorProd);
}

// product operator overload
BaseVectorBTE BaseVectorBTE::operator*(const double &scalar) {
  BaseVectorBTE newPopulation = *this;
  for (long i = 0; i < numCalcs; i++) {
    newPopulation.data.row(i) = this->data.row(i) * scalar;
  }
  return newPopulation;
}

// product operator overload
BaseVectorBTE BaseVectorBTE::operator*(const Eigen::VectorXd &vector) {
  BaseVectorBTE newPopulation = *this;
  for (long i = 0; i < numCalcs; i++) {
    newPopulation.data.row(i) = this->data.row(i) * vector(i);
  }
  return newPopulation;
}

// product operator overload
BaseVectorBTE BaseVectorBTE::operator*(ParallelMatrix<double> &matrix) {
  if (numCalcs != dimensionality) {
    // I mean, you'd need to keep in memory a lot of matrices.
    Error e("We didn't implement BaseVectorBTE * matrix for numCalcs > 1");
  }
  if (matrix.rows() != numStates) {
    Error e("BaseVectorBTE and Matrix not aligned");
  }
  BaseVectorBTE newPopulation = *this;
  newPopulation.data.setZero();
  for (long iCalc = 0; iCalc < numCalcs; iCalc++) {
    for (auto tup : matrix.getAllLocalStates()) {
      auto i = std::get<0>(tup);
      auto j = std::get<1>(tup);
      newPopulation.data(iCalc, i) = matrix(i, j) * data(iCalc, j);
    }
  }
  mpi->allReduceSum(&newPopulation.data);
  return newPopulation;
}

// sum operator overload
BaseVectorBTE BaseVectorBTE::operator+(BaseVectorBTE &that) {
  return baseOperator(that, operatorSums);
}

// product operator overload
BaseVectorBTE BaseVectorBTE::operator-(BaseVectorBTE &that) {
  return baseOperator(that, operatorDiff);
}

// inversion operator overload
BaseVectorBTE BaseVectorBTE::operator-() {
  BaseVectorBTE newPopulation = *this;
  newPopulation.data = -this->data;
  return newPopulation;
}

// division operator overload
BaseVectorBTE BaseVectorBTE::operator/(BaseVectorBTE &that) {
  return baseOperator(that, operatorDivs);
}

BaseVectorBTE BaseVectorBTE::sqrt() {
  BaseVectorBTE newPopulation = *this;
  newPopulation.data << this->data.array().sqrt();
  return newPopulation;
}

BaseVectorBTE BaseVectorBTE::reciprocal() {
  BaseVectorBTE newPopulation = *this;
  newPopulation.data << 1. / this->data.array();
  return newPopulation;
}

void BaseVectorBTE::setConst(const double &constant) {
  data.setConstant(constant);
}

long BaseVectorBTE::glob2Loc(const ChemPotIndex &imu, const TempIndex &it,
                             const DimIndex &idim) {
  long i = compress3Indeces(imu.get(), it.get(), idim.get(), numChemPots,
                            numTemps, dimensionality);
  return i;
}

std::tuple<ChemPotIndex, TempIndex, DimIndex> BaseVectorBTE::loc2Glob(
    const long &i) {
  auto tup =      decompress3Indeces(i, numChemPots, numTemps, dimensionality);
  auto imu = std::get<0>(tup);
  auto it = std::get<1>(tup);
  auto idim = std::get<2>(tup);
  return {ChemPotIndex(imu), TempIndex(it), DimIndex(idim)};
}
