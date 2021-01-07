#include "vector_bte.h"

// default constructor
BaseVectorBTE::BaseVectorBTE(StatisticsSweep &statisticsSweep_,
                             const int &numStates_,
                             const int &dimensionality_)
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
  data.resize(numCalcs, numStates);
  data.setZero();
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
  excludeIndices = that.excludeIndices;
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
    excludeIndices = that.excludeIndices;
  }
  return *this;
}

// product operator overload
Eigen::MatrixXd BaseVectorBTE::dot(const BaseVectorBTE &that) {
  if (that.numCalcs != numCalcs || that.numStates != numStates) {
    Error e("The 2 BaseVectorBTE must be aligned for dot() to work.");
  }
  if (that.dimensionality != 3 ) {
    Error("VectorBTE dot is implemented for 3D vectors only");
  }
  Eigen::VectorXd result(numCalcs,3);
  result.setZero();
  for (int is : mpi->divideWorkIter(numStates)) {
    for (int iCalc = 0; iCalc < statisticsSweep.getNumCalcs(); iCalc++) {
      for (int i : {0,1,2}) {
        result(iCalc, i) += operator()(iCalc, i, is) * that(iCalc, i, is);
      }
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
    for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
      auto tup = loc2Glob(iCalc);
      auto imu = std::get<0>(tup);
      auto it = std::get<1>(tup);
      auto i2 = that.glob2Loc(imu, it, CartIndex(0));

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
  for (auto is : excludeIndices) {
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
  for (int i = 0; i < numCalcs; i++) {
    newPopulation.data.row(i) = this->data.row(i) * scalar;
  }
  return newPopulation;
}

// product operator overload
BaseVectorBTE BaseVectorBTE::operator*(const Eigen::VectorXd &vector) {
  BaseVectorBTE newPopulation = *this;
  for (int i = 0; i < numCalcs; i++) {
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
  for (int iCalc = 0; iCalc < numCalcs; iCalc++) {
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

int BaseVectorBTE::glob2Loc(const ChemPotIndex &imu, const TempIndex &it,
                             const CartIndex &idim) {
  int i = compress3Indices(imu.get(), it.get(), idim.get(), numChemPots,
                            numTemps, dimensionality);
  return i;
}

std::tuple<ChemPotIndex, TempIndex, CartIndex> BaseVectorBTE::loc2Glob(
    const int &i) {
  auto tup = decompress3Indices(i, numChemPots, numTemps, dimensionality);
  auto imu = std::get<0>(tup);
  auto it = std::get<1>(tup);
  auto idim = std::get<2>(tup);
  return {ChemPotIndex(imu), TempIndex(it), CartIndex(idim)};
}

// get/set operator
double &BaseVectorBTE::operator()(const int iCalc, const int iDim,
                                  const int iState) {
  return data(iCalc * dimensionality + iDim, iState);
}

// const get/set operator
const double &BaseVectorBTE::operator()(const int iCalc, const int iDim,
                                        const int iState) const {
  return data(iCalc * dimensionality + iDim, iState);
}
