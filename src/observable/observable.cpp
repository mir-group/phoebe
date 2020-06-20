#include "observable.h"
#include "constants.h"
#include <cmath>

Observable::Observable(StatisticsSweep &statisticsSweep_, Crystal &crystal_)
    : statisticsSweep(statisticsSweep_), crystal(crystal_) {
  numCalcs = statisticsSweep.getNumCalcs();
  numChemPots = statisticsSweep.getNumChemicalPotentials();
  numTemps = statisticsSweep.getNumTemperatures();
  dimensionality = crystal.getDimensionality();
}

// copy constructor
Observable::Observable(const Observable &that)
    : statisticsSweep(that.statisticsSweep), crystal(that.crystal) {
  dimensionality = that.dimensionality;
  numChemPots = that.numChemPots;
  numTemps = that.numTemps;
  numCalcs = that.numCalcs;
  scalar = that.scalar;
  vectord = that.vectord;
  tensordxd = that.tensordxd;
  tensordxdxdxd = that.tensordxdxdxd;
}

// copy assigmnent
Observable &Observable::operator=(const Observable &that) {
  if (this != &that) {
    statisticsSweep = that.statisticsSweep;
    crystal = that.crystal;
    numChemPots = that.numChemPots;
    numTemps = that.numTemps;
    numCalcs = that.numCalcs;
    scalar = that.scalar;
    vectord = that.vectord;
    tensordxd = that.tensordxd;
    tensordxdxdxd = that.tensordxdxdxd;
  }
  return *this;
}

long Observable::glob2Loc(const ChemPotIndex &imu, const TempIndex &it) {
  return compress2Indeces(imu.get(), it.get(), numChemPots, numTemps);
}

std::tuple<ChemPotIndex, TempIndex> Observable::loc2Glob(const long &i) {
  auto tup = decompress2Indeces(i, numChemPots, numTemps);
  auto imu = std::get<0>(tup);
  auto it = std::get<1>(tup);
  return {ChemPotIndex(imu), TempIndex(it)};
}

Observable Observable::operator-(const Observable &that) {
  Observable newObservable(statisticsSweep, crystal);
  baseOperatorMinus(newObservable, that);
  return newObservable;
}

void Observable::baseOperatorMinus(Observable &newObservable,
                                   const Observable &that) {
  if (whichType() == isScalar) {
    for (long is = 0; is < numCalcs; is++) {
      newObservable.scalar(is) = scalar(is) - that.scalar(is);
    }
  } else if (whichType() == isVector) {
    for (long is = 0; is < numCalcs; is++) {
      for (int i = 0; i < dimensionality; i++) {
        newObservable.vectord(is, i) = vectord(is, i) - that.vectord(is, i);
      }
    }
  } else if (whichType() == is2Tensor) {
    for (long is = 0; is < numCalcs; is++) {
      for (int i = 0; i < dimensionality; i++) {
        for (int j = 0; j < dimensionality; j++) {
          newObservable.tensordxd(is, i, j) =
              tensordxd(is, i, j) - that.tensordxd(is, i, j);
        }
      }
    }
  } else if (whichType() == is4Tensor) {
    for (long is = 0; is < numCalcs; is++) {
      for (int i = 0; i < dimensionality; i++) {
        for (int j = 0; j < dimensionality; j++) {
          for (int k = 0; k < dimensionality; k++) {
            for (int l = 0; l < dimensionality; l++) {
              newObservable.tensordxdxdxd(is, i, j, k, l) =
                  tensordxdxdxd(is, i, j, k, l) -
                  that.tensordxdxdxd(is, i, j, k, l);
            }
          }
        }
      }
    }
  }
}

int Observable::whichType() { return isScalar; }

Eigen::VectorXd Observable::getNorm() {
  Eigen::VectorXd norm(numCalcs);
  norm.setZero();
  if (whichType() == isScalar) {
    for (long is = 0; is < numCalcs; is++) {
      norm(is) = abs(scalar(is));
    }
  } else if (whichType() == isVector) {
    for (long is = 0; is < numCalcs; is++) {
      for (int i = 0; i < dimensionality; i++) {
        norm(is) += vectord(is, i) * vectord(is, i);
      }
      norm(is) = sqrt(norm(is)) / double(dimensionality);
    }
  } else if (whichType() == is2Tensor) {
    for (long is = 0; is < numCalcs; is++) {
      for (int i = 0; i < dimensionality; i++) {
        for (int j = 0; j < dimensionality; j++) {
          norm(is) += tensordxd(is, i, j) * tensordxd(is, i, j);
        }
      }
      norm(is) = sqrt(norm(is)) / double(dimensionality * dimensionality);
    }
  } else if (whichType() == is4Tensor) {
    for (long is = 0; is < numCalcs; is++) {
      for (int i = 0; i < dimensionality; i++) {
        for (int j = 0; j < dimensionality; j++) {
          for (int k = 0; k < dimensionality; k++) {
            for (int l = 0; l < dimensionality; l++) {
              norm(is) +=
                  tensordxdxdxd(is, i, j, k, l) * tensordxdxdxd(is, i, j, k, l);
            }
          }
        }
      }
      norm(is) = sqrt(norm(is)) / double(pow(dimensionality, 4));
    }
  }
  return norm;
}
