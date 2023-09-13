#include "phonon_viscosity.h"
#include "constants.h"
#include "mpiHelper.h"
//#include "viscosity_io.h"
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>

PhononViscosity::PhononViscosity(Context &context_,
                                 StatisticsSweep &statisticsSweep_,
                                 Crystal &crystal_,
                                 BaseBandStructure &bandStructure_)
    : Observable(context_, statisticsSweep_, crystal_),
      bandStructure(bandStructure_) {

  tensordxdxdxd = Eigen::Tensor<double, 5>(numCalculations, dimensionality, dimensionality, dimensionality, dimensionality);
  tensordxdxdxd.setZero();
}

void PhononViscosity::calcRTA(VectorBTE &tau) {

  double norm = 1. / context.getQMesh().prod() /
                crystal.getVolumeUnitCell(dimensionality);

  //if(mpi->mpiHead()) std::cout << "volume " << crystal.getVolumeUnitCell(dimensionality) << std::endl;

  auto particle = bandStructure.getParticle();
  tensordxdxdxd.setZero();

  auto excludeIndices = tau.excludeIndices;

  std::vector<int> iss = bandStructure.parallelIrrStateIterator();
  int niss = iss.size();

  Kokkos::View<double*****, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> tensordxdxdxd_k(tensordxdxdxd.data(), numCalculations, dimensionality, dimensionality, dimensionality, dimensionality);
  Kokkos::Experimental::ScatterView<double*****, Kokkos::LayoutLeft, Kokkos::HostSpace> scatter_tensordxdxdxd(tensordxdxdxd_k);
  Kokkos::parallel_for("phonon_viscosity", Kokkos::RangePolicy<Kokkos::HostSpace::execution_space>(0, niss), [&] (int iis){
      auto tmpTensor = scatter_tensordxdxdxd.access();
      int is = iss[iis];
      auto isIdx = StateIndex(is);
      int iBte = bandStructure.stateToBte(isIdx).get();

      // skip the acoustic phonons
      if (std::find(excludeIndices.begin(), excludeIndices.end(), iBte) != excludeIndices.end())
        return;

      auto en = bandStructure.getEnergy(isIdx);
      auto velIrr = bandStructure.getGroupVelocity(isIdx);
      auto qIrr = bandStructure.getWavevector(isIdx);

      auto rotations = bandStructure.getRotationsStar(isIdx);
      for (const Eigen::Matrix3d& rotation : rotations) {

      Eigen::Vector3d q = rotation * qIrr;
      q = bandStructure.getPoints().foldToBz(q,Points::cartesianCoordinates);
      q = bandStructure.getPoints().bzToWs(q,Points::cartesianCoordinates);
      Eigen::Vector3d vel = rotation * velIrr;

      for (int iCalc = 0; iCalc < numCalculations; iCalc++) {

        auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
        double kBT = calcStat.temperature;
        double chemPot = 0; //calcStat.chemicalPotential; For phel this will cause a bug! Must be zero.
        double boseP1 = particle.getPopPopPm1(en, kBT, chemPot);

        for (int i = 0; i < dimensionality; i++) {
          for (int j = 0; j < dimensionality; j++) {
            for (int k = 0; k < dimensionality; k++) {
              for (int l = 0; l < dimensionality; l++) {
                tmpTensor(iCalc, i, j, k, l) +=
                  q(i) * vel(j) * q(k) * vel(l) * boseP1 * tau(iCalc, 0, iBte) / kBT * norm;
              }
            }
          }
        }
      }
    }
  });
  Kokkos::Experimental::contribute(tensordxdxdxd_k, scatter_tensordxdxdxd);

  /*
#pragma omp parallel default(none) shared(tensordxdxdxd,bandStructure,excludeIndices,numCalculations,statisticsSweep,particle,norm,tau)
  {
    Eigen::Tensor<double, 5> tmpTensor = tensordxdxdxd.constant(0.);

#pragma omp for nowait
    for (int is : bandStructure.parallelIrrStateIterator()) {

    }
#pragma omp critical
    for (int iCalc = 0; iCalc < numCalculations; iCalc++) {
      for (int i = 0; i < dimensionality; i++) {
        for (int j = 0; j < dimensionality; j++) {
          for (int k = 0; k < dimensionality; k++) {
            for (int l = 0; l < dimensionality; l++) {
              tensordxdxdxd(iCalc, i, j, k, l) += tmpTensor(iCalc, i, j, k, l);
            }
          }
        }
      }
    }
  }
  */
  mpi->allReduceSum(&tensordxdxdxd);
}

void PhononViscosity::calcFromRelaxons(Eigen::VectorXd &eigenvalues,
                                       ParallelMatrix<double> &eigenvectors) {

  // to simplify, here I do everything considering there is a single
  // temperature (due to memory constraints)
  if (numCalculations > 1) {
    Error("Developer error: Viscosity for relaxons only for 1 temperature.");
  }

  double volume = crystal.getVolumeUnitCell(dimensionality);
  int numStates = bandStructure.getNumStates();
  int numRelaxons = eigenvalues.size();
  auto particle = bandStructure.getParticle();

  Eigen::VectorXd A(dimensionality);
  A.setZero();

  int iCalc = 0;
  auto calcStat = statisticsSweep.getCalcStatistics(iCalc);
  double temp = calcStat.temperature;
  double chemPot = 0.; //calcStat.chemicalPotential; // will cause a bug if phel is added in

  // Code by Andrea, annotation by Jenny
  // Here we are calculating Eq. 9 from the PRX Simoncelli 2020
  //    mu_ijkl = (eta_ijkl + eta_ilkj)/2
  // where
  //   eta_ijkl =
  //    sqrt(A_i A_j) \sum_(alpha > 0) w_i,alpha^j * w_k,alpha^l * tau_alpha
  //
  // For this, we will need three quantities:
  //    A_i = specific momentum in direction i
  //        = (1/kBT*vol) \sum_nu n_nu (n_nu+1) (hbar q_i)^2
  //    w^j_i,alpha = velocity tensor
  //        = 1/vol * \sum_nu phi_nu^i v_nu^j theta_nu^alpha
  //    phi^i_nu = drift eigenvectors  (appendix eq A12)
  //             = zero eigenvectors linked to momentum conservation
  //             = sqrt(n(n+1)/kbT * A_i) * hbar * q_i
  //
  // And here the definitions are:
  //    alpha = relaxons eigenlabel index
  //    ijkl = cartesian direction indices
  //    nu  = phonon mode index
  //    n   = bose factor
  //    phi = special zero eigenvectors linked to momentum conservation
  //        = drift eigenvectors -- see appendix eq A12
  //    q   = wavevector
  //    theta = relaxons eigenvector
  //    tau = relaxons eigenvalues/relaxation times

  // calculate first A_i
  for (int is : bandStructure.parallelStateIterator()) {
    auto isIdx = StateIndex(is);
    auto en = bandStructure.getEnergy(isIdx);
    double boseP1 = particle.getPopPopPm1(en, temp, chemPot); // = n(n+1)
    auto q = bandStructure.getWavevector(isIdx);
    q = bandStructure.getPoints().foldToBz(q,Points::cartesianCoordinates);
    q = bandStructure.getPoints().bzToWs(q,Points::cartesianCoordinates);
    for (int iDim = 0; iDim < dimensionality; iDim++) {
      A(iDim) += boseP1 * q(iDim) * q(iDim);
    }
  }
  A /= temp * context.getQMesh().prod() * volume;
  mpi->allReduceSum(&A);

/*
if(mpi->mpiHead()) {
  std::cout << "num states " << bandStructure.irrStateIterator().size() << std::endl;
  for (int is : bandStructure.irrStateIterator()) {
    auto isIdx = StateIndex(is);
    auto en = bandStructure.getEnergy(isIdx);
    double boseP1 = particle.getPopPopPm1(en, temp, chemPot); // = n(n+1)
    auto q = bandStructure.getWavevector(isIdx);
    q = bandStructure.getPoints().foldToBz(q,Points::cartesianCoordinates);
    q = bandStructure.getPoints().bzToWs(q,Points::cartesianCoordinates);
    std::cout << "idx en q " << is << " " << en << " " << q.transpose() << std::endl;
  }
}
mpi->barrier();
if(mpi->mpiHead()) std::cout << "Ai " << A.transpose() << std::endl;
*/
  // then calculate the drift eigenvectors, phi (eq A12)
  VectorBTE driftEigenvector(statisticsSweep, bandStructure, 3);
  for (int is : bandStructure.parallelStateIterator()) {
    auto isIdx = StateIndex(is);
    auto en = bandStructure.getEnergy(isIdx);
    double boseP1 = particle.getPopPopPm1(en, temp, chemPot); // = n(n+1)
    auto q = bandStructure.getWavevector(isIdx);
    q = bandStructure.getPoints().foldToBz(q,Points::cartesianCoordinates);
    q = bandStructure.getPoints().bzToWs(q,Points::cartesianCoordinates);
    for (auto iDim : {0, 1, 2}) {
      if (A(iDim) != 0.) {
        driftEigenvector(0, iDim, is) = q(iDim) * sqrt(boseP1 / (temp * A(iDim)));
      }
    }
  }
  mpi->allReduceSum(&driftEigenvector.data);

  // calculate the first part of w^j_i,alpha
  Eigen::Tensor<double, 3> tmpDriftEigvecs(3, 3, numStates);
  tmpDriftEigvecs.setZero();
  for (int is : bandStructure.parallelStateIterator()) {
    auto isIdx = StateIndex(is);
    auto v = bandStructure.getGroupVelocity(isIdx);
    for (int i : {0, 1, 2}) {
      for (int j : {0, 1, 2}) {
        tmpDriftEigvecs(i, j, is) = driftEigenvector(0, j, is) * v(i);
      }
    }
  }
  mpi->allReduceSum(&tmpDriftEigvecs);

  // now we're calculating w
  Eigen::Tensor<double, 3> w(3, 3, numStates);
  w.setZero();
  for (auto i : {0, 1, 2}) {
    for (auto j : {0, 1, 2}) {
      // drift eigenvectors * v -- only have phonon state indices
      // and cartesian directions
      Eigen::VectorXd x(numStates);
      for (int is = 0; is < numStates; is++) {
        x(is) = tmpDriftEigvecs(i, j, is);
      }

      // w^j_i,alpha = sum_is1 phi*v*theta
      std::vector<double> x2(numRelaxons, 0.);
      for (auto tup : eigenvectors.getAllLocalStates()) {
        auto is1 = std::get<0>(tup);
        auto alpha = std::get<1>(tup);
        if(alpha >= numRelaxons) continue; // wasn't calculated
        x2[alpha] += x(is1) * eigenvectors(is1, alpha);
      }
      mpi->allReduceSum(&x2);

      // normalize by 1/(Nq*Volume)
      for (int ialpha = 0; ialpha < numRelaxons; ialpha++) {
        w(i, j, ialpha) = x2[ialpha] / ( sqrt(volume) * sqrt(context.getQMesh().prod()) );
      }
      // Andrea's note: in Eq. 9 of PRX, w is normalized by V*N_q
      // here however I normalize the eigenvectors differently:
      // \sum_state theta_s^2 = 1, instead of 1/VN_q \sum_state theta_s^2 = 1
    }
  }

  // Eq. 9, Simoncelli PRX (2020)
  tensordxdxdxd.setZero();
  // eigenvectors and values may only be calculated up to numRelaxons, < numStates
  std::vector<size_t> iss = mpi->divideWorkIter(numRelaxons);
  int niss = iss.size();

// TODO why would we do this rather than a kokkos parallel for?
  Kokkos::View<double*****, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> tensordxdxdxd_k(tensordxdxdxd.data(), numCalculations, dimensionality, dimensionality, dimensionality, dimensionality);
  Kokkos::Experimental::ScatterView<double*****, Kokkos::LayoutLeft, Kokkos::HostSpace> scatter_tensordxdxdxd(tensordxdxdxd_k);

  Kokkos::parallel_for("phonon_viscosity", Kokkos::RangePolicy<Kokkos::HostSpace::execution_space>(0, niss), [&] (int iis){

      auto tmpTensor = scatter_tensordxdxdxd.access();
      int ialpha = iss[iis];
      if (eigenvalues(ialpha) <= 0.) { // avoid division by zero
        return; // return is here continue for kokkos
      }
      // TODO this must be improved, we're here trying to discard the bose eigenvector contribution,
      // which has a divergent lifetime and should not be counted
      if(ialpha == 3) {
        return;
      }
      for (int i = 0; i < dimensionality; i++) {
        for (int j = 0; j < dimensionality; j++) {
          for (int k = 0; k < dimensionality; k++) {
            for (int l = 0; l < dimensionality; l++) {
            tmpTensor(iCalc, i, j, k, l) += 0.5 *
            (w(i, j, ialpha) * w(k, l, ialpha) + w(i, l, ialpha) * w(k, j, ialpha)) *
            A(i) * A(k) / eigenvalues(ialpha);
            }
          }
        }
      }
    });

  Kokkos::Experimental::contribute(tensordxdxdxd_k, scatter_tensordxdxdxd);
  mpi->allReduceSum(&tensordxdxdxd);
}

void PhononViscosity::print() {
  if (!mpi->mpiHead())
    return;

  std::string units;
  if (dimensionality == 1) {
    units = "Pa s / m^2";
  } else if (dimensionality == 2) {
    units = "Pa s / m";
  } else {
    units = "Pa s";
  }

  std::cout << "\n";
  std::cout << "Thermal Viscosity (" << units << ")\n";
  std::cout << "i, j, k, eta[i,j,k,0], eta[i,j,k,1], eta[i,j,k,2]\n";

  double conversion = pow(hBarSi, 2) // momentum is hBar q
                      / pow(distanceRyToSi, dimensionality) // volume conversion
                      / twoPi // because angular frequencies
                      * rydbergSi /
                      hBarSi       // conversion time (q^2 v^2 tau = [time])
                      / rydbergSi; // temperature conversion

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
                      << tensordxdxdxd(iCalc, i, j, k, l) * conversion;
          }
          std::cout << "\n";
        }
      }
    }
    std::cout << std::endl;
  }
}

void PhononViscosity::outputToJSON(const std::string& outFileName) {

  if (mpi->mpiHead()) {

    std::string units;
    if (dimensionality == 1) {
      units = "Pa s / m^2";
    } else if (dimensionality == 2) {
      units = "Pa s / m";
    } else {
      units = "Pa s";
    }

    // NOTE: conversion checked because brute force check produces the same answer
    //double conversion = 2.9421015697e13 * 2.4188843265857e-17;

    double conversion =
        pow(hBarSi, 2)                        // momentum is hBar q
        / pow(distanceRyToSi, dimensionality) // volume conversion
        / twoPi // because angular frequencies
        * rydbergSi / hBarSi // conversion time (q^2 v^2 tau = [time])
        / rydbergSi;         // temperature conversion

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
              ijkl.push_back(tensordxdxdxd(iCalc, i, j, k, l) * conversion);
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
    output["temperatures"] = temps;
    output["phononViscosity"] = viscosity;
    output["temperatureUnit"] = "K";
    output["phononViscosityUnit"] = units;
    output["particleType"] = "phonon";
    std::ofstream o(outFileName);
    o << std::setw(3) << output << std::endl;
    o.close();
  }
}

int PhononViscosity::whichType() { return is4Tensor; }
