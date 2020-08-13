#include "electron_h0_fourier.h"

#include <cmath>
#include <string>

#include "constants.h"
#include "exceptions.h"
#include "io.h"
#include "particle.h"

ElectronH0Fourier::ElectronH0Fourier(Crystal &crystal_,
                                     FullPoints coarsePoints_,
                                     FullBandStructure coarseBandStructure_,
                                     double cutoff_)
    : crystal(crystal_),
      coarseBandStructure(coarseBandStructure_),
      coarsePoints(coarsePoints_),
      particle(Particle::electron) {
  numBands = coarseBandStructure.getNumBands();
  cutoff = cutoff_;
  numDataPoints = coarseBandStructure.getNumPoints();
  refWavevector =
      coarseBandStructure.getPoint(0).getCoords(Points::cartesianCoords);

  // to do the interpolation, we need the lattice vector basis:
  setPositionVectors();
  // now we look for the expansion coefficients that interpolates the bands
  // note that setPositionVectors must stay above this call
  Eigen::MatrixXcd expansionCoefficients_(numBands, numPositionVectors);
  Eigen::VectorXd energies(numDataPoints);
  expansionCoefficients_.setZero();
  energies.setZero();

  LoopPrint loopPrint("setting up Fourier interpolation", "bands", numBands);
  for (long iBand = 0; iBand < numBands; iBand++) {
    loopPrint.update();
    energies = coarseBandStructure.getBandEnergies(iBand);
    expansionCoefficients_.row(iBand) = getCoefficients(energies);
  }
  loopPrint.close();
  expansionCoefficients = expansionCoefficients_;
}

// Copy constructor
ElectronH0Fourier::ElectronH0Fourier(const ElectronH0Fourier &that)
    : crystal(that.crystal),
      coarseBandStructure(that.coarseBandStructure),
      coarsePoints(that.coarsePoints),
      particle(that.particle),
      expansionCoefficients(that.expansionCoefficients),
      numBands(that.numBands),
      cutoff(that.cutoff),
      numDataPoints(that.numDataPoints),
      numPositionVectors(that.numPositionVectors),
      minDistance(that.minDistance),
      positionDegeneracies(that.positionDegeneracies),
      positionVectors(that.positionVectors),
      refWavevector(that.refWavevector) {}

// Copy assignment
ElectronH0Fourier &ElectronH0Fourier::operator=(const ElectronH0Fourier &that) {
  if (this != &that) {
    crystal = that.crystal;
    coarseBandStructure = that.coarseBandStructure;
    coarsePoints = that.coarsePoints;
    particle = that.particle;
    expansionCoefficients = that.expansionCoefficients;
    numBands = that.numBands;
    cutoff = that.cutoff;
    numDataPoints = that.numDataPoints;
    numPositionVectors = that.numPositionVectors;
    minDistance = that.minDistance;
    positionDegeneracies = that.positionDegeneracies;
    positionVectors = that.positionVectors;
    refWavevector = that.refWavevector;
  }
  return *this;
}

Particle ElectronH0Fourier::getParticle() { return particle; }

long ElectronH0Fourier::getNumBands() { return numBands; }

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> ElectronH0Fourier::diagonalize(
    Point &point) {
  Eigen::Vector3d coords = point.getCoords(Points::cartesianCoords);
  auto tup = diagonalizeFromCoords(coords);
  auto energies = std::get<0>(tup);
  auto eigvecs = std::get<1>(tup);
  return {energies, eigvecs};
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
ElectronH0Fourier::diagonalizeFromCoords(Eigen::Vector3d &wavevector) {
  Eigen::MatrixXcd eigvecs(1, 1);
  eigvecs.setZero();
  Eigen::VectorXd energies(numBands);
  for (long ib = 0; ib < numBands; ib++) {
    energies(ib) = getEnergyFromCoords(wavevector, ib);
  }
  return {energies, eigvecs};
}

Eigen::Tensor<std::complex<double>, 3> ElectronH0Fourier::diagonalizeVelocity(
    Point &point) {
  Eigen::Vector3d coords = point.getCoords(Points::cartesianCoords);
  return diagonalizeVelocityFromCoords(coords);
}

Eigen::Tensor<std::complex<double>, 3>
ElectronH0Fourier::diagonalizeVelocityFromCoords(Eigen::Vector3d &coords) {
  Eigen::Tensor<std::complex<double>, 3> velocity(numBands, numBands, 3);
  velocity.setZero();
  for (long ib = 0; ib < numBands; ib++) {
    Eigen::Vector3d v = getGroupVelocityFromCoords(coords, ib);
    for (long i = 0; i < 3; i++) {
      velocity(ib, ib, i) = v(i);
    }
  }
  return velocity;
}

FullBandStructure ElectronH0Fourier::populate(Points &fullPoints,
                                              bool &withVelocities,
                                              bool &withEigenvectors,
                                              bool isDistributed) {
  FullBandStructure fullBandStructure(numBands, particle, withVelocities,
                                      withEigenvectors, fullPoints, isDistributed);

  for (auto ik : fullBandStructure.getWavevectorIndices()) {
    Point point = fullBandStructure.getPoint(ik);
    auto tup = diagonalize(point);
    auto ens = std::get<0>(tup);
    auto eigvecs = std::get<1>(tup);
    fullBandStructure.setEnergies(point, ens);
    if (withVelocities) {
      auto vels = diagonalizeVelocity(point);
      fullBandStructure.setVelocities(point, vels);
    }
  }
  return fullBandStructure;
}

//////////////////////////////////////////////////////////////////////////

double ElectronH0Fourier::getRoughnessFunction(Eigen::Vector3d position) {
  double norm = position.norm();
  return pow(1. - coeff1 * norm / minDistance, 2) +
         coeff2 * pow(norm / minDistance, 6);
}

std::complex<double> ElectronH0Fourier::getStarFunction(
    Eigen::Vector3d &wavevector, long &iR) {
  std::complex<double> phase =
      complexI * wavevector.dot(positionVectors.col(iR));
  std::complex<double> starFunction = exp(phase) / positionDegeneracies(iR);
  return starFunction;
}

Eigen::Vector3cd ElectronH0Fourier::getDerivativeStarFunction(
    Eigen::Vector3d &wavevector, long &iR) {
  std::complex<double> phase =
      complexI * wavevector.dot(positionVectors.col(iR));
  Eigen::Vector3cd starFunctionDerivative = complexI * positionVectors.col(iR) *
                                            exp(phase) /
                                            positionDegeneracies(iR);
  return starFunctionDerivative;
}

void ElectronH0Fourier::setPositionVectors() {
  // load the info on the supercell that we will use for the interpolation
  Eigen::Matrix3d directUnitCell = crystal.getDirectUnitCell();
  auto tup = coarsePoints.getMesh();
  auto grid = std::get<0>(tup);
  auto offset = std::get<1>(tup);

  // the cutoff specifies the grid on which lattice vectors are searched.
  // Given a lattice vector in integer coordinates R = (n0,n1,n2) in
  // integer coordinates (i.e. can be brought in cartesian coordinates as
  // R=n0*a0+n1*a1+n2*a2 where a* are the lattice vectors of the primitive
  // unit cell)
  // then, the coordinates of the lattice vector R vary over a range between
  // -Nk*cutoff and Nk*cutoff, where Nk is the coarse grid of wavevectors.
  // so, with this definition, cutoff is a positive integer.
  long searchSize0 = long(cutoff);
  long searchSize1 = long(cutoff);
  long searchSize2 = long(cutoff);
  // cutoff is an integer
  if (cutoff < 1.) {
    Error e("Fourier cutoff is too small: set it >=1");
  }

  // distDim is the size of the
  long distDim = 1;
  distDim *= ((searchSize0 + 1) * 2 + 1);
  distDim *= ((searchSize1 + 1) * 2 + 1);
  distDim *= ((searchSize2 + 1) * 2 + 1);

  std::vector<Eigen::Vector3d> tmpVectors;
  std::vector<double> tmpDegeneracies;
  long tmpNumPoints = 0;

  // what are we doing:
  // the "n" specifies a grid of lattice vectors, that are (2*(cutoff+1)+1)^3
  // times bigger than the grid of wavevectors. a grid of lattice vectors
  // equal to the grid of wavevectors would be the bare minimum for the
  // interpolation to work, and wouldn't be enough for anything good.

  // now, we loop over the vectors of supercell A.
  for (long n0 = -searchSize0 * grid(0); n0 <= searchSize0 * grid(0); n0++) {
    for (long n1 = -searchSize1 * grid(1); n1 <= searchSize1 * grid(1); n1++) {
      for (long n2 = -searchSize2 * grid(2); n2 <= searchSize2 * grid(2);
           n2++) {
        // loop over a supercell B of size ((searchSize+1)*2+1)^3
        // bigger than supercell A. We compute the distance between any
        // vector in supercell B w.r.t. the vector "n"

        std::vector<double> distances;
        for (long i0 = -searchSize0 - 1; i0 <= searchSize0 + 1; i0++) {
          for (long i1 = -searchSize1 - 1; i1 <= searchSize1 + 1; i1++) {
            for (long i2 = -searchSize2 - 1; i2 <= searchSize2 + 1; i2++) {
              Eigen::Vector3d dist;
              dist(0) = n0 - i0 * grid(0);
              dist(1) = n1 - i1 * grid(1);
              dist(2) = n2 - i2 * grid(2);
              // distances in cartesian space
              dist = directUnitCell * dist;
              double dist2 = dist.norm();
              distances.push_back(dist2);
            }
          }
        }

        // find the minimum distance out of all distances
        double distMin = distances[0];
        for (auto dist : distances) {
          if (dist < distMin) {
            distMin = dist;
          }
        }

        // the point at distDim/2 is the reference "n" vector
        // i.e. the one generated by i0=i1=i2=0
        // if it's the minimum vector, than it's in the
        // Wigner Seitz zone of supercell A.
        // Therefore we save this vector.
        if (abs(distances[distDim / 2] - distMin) < 1.0e-6) {
          tmpNumPoints += 1;

          // count its degeneracy
          double degeneracy = 0.;
          for (long i = 0; i < distDim; i++) {
            if (abs(distances[i] - distMin) < 1.0e-6) {
              degeneracy += 1.;
            }
          }
          tmpDegeneracies.push_back(degeneracy);

          Eigen::Vector3d thisVec;
          thisVec << n0, n1, n2;
          tmpVectors.push_back(thisVec);
        }
      }
    }
  }

  // now we store the list of these lattice vectors in the class members
  numPositionVectors = tmpNumPoints;
  positionDegeneracies = Eigen::VectorXd::Zero(numPositionVectors);
  positionVectors = Eigen::MatrixXd::Zero(3, numPositionVectors);
  long originIndex = 0;  // to look for R=0 vector
  for (long iR = 0; iR < numPositionVectors; iR++) {
    auto thisVec = tmpVectors[iR];
    // we convert from crystal to cartesian coordinates
    positionVectors.col(iR) = directUnitCell * thisVec;
    positionDegeneracies(iR) = tmpDegeneracies[iR];
    //
    if (thisVec.norm() < 1.0e-6) {
      originIndex = iR;
    }
  }

  // It's essential we keep the origin vectors at the index iR = 0
  // so we swap it. The other vectors don't need to be in any order.
  double tmp1 = positionDegeneracies(originIndex);
  double tmp2 = positionDegeneracies(0);
  positionDegeneracies(0) = tmp1;
  positionDegeneracies(originIndex) = tmp2;
  Eigen::Vector3d tmpv1 = positionVectors.col(originIndex);
  Eigen::Vector3d tmpv2 = positionVectors.col(0);
  positionVectors.col(0) = tmpv1;
  positionVectors.col(originIndex) = tmpv2;

  // the interpolation schemes also requires to know the minimum norm
  // of the lattice vectors
  minDistance = directUnitCell.col(0).norm();         // this is a first guess
  for (long iR = 1; iR < numPositionVectors; iR++) {  // exclude 0 vector!
    double thisNorm = positionVectors.col(iR).norm();
    if (thisNorm < minDistance) {
      minDistance = thisNorm;
    }
  }
}

Eigen::VectorXcd ElectronH0Fourier::getLagrangeMultipliers(
    Eigen::VectorXd energies) {
  Eigen::VectorXcd multipliers(numDataPoints - 1);
  Eigen::VectorXcd deltaEnergies(numDataPoints - 1);
  Eigen::MatrixXcd h(numDataPoints - 1, numDataPoints - 1);
  h.setZero();

  for (long i = 1; i < numDataPoints; i++) {
    deltaEnergies(i - 1) = energies(i) - energies(0);
  }

  Eigen::Vector3d iWavevector;
  std::complex<double> smki, smkj;
  Eigen::MatrixXcd smk(numDataPoints, numPositionVectors);
  for (long iR = 0; iR < numPositionVectors; iR++) {
    for (long i = 0; i < numDataPoints; i++) {
      iWavevector =
          coarseBandStructure.getPoint(i).getCoords(Points::cartesianCoords);
      std::complex<double> smki = getStarFunction(iWavevector, iR);
      smk(i, iR) = smki;
    }
  }

  for (long iR = 0; iR < numPositionVectors; iR++) {
    Eigen::Vector3d position = positionVectors.col(iR);
    double rho = getRoughnessFunction(position);
    std::complex<double> smk0 = getStarFunction(refWavevector, iR);
    for (long i = 0; i < numDataPoints - 1; i++) {
      smki = smk(i + 1, iR);
      for (long j = 0; j < numDataPoints - 1; j++) {
        smkj = smk(j + 1, iR);
        h(i, j) += (smki - smk0) * std::conj(smkj - smk0) / rho;
      }
    }
  }

  // solve h*multipliers = deltaEnergies
  multipliers = h.ldlt().solve(deltaEnergies);
  return multipliers;
}

Eigen::VectorXcd ElectronH0Fourier::getCoefficients(Eigen::VectorXd energies) {
  Eigen::VectorXcd multipliers = getLagrangeMultipliers(energies);

  Eigen::VectorXcd coefficients(numPositionVectors);
  coefficients.setZero();
  Eigen::Vector3d wavevector;

  for (long m = 1; m < numPositionVectors; m++) {
    std::complex<double> smk0 = getStarFunction(refWavevector, m);
    for (long i = 1; i < numDataPoints; i++) {
      wavevector =
          coarseBandStructure.getPoint(i).getCoords(Points::cartesianCoords);
      coefficients(m) +=
          multipliers(i - 1) * std::conj(getStarFunction(wavevector, m) - smk0);
    }
    Eigen::Vector3d position = positionVectors.col(m);
    coefficients(m) /= getRoughnessFunction(position);
  }

  // special case for R=0
  coefficients(0) = energies(0);
  for (long m = 1; m < numPositionVectors; m++) {
    coefficients(0) -= coefficients(m) * getStarFunction(refWavevector, m);
  }

  return coefficients;
}

double ElectronH0Fourier::getEnergyFromCoords(Eigen::Vector3d &wavevector,
                                              long &bandIndex) {
  double energy = 0.;
  for (long m = 0; m < numPositionVectors; m++) {
    std::complex<double> c =
        expansionCoefficients(bandIndex, m) * getStarFunction(wavevector, m);
    energy += c.real();
  }
  return energy;
}

Eigen::Vector3d ElectronH0Fourier::getGroupVelocityFromCoords(
    Eigen::Vector3d &wavevector, long &bandIndex) {
  Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
  for (long m = 0; m < numPositionVectors; m++) {
    Eigen::Vector3cd c = expansionCoefficients(bandIndex, m) *
                         getDerivativeStarFunction(wavevector, m);
    velocity += c.real();
  }
  return velocity;
}
