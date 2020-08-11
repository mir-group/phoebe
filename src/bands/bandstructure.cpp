#include "bandstructure.h"

#include "exceptions.h"
#include "mpiHelper.h"
#include "particle.h"
#include "points.h"
#include "utilities.h"

std::vector<int> BaseBandStructure::parallelStateIterator() {
  int numStates = getNumStates();
  return mpi->divideWorkIter(numStates);
}

//-----------------------------------------------------------------------------

FullBandStructure::FullBandStructure(long numBands_, Particle &particle_,
                                     bool withVelocities, bool withEigenvectors,
                                     Points &points_)
    : particle(particle_), points(points_) {
  //	I need to crash if I'm using active points

  numBands = numBands_;
  numAtoms = numBands_ / 3;

  if (withVelocities) {
    hasVelocities = true;
    velocities =
        Eigen::MatrixXcd::Zero(numBands * numBands * 3, getNumPoints());
  }

  if (withEigenvectors) {
    hasEigenvectors = true;
    eigenvectors =
        Eigen::MatrixXcd::Zero(3 * numAtoms * numBands, getNumPoints());
  }

  energies = Eigen::MatrixXd::Zero(numBands, getNumPoints());

  // now, I want to manipulate the Eigen matrices at lower level
  // I create this pointer to data, so I can move it around
  rawEnergies = energies.data();
  if (hasVelocities) {
    rawVelocities = velocities.data();
  }
  if (hasEigenvectors) {
    rawEigenvectors = eigenvectors.data();
  }

  energiesCols = numBands;
  velocitiesCols = numBands * numBands * 3;
  eigenvectorsCols = numBands * numAtoms * 3;
}

// copy constructor
FullBandStructure::FullBandStructure(const FullBandStructure &that)
    : particle(that.particle),
      points(that.points),
      energies(that.energies),
      velocities(that.velocities),
      eigenvectors(that.eigenvectors),
      rawEnergies(that.rawEnergies),
      rawVelocities(that.rawVelocities),
      rawEigenvectors(that.rawEigenvectors),
      energiesCols(that.energiesCols),
      velocitiesCols(that.velocitiesCols),
      eigenvectorsCols(that.eigenvectorsCols),
      numBands(that.numBands),
      numAtoms(that.numAtoms),
      hasEigenvectors(that.hasEigenvectors),
      hasVelocities(that.hasVelocities) {}

FullBandStructure &FullBandStructure::operator=(  // copy assignment
    const FullBandStructure &that) {
  if (this != &that) {
    particle = that.particle;
    points = that.points;
    energies = that.energies;
    velocities = that.velocities;
    eigenvectors = that.eigenvectors;
    rawEnergies = that.rawEnergies;
    rawVelocities = that.rawVelocities;
    rawEigenvectors = that.rawEigenvectors;
    energiesCols = that.energiesCols;
    velocitiesCols = that.velocitiesCols;
    eigenvectorsCols = that.eigenvectorsCols;
    numBands = that.numBands;
    numAtoms = that.numAtoms;
    hasEigenvectors = that.hasEigenvectors;
    hasVelocities = that.hasVelocities;
  }
  return *this;
}

Particle FullBandStructure::getParticle() { return particle; }

Points FullBandStructure::getPoints() { return points; }

Point FullBandStructure::getPoint(const long &pointIndex) {
  return points.getPoint(pointIndex);
}

long FullBandStructure::getNumPoints(const bool &useFullGrid) {
  (void)useFullGrid;
  return points.getNumPoints();
}

long FullBandStructure::getNumBands() { return numBands; }

long FullBandStructure::hasWindow() { return 0; }

long FullBandStructure::getIndex(const WavevectorIndex &ik,
                                 const BandIndex &ib) {
  return ik.get() * numBands + ib.get();
}

std::tuple<WavevectorIndex, BandIndex> FullBandStructure::getIndex(
    const long &is) {
  int ik = is / numBands;
  int ib = is - ik * numBands;
  auto ikk = WavevectorIndex(ik);
  auto ibb = BandIndex(ib);
  return {ikk, ibb};
}

std::tuple<WavevectorIndex, BandIndex> FullBandStructure::getIndex(
    StateIndex &is) {
  return getIndex(is.get());
}

long FullBandStructure::getNumStates() { return numBands * getNumPoints(); }

const double &FullBandStructure::getEnergy(const long &stateIndex) {
  auto tup = decompress2Indeces(stateIndex, getNumPoints(), numBands);
  auto ik = std::get<0>(tup);
  auto ib = std::get<1>(tup);
  return energies(ib, ik);
}

const double &FullBandStructure::getEnergy(StateIndex &is) {
  return getEnergy(is.get());
}

Eigen::VectorXd FullBandStructure::getEnergies(WavevectorIndex &ik) {
  Eigen::VectorXd x = energies.col(ik.get());
  return x;
}

Eigen::Vector3d FullBandStructure::getGroupVelocity(const long &stateIndex) {
  auto tup = decompress2Indeces(stateIndex, getNumPoints(), numBands);
  auto ik = std::get<0>(tup);
  auto ib = std::get<1>(tup);
  Eigen::Vector3d vel;
  for (int i : {0, 1, 2}) {
    long ind = compress3Indeces(ib, ib, i, numBands, numBands, 3);
    vel(i) = velocities(ind, ik).real();
  }
  return vel;
}

Eigen::Vector3d FullBandStructure::getGroupVelocity(StateIndex &is) {
  return getGroupVelocity(is.get());
}

Eigen::MatrixXd FullBandStructure::getGroupVelocities(WavevectorIndex &ik) {
  long ikk = ik.get();
  Eigen::MatrixXd vel(numBands,3);
  for (int ib=0; ib<numBands; ib++ ) {
    for (int i : {0, 1, 2}) {
      int ind = compress3Indeces(ib, ib, i, numBands, numBands, 3);
      vel(ib,i) = velocities(ind, ikk).real();
    }
  }
  return vel;
}

Eigen::Tensor<std::complex<double>, 3> FullBandStructure::getVelocities(
    WavevectorIndex &ik) {
  long ikk = ik.get();
  Eigen::Tensor<std::complex<double>, 3> vel(numBands, numBands, 3);
  for (int ib1 = 0; ib1 < numBands; ib1++) {
    for (int ib2 = 0; ib2 < numBands; ib2++) {
      for (int i : {0, 1, 2}) {
        int ind = compress3Indeces(ib1, ib2, i, numBands, numBands, 3);
        vel(ib1,ib2,i) = velocities(ind, ikk);
      }
    }
  }
  return vel;
}

Eigen::MatrixXcd FullBandStructure::getEigenvectors(WavevectorIndex &ik) {
  long ikk = ik.get();
  Eigen::MatrixXcd eigs(numBands, numBands);
  eigs.setZero();
  for (int ib1 = 0; ib1 < numBands; ib1++) {
    for (int ib2 = 0; ib2 < numBands; ib2++) {
      int ind = compress2Indeces(ib1, ib2, numBands, numBands);
      eigs(ib1, ib2) = eigenvectors(ind, ikk);
    }
  }
  return eigs;
}

Eigen::Tensor<std::complex<double>, 3> FullBandStructure::getPhEigenvectors(
    WavevectorIndex &ik) {
  long ikk = ik.get();
  Eigen::Tensor<std::complex<double>, 3> eigs_(3, numAtoms, numBands);
  for (int ib = 0; ib < numBands; ib++) {
    for (int ia = 0; ia < numAtoms; ia++) {
      for (int ic : {0, 1, 2}) {
        long ind = compress3Indeces(ia, ic, ib, numAtoms, 3, numBands);
        eigs_(ic, ia, ib) = eigenvectors(ind, ikk);
      }
    }
  }
  return eigs_;
}

Eigen::Vector3d FullBandStructure::getWavevector(const long &stateIndex) {
  auto tup = decompress2Indeces(stateIndex, getNumPoints(), numBands);
  auto ik = std::get<0>(tup);
  return points.getPoint(ik).getCoords(Points::cartesianCoords, true);
}

Eigen::Vector3d FullBandStructure::getWavevector(StateIndex &is) {
  return getWavevector(is.get());
}

Eigen::Vector3d FullBandStructure::getWavevector(WavevectorIndex &ik) {
  return points.getPoint(ik.get()).getCoords(Points::cartesianCoords, true);
}

double FullBandStructure::getWeight(const long &stateIndex) {
  auto tup = getIndex(stateIndex);
  auto ik = std::get<0>(tup);
  return getWeight(ik);
}

double FullBandStructure::getWeight(StateIndex &is) {
  return getWeight(is.get());
}

double FullBandStructure::getWeight(WavevectorIndex &ik) {
  return points.getWeight(ik.get());
}

void FullBandStructure::setEnergies(Eigen::Vector3d &coords,
                                    Eigen::VectorXd &energies_) {
  long ik = points.getIndex(coords);
  energies.col(ik) = energies_;
}

void FullBandStructure::setEnergies(Point &point, Eigen::VectorXd &energies_) {
  long ik = point.getIndex();
  energies.col(ik) = energies_;
}

void FullBandStructure::setVelocities(
    Point &point, Eigen::Tensor<std::complex<double>, 3> &velocities_) {
  if (!hasVelocities) {
    Error e("FullBandStructure was initialized without velocities", 1);
  }
  // we convert from a tensor to a vector (how it's stored in memory)
  Eigen::VectorXcd tmpVelocities_(numBands * numBands * 3);
  for (long i = 0; i < numBands; i++) {
    for (long j = 0; j < numBands; j++) {
      for (long k = 0; k < 3; k++) {
        // Note: State must know this order of index compression
        long idx = compress3Indeces(i, j, k, numBands, numBands, 3);
        tmpVelocities_(idx) = velocities_(i, j, k);
      }
    }
  }
  long ik = point.getIndex();
  velocities.col(ik) = tmpVelocities_;
}

void FullBandStructure::setEigenvectors(Point &point,
                                        Eigen::MatrixXcd &eigenvectors_) {
  if (!hasEigenvectors) {
    Error e("FullBandStructure was initialized without eigvecs", 1);
  }
  // we convert from a matrix to a vector (how it's stored in memory)
  Eigen::VectorXcd tmp(numBands * numBands);
  for (long i = 0; i < numBands; i++) {
    for (long j = 0; j < numBands; j++) {
      // Note: State must know this order of index compression
      long idx = compress2Indeces(i, j, numBands, numBands);
      tmp(idx) = eigenvectors_(i, j);
    }
  }
  long ik = point.getIndex();
  eigenvectors.col(ik) = tmp;
}

Eigen::VectorXd FullBandStructure::getBandEnergies(long &bandIndex) {
  Eigen::VectorXd bandEnergies = energies.row(bandIndex);
  return bandEnergies;
}
