#include "bandstructure.h"

#include "exceptions.h"
#include "mpiHelper.h"
#include "particle.h"
#include "points.h"
#include "utilities.h"
#include "Matrix.h"

std::vector<long> BaseBandStructure::parallelStateIterator() {
    int numStates = getNumStates();
    return mpi->divideWorkIter(numStates);
}
//-----------------------------------------------------------------------------

FullBandStructure::FullBandStructure(long numBands_, Particle &particle_,
                                     bool withVelocities, bool withEigenvectors,
                                     Points &points_, bool isDistributed_)
    : particle(particle_), points(points_), isDistributed(isDistributed_) {

  numBands = numBands_;
  numAtoms = numBands_ / 3;
  numPoints = points.getNumPoints();
  hasVelocities = withVelocities;
  hasEigenvectors = withEigenvectors;

  // Initialize data structures depending on memory distribution.
  // If is distributed is true, numBlockCols is used to column/wavevector
  // distribute the internal matrices
  int numBlockCols = std::min((long)mpi->getSize(), numPoints);
  energies = Matrix<double>(numBands, numPoints, 1, numBlockCols, isDistributed);
  numLocalPoints = energies.localCols();
  if (hasVelocities) {
    velocities = Matrix<std::complex<double>>(
        numBands * numBands * 3, numPoints, 1, numBlockCols, isDistributed);
  }
  if (hasEigenvectors) {
    if ( particle.isPhonon() ) {
      eigenvectors = Matrix<std::complex<double>>(
          3 * numAtoms * numBands, numPoints, 1, numBlockCols, isDistributed);
    } else {
      eigenvectors = Matrix<std::complex<double>>(
          numBands * numBands, numPoints, 1, numBlockCols, isDistributed);
    }
  }
}

// copy constructor
FullBandStructure::FullBandStructure(const FullBandStructure &that)
    : particle(that.particle),
      points(that.points),
      isDistributed(that.isDistributed),
      energies(that.energies),
      velocities(that.velocities),
      eigenvectors(that.eigenvectors),
      numBands(that.numBands),
      numAtoms(that.numAtoms),
      numPoints(that.numPoints),
      numLocalPoints(that.numLocalPoints),
      hasEigenvectors(that.hasEigenvectors),
      hasVelocities(that.hasVelocities) {}

FullBandStructure &FullBandStructure::operator=(  // copy assignment
    const FullBandStructure &that) {
  if (this != &that) {
    particle = that.particle;
    points = that.points;
    isDistributed = that.isDistributed;
    energies = that.energies;
    velocities = that.velocities;
    eigenvectors = that.eigenvectors;
    numBands = that.numBands;
    numAtoms = that.numAtoms;
    numPoints = that.numPoints;
    numLocalPoints = that.numLocalPoints;
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
  if ( useFullGrid) {
    return points.getNumPoints();
  } else {
    return numLocalPoints;
  }
}

long FullBandStructure::getNumBands() { return numBands; }

long FullBandStructure::hasWindow() { return 0; }

bool FullBandStructure::getIsDistributed() { return isDistributed; }

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

// TODO this might be something I can simplify
std::vector<long> FullBandStructure::getWavevectorIndices() {
    std::vector<long> kptsList;
    // loop over local states
    for ( auto tup : energies.getAllLocalStates()) {
        // returns global indices for local index
        auto ik = std::get<1>(tup);
        kptsList.push_back(ik);
    }
    sort( kptsList.begin(), kptsList.end() );
    kptsList.erase( unique( kptsList.begin(), kptsList.end() ), kptsList.end() );
    return kptsList;
}

std::vector<std::tuple<WavevectorIndex,BandIndex>> FullBandStructure::getStateIndices() {
  auto allLocalStates = energies.getAllLocalStates();
  std::vector<std::tuple<WavevectorIndex, BandIndex>> idxs;
  for ( auto t : allLocalStates ) {
    auto ib = BandIndex(std::get<0>(t));
    auto ik = WavevectorIndex(std::get<1>(t));
    auto p = std::make_tuple(ik, ib);
    idxs.push_back(p);
  }
  return idxs;
}

std::vector<long> FullBandStructure::getBandIndices() {
  std::vector<long> bandsList;
  for(int ib = 0; ib < numBands; ib++) {
      bandsList.push_back(ib);
  }
  return bandsList;
}

const double &FullBandStructure::getEnergy(const long &stateIndex) {
  auto tup = decompress2Indeces(stateIndex, numPoints, numBands);
  auto ik = std::get<0>(tup);
  auto ib = std::get<1>(tup);
  if (!energies.indecesAreLocal(ib,ik)) {
    Error e("Cannot access a non-local energy.");
  }
  return energies(ib, ik);
}

const double &FullBandStructure::getEnergy(WavevectorIndex &ik, BandIndex &ib) {
  long ibb = ib.get();
  long ikk = ik.get();
  if (!energies.indecesAreLocal(ibb,ikk)) {
    Error e("Cannot access a non-local energy.");
  }
  return energies(ibb, ikk);
}

const double &FullBandStructure::getEnergy(StateIndex &is) {
  return getEnergy(is.get());
}

Eigen::VectorXd FullBandStructure::getEnergies(WavevectorIndex &ik) {
  Eigen::VectorXd x(numBands);
  if (!energies.indecesAreLocal(0,ik.get())) {
    Error e("Cannot access a non-local energy.");
  }
  for (int ib=0; ib<numBands; ib++) {
    x(ib) = energies(ib,ik.get());
  }
  return x;
}

Eigen::Vector3d FullBandStructure::getGroupVelocity(const long &stateIndex) {
  auto tup = decompress2Indeces(stateIndex, numPoints, numBands);
  auto ik = std::get<0>(tup);
  auto ib = std::get<1>(tup);
  if (!velocities.indecesAreLocal(ib,ik)) { // note ib is smaller than nRows
    Error e("Cannot access a non-local velocity.");
  }
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
  if (!velocities.indecesAreLocal(0,ikk)) {
    Error e("Cannot access a non-local velocity.");
  }
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
  if (!velocities.indecesAreLocal(0,ikk)) {
    Error e("Cannot access a non-local velocity.");
  }
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
  if (!eigenvectors.indecesAreLocal(0,ikk)) {
    Error e("Cannot access a non-local velocity.");
  }

  Eigen::MatrixXcd eigs(numBands, numBands);
  eigs.setZero();
  for (long ib1 = 0; ib1 < numBands; ib1++) {
    for (long ib2 = 0; ib2 < numBands; ib2++) {
      long ind = compress2Indeces(ib1, ib2, numBands, numBands);
      eigs(ib1, ib2) = eigenvectors(ind, ikk);
    }
  }
  return eigs;
}

Eigen::Tensor<std::complex<double>, 3> FullBandStructure::getPhEigenvectors(
    WavevectorIndex &ik) {
  long ikk = ik.get();
  if (!eigenvectors.indecesAreLocal(0,ikk)) {
    Error e("Cannot access a non-local velocity.");
  }
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
  auto tup = decompress2Indeces(stateIndex, numPoints, numBands);
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
  if (!energies.indecesAreLocal(0,ik)) {
    // col distributed, only need to check ik
    Error e("Cannot access a non-local energy");
  }
  for (int ib = 0; ib < energies.localRows(); ib++) {
    energies(ib, ik) = energies_(ib);
  }
}

void FullBandStructure::setEnergies(Point &point, Eigen::VectorXd &energies_) {
  long ik = point.getIndex();
  if (!energies.indecesAreLocal(0,ik)) {
    // col distributed, only need to check ik
    Error e("Cannot access a non-local energy");
  }
  for (int ib = 0; ib < energies.localRows(); ib++) {
    energies(ib, ik) = energies_(ib);
  }
}

void FullBandStructure::setVelocities(
    Point &point, Eigen::Tensor<std::complex<double>, 3> &velocities_) {
  if (!hasVelocities) {
    Error e("FullBandStructure was initialized without velocities");
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
  if (!velocities.indecesAreLocal(0,ik)) {
    // col distributed, only need to check ik
    Error e("Cannot access a non-local velocity");
  }
  for (int ib = 0; ib < velocities.localRows(); ib++) {
    velocities(ib, ik) = tmpVelocities_(ib);
  }
}

void FullBandStructure::setEigenvectors(Point &point,
                                        Eigen::MatrixXcd &eigenvectors_) {
  if (!hasEigenvectors) {
    Error e("FullBandStructure was initialized without eigvecs");
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
  if (!eigenvectors.indecesAreLocal(0,ik)) {
    // col distributed, only need to check ik
    Error e("Cannot access a non-local eigenvector.");
  }
  for (int ib = 0; ib < eigenvectors.localRows(); ib++) {
    eigenvectors(ib, ik) = tmp(ib);
  }
}

Eigen::VectorXd FullBandStructure::getBandEnergies(long &bandIndex) {
  // note: here we use the getWavevectorIndices function because if the
  // energies are distributed, we need to use global k indices
  // when calling energies(ib,ik)
  Eigen::VectorXd bandEnergies(energies.localCols());
  std::vector<long> wavevectorIndices = getWavevectorIndices();
  for (int i = 0; i < energies.localCols(); i++) {
    long ik = wavevectorIndices[i];  // global wavevector index
    bandEnergies(i) = energies(bandIndex, ik);
  }
  return bandEnergies;
}
