#include "bandstructure.h"
#include "points.h"
#include "state.h"
#include "particle.h"
#include "exceptions.h"
#include "utilities.h"
#include "mpiHelper.h"
#include "PMatrix.h"

Particle BaseBandStructure::getParticle() {
    Error e("BaseBandStructure method not implemented");
    return Particle(Particle::phonon);
}

Points BaseBandStructure::getPoints() {
    Error e("BaseBandStructure method not implemented");
}

Point BaseBandStructure::getPoint(const long &pointIndex) {
    Error e("BaseBandStructure method not implemented");
    (void) pointIndex;
}

long BaseBandStructure::getNumPoints(const bool &useFullGrid) {
    (void) useFullGrid;
    Error e("BaseBandStructure method not implemented");
    return 0;
}

long BaseBandStructure::getNumBands() {
    Error e("BaseBandStructure method not implemented");
    return 0;
}

long BaseBandStructure::hasWindow() {
    Error e("BaseBandStructure method not implemented");
    return 0;
}

State BaseBandStructure::getState(Point &point) {
    Error e("BaseBandStructure method not implemented");
    (void) point;
}

State BaseBandStructure::getState(const long &pointIndex) {
    Point point = getPoint(pointIndex);
    return getState(point);
}

long BaseBandStructure::getIndex(const WavevectorIndex &ik,
        const BandIndex &ib) {
    Error e("BaseBandStructure method not implemented");
    (void) ik;
    (void) ib;
    return 0;
}

std::tuple<WavevectorIndex,BandIndex> BaseBandStructure::getIndex(
        const long &is) {
    Error e("BaseBandStructure method not implemented");
    (void) is;
    auto ikk = WavevectorIndex(-1);
    auto ibb = BandIndex(-1);
    return {ikk, ibb};
}

long BaseBandStructure::getNumStates() {
    Error e("BaseBandStructure method not implemented");
    return 0;
}

std::vector<long> BaseBandStructure::getWavevectorIndices() {
    Error e("BaseBandStructure method not implemented");
    std::vector<long> wavevectorIndices; 
    return wavevectorIndices;
} 

std::vector<int> BaseBandStructure::parallelStateIterator() {
    int numStates = getNumStates();
    return mpi->divideWorkIter(numStates);
}

const double& BaseBandStructure::getEnergy(const long &stateIndex) {
    Error e("BaseBandStructure method not implemented");
    (void) stateIndex;
    return 0.;
}

Eigen::Vector3d BaseBandStructure::getGroupVelocity(const long &stateIndex) {
    Error e("BaseBandStructure method not implemented");
    (void) stateIndex;
    return Eigen::Vector3d::Zero();
}

Eigen::Vector3d BaseBandStructure::getWavevector(const long &stateIndex) {
    Error e("BaseBandStructure method not implemented");
    (void) stateIndex;
    return Eigen::Vector3d::Zero();
}

double BaseBandStructure::getWeight(const long &stateIndex) {
  (void) stateIndex;
  return 0.;
}

void BaseBandStructure::setEnergies(Point &point, Eigen::VectorXd &energies_) {
    Error e("BaseBandStructure method not implemented");
    (void) point;
    (void) energies_;
}

void BaseBandStructure::setVelocities(Point &point,
        Eigen::Tensor<std::complex<double>, 3> &velocities_) {
    Error e("BaseBandStructure method not implemented");
    (void) point;
    (void) velocities_;
}

void BaseBandStructure::setEigenvectors(Point &point,
        Eigen::MatrixXcd &eigenvectors_) {
    Error e("BaseBandStructure method not implemented");
    (void) point;
    (void) eigenvectors_;
}

//-----------------------------------------------------------------------------

FullBandStructure::FullBandStructure(long numBands_, Particle &particle_,
        bool hasVelocities_, bool hasEigenvectors_, Points &points_, bool isDistributed_) :
        particle(particle_), hasVelocities(hasVelocities_), hasEigenvectors(hasEigenvectors_), points(points_), isDistributed(isDistributed_)  {

    // I need to crash if I'm using active points

    numBands = numBands_;
    numAtoms = numBands_ / 3;
    numPoints = points.getNumPoints(); 

    // Initialize data structures depending on memory distribution
    if(isDistributed) { 
        int numBlockCols = std::min((long)mpi->getSize(), numPoints );
        energies =  ParallelMatrix<double>(numBands,numPoints,1,numBlockCols); 
        numLocalPoints = energies.localRows(); 

        if(hasVelocities) { 
          velocities = ParallelMatrix<std::complex<double>>(
                numBands * numBands * 3, numPoints,1,numBlockCols);
        }
        if(hasEigenvectors) {
          eigenvectors = ParallelMatrix<std::complex<double>>(
                3 * numAtoms * numBands, numPoints,1,numBlockCols);
        }
    }
    else {
        energies = Matrix<double>(numBands,numPoints);
        numLocalPoints = numPoints; 

        if(hasVelocities) {
          velocities = Matrix<std::complex<double>>(
              numBands * numBands * 3, numPoints);
        }
        if(hasEigenvectors) {
          eigenvectors = Matrix<std::complex<double>>(
                3 * numAtoms * numBands, numPoints);
        }
    }
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
FullBandStructure::FullBandStructure(const FullBandStructure &that) :
        particle(that.particle), points(that.points), energies(that.energies),
        velocities(that.velocities), eigenvectors(that.eigenvectors),
        rawEnergies(that.rawEnergies), rawVelocities(that.rawVelocities),
        rawEigenvectors(that.rawEigenvectors), energiesCols(that.energiesCols),
        velocitiesCols(that.velocitiesCols),
        eigenvectorsCols(that.eigenvectorsCols), numBands(that.numBands),
        numAtoms(that.numAtoms), hasEigenvectors(that.hasEigenvectors),
        hasVelocities(that.hasVelocities), isDistributed(that.isDistributed),
        numLocalPoints(that.numLocalPoints) {
}

FullBandStructure& FullBandStructure::operator =( // copy assignment
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
        numPoints = that.numPoints; 
        numLocalPoints = that.numLocalPoints; 
        hasEigenvectors = that.hasEigenvectors;
        hasVelocities = that.hasVelocities;
        isDistributed = that.isDistributed; 
    }
    return *this;
}

Particle FullBandStructure::getParticle() {
    return particle;
}

Point FullBandStructure::getPoint(const long &pointIndex) {
    // NOTE: points object is not distributed, so we *can* return this info, 
    // but it might be better to give an error if we access one which is 
    // not stored in this class in the distributed case
    return points.getPoint(pointIndex);
}

// NOTE: behavior has changed -- it now gives the
// local number of points in the distributed case. 
long FullBandStructure::getNumPoints(const bool &useFullGrid) {
    (void) useFullGrid;
    return numLocalPoints;
}

long FullBandStructure::getNumBands() {
    // If matrix is distributed over k cols, numBands is unaffected
    return numBands;
}

long FullBandStructure::getNumStates() {
    // because getNumPoints returns local numPoints, this applies in either case
    return numBands * getNumPoints();
}

// Returns a list of global wavevector indices corresponding 
// to the local wavevectors available to this process. 
std::vector<long> FullBandStructure::getWavevectorIndices() {
    std::vector<long> kptsList;
    // loop over local states
    for ( auto tup : energies.getAllLocalStates()) { 
        // returns global indices for local index
        auto ik = std::get<1>(tup);
        kptsList.push_back(ik);
    }
    return kptsList;
}

std::vector<long> FullBandStructure::getBandIndices() {
    std::vector<long> bandsList;
    // loop over local states
    for ( auto tup : energies.getAllLocalStates()) {
        auto ib = std::get<0>(tup);
        bandsList.push_back(ib);
    }
    return bandsList;
}

Points FullBandStructure::getPoints() {
   // NOTE: as there's no way to make a partial points object, 
   // this function returns all the points -- not just the local ones. 
   // If we want to return only local points, we should make a function to 
   // return a vector of only the local points -- a getLocalPoints method. 
   return points;
}

long FullBandStructure::hasWindow() {
    return 0;
}

State FullBandStructure::getState(Point &point) {
    long pointIndex = point.getIndex();
    return getState(pointIndex);
}

State FullBandStructure::getState(const long &pointIndex) {
    Point point = getPoint(pointIndex);
    // we construct the vector by defining begin() and end()
    // note that rawEnergies points at the start of the matrix
    // and pointIndex*energiesCols skips the first pointIndex-1 wavevectors
    // we are assuming column-major ordering!
    // thisEn points at the location with where the next numBands elements have
    // the desired energies for this wavevector.
    double *thisEn;
    thisEn = rawEnergies + pointIndex * energiesCols;

    // note: in some cases these are nullptr
    std::complex<double> *thisVel = nullptr;
    std::complex<double> *thisEig = nullptr;

    if (hasVelocities) {
        thisVel = rawVelocities + pointIndex * velocitiesCols;
    }
    if (hasEigenvectors) {
        thisEig = rawEigenvectors + pointIndex * eigenvectorsCols;
    }

    State s(point, thisEn, numBands, numBands, thisVel, thisEig);
    return s;
}

long FullBandStructure::getIndex(const WavevectorIndex &ik,
        const BandIndex &ib) {
    return ik.get() * numBands + ib.get();
}

std::tuple<WavevectorIndex,BandIndex> FullBandStructure::getIndex(
        const long &is) {
    int ik = is / numBands;
    int ib = is - ik * numBands;
    auto ikk = WavevectorIndex(ik);
    auto ibb = BandIndex(ib);
    return {ikk, ibb};
}

const double& FullBandStructure::getEnergy(const long &stateIndex) {
    auto tup = decompress2Indeces(stateIndex,numPoints,numBands);
    auto ik = std::get<0>(tup);
    auto ib = std::get<1>(tup);
    return energies(ib, ik);
}

Eigen::Vector3d FullBandStructure::getGroupVelocity(const long &stateIndex) {
    auto tup = decompress2Indeces(stateIndex,numPoints,numBands);
    auto ik = std::get<0>(tup);
    auto ib = std::get<1>(tup);
    Eigen::Vector3d vel;
    for (long i = 0; i < 3; i++) {
        long ind = compress3Indeces(ib, ib, i, numBands, numBands, 3);
        vel(i) = velocities(ind, ik).real();
    }
    return vel;
}

Eigen::Vector3d FullBandStructure::getWavevector(const long &stateIndex) {
    auto tup = decompress2Indeces(stateIndex,numPoints,numBands);
    auto ik = std::get<0>(tup);
    return points.getPoint(ik).getCoords(Points::cartesianCoords,true);
}

double FullBandStructure::getWeight(const long &stateIndex) {
    auto tup = getIndex(stateIndex);
    auto ik = std::get<0>(tup);
    return points.getWeight(ik.get());
}

Eigen::VectorXd FullBandStructure::getBandEnergies(long &bandIndex) { 
    Eigen::VectorXd bandEnergies(energies.cols()); 
    for(int ik=0;ik<energies.localCols();ik++) bandEnergies(ik) =  energies(bandIndex,ik);
    return bandEnergies;
}

void FullBandStructure::setEnergies(Eigen::Vector3d &coords, Eigen::VectorXd &energies_) {
    long ik = points.getIndex(coords); // global kidx
    for(int ib=0;ib<energies.localRows();ib++) energies(ib,ik) = energies_(ib);
}

void FullBandStructure::setEnergies(Point &point, Eigen::VectorXd &energies_) {
    long ik = point.getIndex();
    for(int ib=0;ib<energies.localRows();ib++) energies(ib,ik) = energies_(ib);
}

void FullBandStructure::setVelocities(Point &point,
        Eigen::Tensor<std::complex<double>, 3> &velocities_) {
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
    for(int ib=0;ib<velocities.localRows();ib++) velocities(ib,ik) = tmpVelocities_(ib);
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
    for(int ib=0;ib<eigenvectors.localRows();ib++) eigenvectors(ib,ik) = tmp(ib);
}
