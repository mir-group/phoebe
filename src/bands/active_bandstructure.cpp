#include "active_bandstructure.h"

#include <cstdlib>

#include "bandstructure.h"
#include "exceptions.h"
#include "mpiHelper.h"
#include "window.h"

ActiveBandStructure::ActiveBandStructure(Particle &particle_,
                                         ActivePoints &activePoints)
    : particle(particle_), activePoints(activePoints) {}

// copy constructor
ActiveBandStructure::ActiveBandStructure(const ActiveBandStructure &that)
    : particle(that.particle),
      activePoints(that.activePoints),
      energies(that.energies),
      velocities(that.velocities),
      eigenvectors(that.eigenvectors),
      hasEigenvectors(that.hasEigenvectors),
      numStates(that.numStates),
      numPoints(that.numPoints),
      numBands(that.numBands),
      numFullBands(that.numFullBands),
      windowMethod(that.windowMethod),
      auxBloch2Comb(that.auxBloch2Comb),
      cumulativeKbOffset(that.cumulativeKbOffset),
      cumulativeKbbOffset(that.cumulativeKbbOffset) {}

ActiveBandStructure &ActiveBandStructure::operator=(
    const ActiveBandStructure &that) {  // assignment operator
  if (this != &that) {
    particle = that.particle;
    activePoints = that.activePoints;
    energies = that.energies;
    velocities = that.velocities;
    eigenvectors = that.eigenvectors;
    hasEigenvectors = that.hasEigenvectors;
    numStates = that.numStates;
    numPoints = that.numPoints;
    numBands = that.numBands;
    numFullBands = that.numFullBands;
    windowMethod = that.windowMethod;
    auxBloch2Comb = that.auxBloch2Comb;
    cumulativeKbOffset = that.cumulativeKbOffset;
    cumulativeKbbOffset = that.cumulativeKbbOffset;
  }
  return *this;
}

ActiveBandStructure::ActiveBandStructure(const ActivePoints &activePoints_,
                                         HarmonicHamiltonian *h0,
                                         const bool &withEigenvectors,
                                         const bool &withVelocities)
    : particle(Particle(h0->getParticle().getParticleKind())),
      activePoints(activePoints_) {
  // TODO might want to check on this section to make sure this works long term?
  // not sure what activePoints will return here
  numPoints = activePoints.getNumPoints();
  numFullBands = h0->getNumBands();
  numBands = Eigen::VectorXi::Zero(numPoints);
  for (int ik = 0; ik < numPoints; ik++) {
    numBands(ik) = numFullBands;
  }
  numStates = numFullBands * numPoints;
  hasEigenvectors = withEigenvectors;

  energies.resize(numPoints * numFullBands, 0.);
  velocities.resize(numPoints * numFullBands * numFullBands * 3, complexZero);
  eigenvectors.resize(numPoints * numFullBands * numFullBands, complexZero);

  windowMethod = Window::nothing;

  buildIndeces();

  // now we can loop over the trimmed list of points
  for (long ik : mpi->divideWorkIter(numPoints)) {
    Point point = activePoints.getPoint(ik);
    auto tup = h0->diagonalize(point);
    auto theseEnergies = std::get<0>(tup);
    auto theseEigenvectors = std::get<1>(tup);
    setEnergies(point, theseEnergies);

    if (withEigenvectors) {
      setEigenvectors(point, theseEigenvectors);
    }

    if (withVelocities) {
      auto thisVelocity = h0->diagonalizeVelocity(point);
      setVelocities(point, thisVelocity);
    }
  }
  mpi->allReduceSum(&energies);
  mpi->allReduceSum(&velocities);
  mpi->allReduceSum(&eigenvectors);
}

Particle ActiveBandStructure::getParticle() { return particle; }

bool ActiveBandStructure::hasPoints() {
  //	if ( activePoints != nullptr ) {
  return true;
  //	} else {
  //		return false;
  //	}
}

Points ActiveBandStructure::getPoints() {
  if (!hasPoints()) {
    Error e("ActiveBandStructure hasn't been populated yet");
  }
  return activePoints;
}

Point ActiveBandStructure::getPoint(const long &pointIndex) {
  if (!hasPoints()) {
    Error e("ActiveBandStructure hasn't been populated yet");
  }
  return activePoints.getPoint(pointIndex);
}

long ActiveBandStructure::getNumPoints(const bool &useFullGrid) {
  if (!hasPoints()) {
    Error e("ActiveBandStructure hasn't been populated yet");
  }
  if (useFullGrid) {
    return activePoints.getParentPoints().getNumPoints();
  } else {  // default
    return numPoints;
  }
}

long ActiveBandStructure::getNumBands() {
  if (windowMethod == Window::nothing) {
    return numFullBands;
  } else {
    Error e("ActiveBandStructure doesn't have constant number of bands");
    return 0;
  }
}

long ActiveBandStructure::hasWindow() { return windowMethod; }

bool ActiveBandStructure::getIsDistributed() { return false; }

long ActiveBandStructure::getIndex(const WavevectorIndex &ik,
                                   const BandIndex &ib) {
  return bloch2Comb(ik.get(), ib.get());
}

std::tuple<WavevectorIndex, BandIndex> ActiveBandStructure::getIndex(
    const long &is) {
  auto tup = comb2Bloch(is);
  auto ik = std::get<0>(tup);
  auto ib = std::get<1>(tup);
  auto ikk = WavevectorIndex(ik);
  auto ibb = BandIndex(ib);
  return {ikk, ibb};
}

std::tuple<WavevectorIndex, BandIndex> ActiveBandStructure::getIndex(
    StateIndex &is) {
  long iss = is.get();
  return getIndex(iss);
}

long ActiveBandStructure::getNumStates() {
  return numStates;
}

const double &ActiveBandStructure::getEnergy(const long &stateIndex) {
  if (energies.size() == 0) {
    Error e("ActiveBandStructure energies haven't been populated");
  }
  return energies[stateIndex];
}

const double &ActiveBandStructure::getEnergy(StateIndex &is) {
  return getEnergy(is.get());
}

Eigen::VectorXd ActiveBandStructure::getEnergies(WavevectorIndex &ik) {
  long ikk = ik.get();
  long nb = numBands(ikk);
  Eigen::VectorXd x(nb);
  for (long ib = 0; ib < nb; ib++) {
    long ind = bloch2Comb(ikk, ib);
    x(ib) = energies[ind];
  }
  return x;
}

Eigen::Vector3d ActiveBandStructure::getGroupVelocity(const long &stateIndex) {
  if (velocities.size() == 0) {
    Error e("ActiveBandStructure velocities haven't been populated");
  }
  auto tup = comb2Bloch(stateIndex);
  auto ik = std::get<0>(tup);
  auto ib = std::get<1>(tup);
  Eigen::Vector3d vel;
  vel(0) = velocities[velBloch2Comb(ik, ib, ib, 0)].real();
  vel(1) = velocities[velBloch2Comb(ik, ib, ib, 1)].real();
  vel(2) = velocities[velBloch2Comb(ik, ib, ib, 2)].real();
  return vel;
}

Eigen::Vector3d ActiveBandStructure::getGroupVelocity(StateIndex &is) {
  return getGroupVelocity(is.get());
}

Eigen::MatrixXd ActiveBandStructure::getGroupVelocities(WavevectorIndex &ik) {
  long ikk = ik.get();
  long nb = numBands(ikk);
  Eigen::MatrixXd vel(nb, 3);
  for (long ib = 0; ib < nb; ib++) {
    for (long i : {0, 1, 2}) {
      vel(ib, i) = velocities[velBloch2Comb(ikk, ib, ib, i)].real();
    }
  }
  return vel;
}

Eigen::Tensor<std::complex<double>, 3> ActiveBandStructure::getVelocities(
    WavevectorIndex &ik) {
  long ikk = ik.get();
  long nb = numBands(ikk);
  Eigen::Tensor<std::complex<double>, 3> vel(nb, nb, 3);
  for (long ib1 = 0; ib1 < nb; ib1++) {
    for (long ib2 = 0; ib2 < nb; ib2++) {
      for (long i : {0, 1, 2}) {
        vel(ib1, ib2, i) = velocities[velBloch2Comb(ikk, ib1, ib2, i)];
      }
    }
  }
  return vel;
}

Eigen::MatrixXcd ActiveBandStructure::getEigenvectors(WavevectorIndex &ik) {
  long ikk = ik.get();
  long nb = numBands(ikk);
  Eigen::MatrixXcd eigs(numFullBands, nb);
  eigs.setZero();
  for (int ib1 = 0; ib1 < numFullBands; ib1++) {
    for (int ib2 = 0; ib2 < nb; ib2++) {
      long ind = eigBloch2Comb(ikk, ib1, ib2);
      eigs(ib1, ib2) = eigenvectors[ind];
    }
  }
  return eigs;
}

Eigen::Tensor<std::complex<double>, 3> ActiveBandStructure::getPhEigenvectors(
    WavevectorIndex &ik) {
  Eigen::MatrixXcd eigsMatrix = getEigenvectors(ik);
  int ikk = ik.get();
  int numAtoms = numFullBands / 3;
  Eigen::Tensor<std::complex<double>, 3> eigs(3, numAtoms, numBands(ikk));
  for (long i = 0; i < numFullBands; i++) {
    auto tup = decompress2Indeces(i, numAtoms, 3);
    auto iat = std::get<0>(tup);
    auto ic = std::get<1>(tup);
    for (long ib2 = 0; ib2 < numBands(ikk); ib2++) {
      eigs(ic, iat, ib2) = eigsMatrix(i, ib2);
    }
  }
  return eigs;
}

Eigen::Vector3d ActiveBandStructure::getWavevector(const long &stateIndex) {
  auto tup = getIndex(stateIndex);
  auto ik = std::get<0>(tup);
  return getWavevector(ik);
}

Eigen::Vector3d ActiveBandStructure::getWavevector(StateIndex &is) {
  return getWavevector(is.get());
}

Eigen::Vector3d ActiveBandStructure::getWavevector(WavevectorIndex &ik) {
  Point p = activePoints.getPoint(ik.get());
  return p.getCoords(Points::cartesianCoords, true);
}

double ActiveBandStructure::getWeight(const long &stateIndex) {
  auto tup = getIndex(stateIndex);
  auto ik = std::get<0>(tup);
  return getWeight(ik);
}

double ActiveBandStructure::getWeight(StateIndex &is) { return getWeight(is.get()); }

double ActiveBandStructure::getWeight(WavevectorIndex &ik) {
  return activePoints.getWeight(ik.get());
}

void ActiveBandStructure::setEnergies(Point &point,
                                      Eigen::VectorXd &energies_) {
  long ik = point.getIndex();
  for (long ib = 0; ib < energies_.size(); ib++) {
    long index = bloch2Comb(ik, ib);
    energies[index] = energies_(ib);
  }
}

void ActiveBandStructure::setEnergies(Point &point,
                                      std::vector<double> &energies_) {
  long ik = point.getIndex();
  for (long unsigned ib = 0; ib < energies_.size(); ib++) {
    long index = bloch2Comb(ik, ib);
    energies[index] = energies_[ib];
  }
}

void ActiveBandStructure::setEigenvectors(Point &point,
                                          Eigen::MatrixXcd &eigenvectors_) {
  long ik = point.getIndex();
  for (long i = 0; i < eigenvectors_.rows(); i++) {
    for (long j = 0; j < eigenvectors_.cols(); j++) {
      long index = eigBloch2Comb(ik, i, j);
      eigenvectors[index] = eigenvectors_(i, j);
    }
  }
}

void ActiveBandStructure::setVelocities(
    Point &point, Eigen::Tensor<std::complex<double>, 3> &velocities_) {
  long ik = point.getIndex();
  for (long ib1 = 0; ib1 < velocities_.dimension(0); ib1++) {
    for (long ib2 = 0; ib2 < velocities_.dimension(1); ib2++) {
      for (long j : {0, 1, 2}) {
        long index = velBloch2Comb(ik, ib1, ib2, j);
        velocities[index] = velocities_(ib1, ib2, j);
      }
    }
  }
}

long ActiveBandStructure::velBloch2Comb(const long &ik, const long &ib1,
                                        const long &ib2, const long &i) {
  return cumulativeKbbOffset(ik) + ib1 * numBands(ik) * 3 + ib2 * 3 + i;
}

long ActiveBandStructure::eigBloch2Comb(const long &ik, const long &ib1,
                                        const long &ib2) {
  return cumulativeKbOffset(ik) * numFullBands + ib1 * numBands(ik) + ib2;
}

long ActiveBandStructure::bloch2Comb(const long &ik, const long &ib) {
  return cumulativeKbOffset(ik) + ib;
}

std::tuple<long, long> ActiveBandStructure::comb2Bloch(const long &is) {
  return {auxBloch2Comb(is, 0), auxBloch2Comb(is, 1)};
}

void ActiveBandStructure::buildIndeces() {
  auxBloch2Comb = Eigen::MatrixXi::Zero(numStates, 2);
  cumulativeKbOffset = Eigen::VectorXi::Zero(numPoints);
  cumulativeKbbOffset = Eigen::VectorXi::Zero(numPoints);

  for (long ik = 1; ik < numPoints; ik++) {
    cumulativeKbOffset(ik) = cumulativeKbOffset(ik - 1) + numBands(ik - 1);
    cumulativeKbbOffset(ik) =
        cumulativeKbbOffset(ik - 1) + 3 * numBands(ik - 1) * numBands(ik - 1);
  }

  long is = 0;
  for (long ik = 0; ik < numPoints; ik++) {
    for (long ib = 0; ib < numBands(ik); ib++) {
      auxBloch2Comb(is, 0) = ik;
      auxBloch2Comb(is, 1) = ib;
      is += 1;
    }
  }
}

std::tuple<ActiveBandStructure, StatisticsSweep> ActiveBandStructure::builder(
    Context &context, HarmonicHamiltonian &h0, Points &points,
    const bool &withEigenvectors, const bool &withVelocities) {
  Particle particle = h0.getParticle();

  ActivePoints bogusPoints(points, Eigen::VectorXi::Zero(1));
  ActiveBandStructure activeBandStructure(particle, bogusPoints);

  if (particle.isPhonon()) {
    Eigen::VectorXd temperatures = context.getTemperatures();
    double temperatureMin = temperatures.minCoeff();
    double temperatureMax = temperatures.maxCoeff();

    Window window(context, particle, temperatureMin, temperatureMax);

    activeBandStructure.buildOnTheFly(window, points, h0, withEigenvectors,
                                      withVelocities);
    StatisticsSweep statisticsSweep(context);
    return {activeBandStructure, statisticsSweep};
  } else {
    StatisticsSweep s = activeBandStructure.buildAsPostprocessing(
        context, h0, points, withEigenvectors, withVelocities);
    return {activeBandStructure, s};
  }
}
// TODO check this out, may or may not need modificiation
void ActiveBandStructure::buildOnTheFly(Window &window, Points &points,
                                        HarmonicHamiltonian &h0,
                                        const bool &withEigenvectors,
                                        const bool &withVelocities) {
  // this function proceeds in three logical blocks:
  // 1- we find out the list of "relevant" points
  // 2- initialize internal raw buffer for energies, velocities, eigvecs
  // 3- populate the raw buffer

  // we have to build this in a way that works in parallel
  // ALGORITHM:
  // - loop over points. Diagonalize, and find if we want this k-point
  //   (while we are at it, we could save energies and the eigenvalues)
  // - find how many points each MPI rank has found
  // - communicate the indeces
  // - loop again over wavevectors to compute energies and velocities

  numFullBands = 0;  // save the unfiltered number of bands
  std::vector<int> myFilteredPoints;
  std::vector<std::vector<int>> myFilteredBands;

  for (long ik : mpi->divideWorkIter(points.getNumPoints())) {
    Point point = points.getPoint(ik);
    // diagonalize harmonic hamiltonian
    auto tup = h0.diagonalize(point);
    auto theseEnergies = std::get<0>(tup);
    auto theseEigenvectors = std::get<1>(tup);
    // ens is empty if no "relevant" energy is found.
    // bandsExtrema contains the lower and upper band index of "relevant"
    // bands at this point
    auto tup1 = window.apply(theseEnergies);
    auto ens = std::get<0>(tup1);
    auto bandsExtrema = std::get<1>(tup1);
    if (ens.empty()) {  // nothing to do
      continue;
    } else {  // save point index and "relevant" band indices
      myFilteredPoints.push_back(ik);
      myFilteredBands.push_back(bandsExtrema);
    }
    numFullBands = theseEnergies.size();
  }

  // now, we let each MPI process now how many points each process has found
  int myNumPts = myFilteredPoints.size();
  int mpiSize = mpi->getSize();
  int *receiveCounts = nullptr;
  receiveCounts = new int[mpiSize];
  for (int i = 0; i < mpiSize; i++) *(receiveCounts + i) = 0;
  if (mpiSize > 1) {
    mpi->gather(&myNumPts, receiveCounts);
#ifdef MPI_AVAIL
    // note: bcast interface doesn't work for objects of unknown size
    // (here, we'd need to pass the size to the interface)
    // convert receiveCounts to std::vector?
    MPI_Bcast(receiveCounts, mpiSize, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  } else {
    receiveCounts[0] = myNumPts;
  }

  // receivecounts now tells how many wavevectors found per MPI process

  // now we count the total number of wavevectors
  numPoints = 0;
  for (int i = 0; i < mpi->getSize(); i++) {
    numPoints += *(receiveCounts + i);
  }

  // now we collect the wavevector indeces
  // first we find the offset to compute global indeces from local indices
  std::vector<int> displacements(mpiSize, 0);
  for (int i = 1; i < mpiSize; i++) {
    displacements[i] = displacements[i - 1] + *(receiveCounts + i - 1);
  }

  // collect all the indeces in the filteredPoints vector
  Eigen::VectorXi filter(numPoints);
  filter.setZero();
  for (int i = 0; i < myNumPts; i++) {
    int index = i + displacements[mpi->getRank()];
    filter(index) = myFilteredPoints[i];
  }
  mpi->allReduceSum(&filter);

  // unfortunately, a vector<vector> isn't contiguous
  // let's use Eigen matrices
  Eigen::MatrixXi filteredBands(numPoints, 2);
  filteredBands.setZero();
  for (int i = 0; i < myNumPts; i++) {
    int index = i + displacements[mpi->getRank()];
    filteredBands(index, 0) = myFilteredBands[i][0];
    filteredBands(index, 1) = myFilteredBands[i][1];
  }
  mpi->allReduceSum(&filteredBands);

  //    free(receiveCounts);
  delete[] receiveCounts;

  //////////////// Done MPI recollection

  // numBands is a book-keeping of how many bands per kpoint there are
  // this isn't a constant number.
  // on top of that, we look for the size of the arrays containing bandstruc.
  numBands = Eigen::VectorXi::Zero(numPoints);
  long numEnStates = 0;
  long numVelStates = 0;
  long numEigStates = 0;
  for (long ik = 0; ik < numPoints; ik++) {
    numBands(ik) = filteredBands(ik, 1) - filteredBands(ik, 0) + 1;
    //
    numEnStates += numBands(ik);
    numVelStates += 3 * numBands(ik) * numBands(ik);
    numEigStates += numBands(ik) * numFullBands;
  }
  numStates = numEnStates;

  // initialize the raw data buffers of the activeBandStructure

  ActivePoints activePoints_(points, filter);
  activePoints = activePoints_;
  // construct the mapping from combined indices to Bloch indices
  buildIndeces();

  energies.resize(numEnStates, 0.);
  if (withVelocities) {
    velocities.resize(numVelStates, complexZero);
  }
  if (withEigenvectors) {
    hasEigenvectors = true;
    eigenvectors.resize(numEigStates, complexZero);
  }

  windowMethod = window.getMethodUsed();

/////////////////

// now we can loop over the trimmed list of points
#pragma omp parallel for
  for (long ik : mpi->divideWorkIter(numPoints)) {
    Point point = activePoints.getPoint(ik);
    auto tup = h0.diagonalize(point);
    auto theseEnergies = std::get<0>(tup);
    auto theseEigenvectors = std::get<1>(tup);
    // eigenvectors(3,numAtoms,numBands)
    auto tup1 = window.apply(theseEnergies);
    auto ens = std::get<0>(tup1);
    auto bandsExtrema = std::get<1>(tup1);

    Eigen::VectorXd eigEns(numBands(ik));
    long ibAct = 0;
    for (long ibFull = filteredBands(ik, 0); ibFull <= filteredBands(ik, 1);
         ibFull++) {
      eigEns(ibAct) = theseEnergies(ibFull);
      ibAct++;
    }
    setEnergies(point, eigEns);

    if (withEigenvectors) {
      // we are reducing the basis size!
      // the first index has the size of the Hamiltonian
      // the second index has the size of the filtered bands
      Eigen::MatrixXcd theseEigvecs(numFullBands, numBands(ik));
      long ibAct = 0;
      for (long ibFull = filteredBands(ik, 0); ibFull <= filteredBands(ik, 1);
           ibFull++) {
        theseEigvecs.col(ibAct) = theseEigenvectors.col(ibFull);
        ibAct++;
      }
      setEigenvectors(point, theseEigvecs);
    }

    if (withVelocities) {
      // thisVelocity is a tensor of dimensions (ib, ib, 3)
      auto thisVelocity = h0.diagonalizeVelocity(point);

      // now we filter it
      Eigen::Tensor<std::complex<double>, 3> thisVels(numBands(ik),
                                                      numBands(ik), 3);
      long ib1New = 0;
      for (long ib1Old = filteredBands(ik, 0);
           ib1Old < filteredBands(ik, 1) + 1; ib1Old++) {
        long ib2New = 0;
        for (long ib2Old = filteredBands(ik, 0);
             ib2Old < filteredBands(ik, 1) + 1; ib2Old++) {
          for (long i = 0; i < 3; i++) {
            thisVels(ib1New, ib2New, i) = thisVelocity(ib1Old, ib2Old, i);
          }
          ib2New++;
        }
        ib1New++;
      }
      setVelocities(point, thisVels);
    }
  }
  mpi->allReduceSum(&energies);
  mpi->allReduceSum(&velocities);
  mpi->allReduceSum(&eigenvectors);
}

/** in this function, useful for electrons, we first compute the bandstructure
 * on a dense grid of wavevectors, then compute chemical potential/temperatures
 * and then filter it
 */
StatisticsSweep ActiveBandStructure::buildAsPostprocessing(
    Context &context, HarmonicHamiltonian &h0, Points &points,
    const bool &withEigenvectors, const bool &withVelocities) {
  bool tmpWithVel_ = false;
  bool tmpWithEig_ = true;
  bool tmpIsDistributed_ = false; // TODO temporary, we need to pass this instead

  FullBandStructure fullBandStructure =
      h0.populate(points, tmpWithVel_, tmpWithEig_, tmpIsDistributed_);

  // This will work even if fullbandstructure is distributed
  StatisticsSweep statisticsSweep(context, &fullBandStructure);

  // find min/max value of temperatures and chemical potentials
  int numCalcs = statisticsSweep.getNumCalcs();
  std::vector<double> chemPots;
  std::vector<double> temps;
  for ( int i=0; i<numCalcs; i++ ) {
    auto calcStat = statisticsSweep.getCalcStatistics(i);
    chemPots.push_back(calcStat.chemicalPotential);
    temps.push_back(calcStat.temperature);
  }

  double temperatureMin = *min_element(temps.begin(), temps.end());
  double temperatureMax = *max_element(temps.begin(), temps.end());
  double chemicalPotentialMin = *min_element(chemPots.begin(), chemPots.end());
  double chemicalPotentialMax = *max_element(chemPots.begin(), chemPots.end());

  // now we can apply the window
  Window window(context, particle, temperatureMin, temperatureMax,
                chemicalPotentialMin, chemicalPotentialMax);

  numFullBands = h0.getNumBands();
  std::vector<int> myFilteredPoints;
  std::vector<std::vector<int>> myFilteredBands;

  // if all processes have the same points, divide up the points across 
  // processes. If this bandstructure was already distributed, then we 
  // can just perform this for the wavevectors belonging to each process's 
  // part of the distributed bandstructure.  
  std::vector<long> parallelIter; 
  if(!tmpIsDistributed_)
    parallelIter = mpi->divideWorkIter(points.getNumPoints());
  else 
    parallelIter = fullBandStructure.getWavevectorIndices();

  // iterate over mpi-parallelized wavevectors
  for (long ik : parallelIter) {

    auto ikIndex = WavevectorIndex(ik);
    Eigen::VectorXd theseEnergies = fullBandStructure.getEnergies(ikIndex);

    // ens is empty if no "relevant" energy is found.
    // bandsExtrema contains the lower and upper band index of "relevant"
    // bands at this point
    auto tup1 = window.apply(theseEnergies);
    auto ens = std::get<0>(tup1);
    auto bandsExtrema = std::get<1>(tup1);
    if (ens.empty()) {  // nothing to do
      continue;
    } else {  // save point index and "relevant" band indices
      myFilteredPoints.push_back(ik);
      myFilteredBands.push_back(bandsExtrema);
    }
    numFullBands = theseEnergies.size();
  }

  // now that we've counted up the selected points and their 
  // indices on each process, we need to reduce 
  int myNumPts = myFilteredPoints.size(); // TODO can we rename this, it's not super clear
  int mpiSize = mpi->getSize();
  int *receiveCounts = nullptr;
  receiveCounts = new int[mpiSize];
  for (int i = 0; i < mpiSize; i++) *(receiveCounts + i) = 0;
  if (mpiSize > 1) { // TODO this might be unnecessary once we add bcast function

    // take the number of kpoints of each process and fill 
    // buffer receiveCounts with these values
    mpi->gather(&myNumPts, receiveCounts);

    #ifdef MPI_AVAIL
    // TODO add a function to mpiController for bcast with a function call
    // taking receive counts
    // note: bcast interface doesn't work for objects of unknown size
    // (here, we'd need to pass the size to the interface)
    // convert receiveCounts to std::vector?

    // are we here now broadcasting the receive array to everyone? 
    // wouldn't it make more sense to do an allGather in the first place? 
    MPI_Bcast(receiveCounts, mpiSize, MPI_INT, 0, MPI_COMM_WORLD);
    #endif
  } else {
    receiveCounts[0] = myNumPts;
  }
  // receivecounts now tells how many wavevectors found per MPI process

  // now we count the total number of wavevectors
  // by summing over receive counts
  numPoints = 0;
  for (int i = 0; i < mpi->getSize(); i++) {
    numPoints += *(receiveCounts + i);
  }

  // now we collect the wavevector indices
  // first we find the offset to compute global indices from local indices
  std::vector<int> displacements(mpiSize, 0);
  for (int i = 1; i < mpiSize; i++) {
    displacements[i] = displacements[i - 1] + *(receiveCounts + i - 1);
  }

  // collect all the indices in the filteredPoints vector
  Eigen::VectorXi filter(numPoints);
  filter.setZero();
  for (int i = 0; i < myNumPts; i++) {
    int index = i + displacements[mpi->getRank()];
    filter(index) = myFilteredPoints[i];
  }
  mpi->allReduceSum(&filter);

  // unfortunately, a vector<vector> isn't contiguous
  // let's use Eigen matrices
  Eigen::MatrixXi filteredBands(numPoints, 2);
  filteredBands.setZero();
  for (int i = 0; i < myNumPts; i++) {
    int index = i + displacements[mpi->getRank()];
    filteredBands(index, 0) = myFilteredBands[i][0];
    filteredBands(index, 1) = myFilteredBands[i][1];
  }
  mpi->allReduceSum(&filteredBands);

  delete[] receiveCounts;

  //////////////// Done MPI recollection

  // numBands is a book-keeping of how many bands per kpoint there are
  // this isn't a constant number.
  // on top of that, we look for the size of the arrays containing bandstruc.
  numBands = Eigen::VectorXi::Zero(numPoints);
  long numEnStates = 0;
  long numVelStates = 0;
  long numEigStates = 0;
  for (long ik = 0; ik < numPoints; ik++) {
    numBands(ik) = filteredBands(ik, 1) - filteredBands(ik, 0) + 1;
    numEnStates += numBands(ik);
    numVelStates += 3 * numBands(ik) * numBands(ik);
    numEigStates += numBands(ik) * numFullBands;
  }
  numStates = numEnStates;

  // initialize the raw data buffers of the activeBandStructure

  ActivePoints activePoints_(points, filter);
  activePoints = activePoints_;

  // construct the mapping from combined indices to Bloch indices
  buildIndeces();

  energies.resize(numEnStates, 0.);
  if (withVelocities) {
    velocities.resize(numVelStates, complexZero);
  }
  if (withEigenvectors) {
    hasEigenvectors = true;
    eigenvectors.resize(numEigStates, complexZero);
  }
  windowMethod = window.getMethodUsed();

  /////////////////

  // Now we can loop over the trimmed list of points.
  // To accomodate the case where FullBS is distributed, 
  // we save the energies related to myFilteredPoints/Bands 
  // and then allReduce or allGather those instead

  for (int i=0; i<myFilteredPoints.size(); i++) {

    long ik = myFilteredPoints[i];
    auto ikIndex = WavevectorIndex(ik);

    // use displacement array to get global idx of this point
    // within the list of activePoints
    int ikg = i + displacements[mpi->getRank()];
    Point point = activePoints.getPoint(ikg);

    Eigen::VectorXd theseEnergies = fullBandStructure.getEnergies(ikIndex);
    Eigen::MatrixXcd theseEigenvectors = fullBandStructure.getEigenvectors(ikIndex);

    // copy energies into internal storage
    Eigen::VectorXd eigEns(numBands(i));
    long ibAct = 0;
    for (long ibFull = filteredBands(ikg, 0); ibFull <= filteredBands(ikg, 1); ibFull++) {
      eigEns(ibAct) = theseEnergies(ibFull);
      ibAct++;
    }
    setEnergies(point, eigEns);

    // copy eigenvectors into internal storage
    if (withEigenvectors) {
      // we are reducing the basis size!
      // the first index has the size of the Hamiltonian
      // the second index has the size of the filtered bands
      Eigen::MatrixXcd theseEigvecs(numFullBands, numBands(i));
      long ibAct = 0;
      for (long ibFull = filteredBands(ikg, 0); ibFull <= filteredBands(ikg, 1); ibFull++) {
        theseEigvecs.col(ibAct) = theseEigenvectors.col(ibFull);
        ibAct++;
      }
      setEigenvectors(point, theseEigvecs);
    }
  }
  // all reduce sum the energies and eigenvectors here  
  mpi->allReduceSum(&energies);
  if (withEigenvectors) mpi->allReduceSum(&eigenvectors);
  
  // compute velocities
  if (withVelocities) {

    // loop over the points available to this process
    #pragma omp parallel for
    for (int i=0; i<myFilteredPoints.size(); i++) {

      long ik = myFilteredPoints[i];
      auto ikIndex = WavevectorIndex(ik);
      int ikg = i + displacements[mpi->getRank()];
      Point point = activePoints.getPoint(ikg);

      // thisVelocity is a tensor of dimensions (ib, ib, 3)
      auto thisVelocity = h0.diagonalizeVelocity(point);

      // now we filter it
      Eigen::Tensor<std::complex<double>, 3> thisVels(numBands(ikg),
                                                      numBands(ikg), 3);
      long ib1New = 0;
      for (long ib1Old = filteredBands(ikg, 0); ib1Old < filteredBands(ikg, 1) + 1; ib1Old++) {
        long ib2New = 0;
        for (long ib2Old = filteredBands(ikg, 0); ib2Old < filteredBands(ikg, 1) + 1; ib2Old++) {
          for (long ic = 0; ic < 3; ic++) {
            thisVels(ib1New, ib2New, ic) = thisVelocity(ib1Old, ib2Old, ic);
          }
          ib2New++;
        }
        ib1New++;
      }
      setVelocities(point, thisVels);
    }
    mpi->allReduceSum(&velocities);
  }
  return statisticsSweep;
}
