#include "active_bandstructure.h"
#include "bandstructure.h"
#include "exceptions.h"
#include "window.h"
#include "mpiHelper.h"
#include <cstdlib>

ActiveBandStructure::ActiveBandStructure(Particle &particle_,
        ActivePoints &activePoints) :
        particle(particle_), activePoints(activePoints) {
}

// copy constructor
ActiveBandStructure::ActiveBandStructure(const ActiveBandStructure &that) :
        particle(that.particle), activePoints(that.activePoints), energies(
                that.energies), velocities(that.velocities), eigenvectors(
                that.eigenvectors), hasEigenvectors(that.hasEigenvectors),
                numStates(that.numStates), numPoints(that.numPoints), numBands(
                that.numBands), numFullBands(that.numFullBands), windowMethod(
                that.windowMethod), auxBloch2Comb(that.auxBloch2Comb),
                cumulativeKbOffset(that.cumulativeKbOffset),
                cumulativeKbbOffset(that.cumulativeKbbOffset) {
}

ActiveBandStructure& ActiveBandStructure::operator=(
        const ActiveBandStructure &that) { // assignment operator
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

Particle ActiveBandStructure::getParticle() {
    return particle;
}

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
    if ( useFullGrid ) {
        return activePoints.getParentPoints().getNumPoints();
    } else { // default
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

long ActiveBandStructure::hasWindow() {
    return windowMethod;
}

State ActiveBandStructure::getState(const long &pointIndex) {
    Point point = getPoint(pointIndex);
    return getState(point);
}

State ActiveBandStructure::getState(Point &point) {
    if (!hasPoints()) {
        Error e("ActiveBandStructure hasn't been populated yet");
    }

    long ik = point.getIndex();

    long ind = bloch2Comb(ik, 0);
    double *thisEn = &energies[ind];

    std::complex<double> *thisVel = nullptr;
    if (velocities.size() != 0) {
        ind = velBloch2Comb(ik, 0, 0, 0);
        thisVel = &velocities[ind];
    }

    std::complex<double> *thisEig = nullptr;
    if (eigenvectors.size() != 0) {
        ind = eigBloch2Comb(ik, 0, 0);
        thisEig = &eigenvectors[ind];
    }

    State s(point, thisEn, numFullBands, numBands(ik), thisVel, thisEig);
    return s;
}

long ActiveBandStructure::getIndex(const WavevectorIndex &ik,
        const BandIndex &ib) {
    if (!hasPoints()) {
        Error e("ActiveBandStructure hasn't been populated yet");
    }
    return bloch2Comb(ik.get(), ib.get());
}

std::tuple<WavevectorIndex,BandIndex> ActiveBandStructure::getIndex(
        const long &is) {
    if (!hasPoints()) {
        Error e("ActiveBandStructure hasn't been populated yet");
    }
    auto [ik,ib] = comb2Bloch(is);
    auto ikk = WavevectorIndex(ik);
    auto ibb = BandIndex(ib);
    return {ikk,ibb};
}

long ActiveBandStructure::getNumStates() {
    if (!hasPoints()) {
        Error e("ActiveBandStructure hasn't been populated yet");
    }
    return numStates;
}

const double& ActiveBandStructure::getEnergy(const long &stateIndex) {
    if (energies.size() == 0) {
        Error e("ActiveBandStructure energies haven't been populated");
    }
    return energies[stateIndex];
}

Eigen::Vector3d ActiveBandStructure::getGroupVelocity(const long &stateIndex) {
    if (velocities.size() == 0) {
        Error e("ActiveBandStructure velocities haven't been populated");
    }
    auto[ik,ib] = comb2Bloch(stateIndex);
    Eigen::Vector3d vel;
    vel(0) = velocities[velBloch2Comb(ik, ib, ib, 0)].real();
    vel(1) = velocities[velBloch2Comb(ik, ib, ib, 1)].real();
    vel(2) = velocities[velBloch2Comb(ik, ib, ib, 2)].real();
    return vel;
}

Eigen::Vector3d ActiveBandStructure::getWavevector(const long &stateIndex) {
    if (!hasPoints()) {
        Error e("ActiveBandStructure hasn't been populated yet");
    }
    auto[ik,ib] = comb2Bloch(stateIndex);
    Point p = activePoints.getPoint(ik);
    return p.getCoords(Points::cartesianCoords);
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

void ActiveBandStructure::setVelocities(Point &point,
        Eigen::Tensor<std::complex<double>, 3> &velocities_) {
    long ik = point.getIndex();
    for (long ib1 = 0; ib1 < velocities_.dimension(0); ib1++) {
        for (long ib2 = 0; ib2 < velocities_.dimension(1); ib2++) {
            for (long j : { 0, 1, 2 }) {
                long index = velBloch2Comb(ik, ib1, ib2, j);
                velocities[index] = velocities_(ib1, ib2, j);
            }
        }
    }
}

//ActivePoints ActiveBandStructure::buildAsPostprocessing(Window & window,
//		FullBandStructure & fullBandStructure) {
//
//	numFullBands = fullBandStructure.getNumBands();
//	auto fullPoints = fullBandStructure.getPoints();
//
//	std::vector<long> filteredPoints;
//	std::vector<std::vector<long>> filteredBands;
//	// first we find the states that can be kept
//	for ( long ik=0; ik<fullBandStructure.getNumPoints(); ik++ ) {
//		Point point = fullBandStructure.getPoint(ik);
//		State s = fullBandStructure.getState(point);
//		auto theseEnergies = s.getEnergies();
//		auto [ens, bandsExtrema] = window.apply(theseEnergies);
//		if ( ens.empty() ) {
//			continue;
//		} else {
//			filteredPoints.push_back(ik);
//			filteredBands.push_back(bandsExtrema);
//		}
//	}
//
//	// Here, we should use MPI_gather(v) to put filteredPoints in common
//
//	numPoints = filteredPoints.size();
//	// this vector mapps the indices of the new point to the old list
//	VectorXl filter(numPoints);
//	for ( long i=0; i<numPoints; i++ ) {
//		filter(i) = filteredPoints[i];
//	}
//
//	// numBands is a book-keeping of how many bands per kpoint there are
//	// this isn't a constant number.
//	// on top of that, we look for the size of the arrays containing bandstruc.
//	VectorXl numBands(numPoints);
//	long numEnStates = 0;
//	long numVelStates = 0;
//	long numEigStates = 0;
//	for ( long ik=0; ik<numPoints; ik++ ) {
//		numBands(ik) = filteredBands[ik][1] - filteredBands[ik][0] + 1;
//		//
//		numEnStates += numBands(ik);
//		numVelStates += 3 * numBands(ik) * numBands(ik);
//		numEigStates += numBands(ik) * numFullBands;
//	}
//
//	// initialize the raw data buffers of the activeBandStructure
//
//	ActivePoints activePoints_(fullPoints, filter);
//	activePoints = activePoints_;
//	// construct the mapping from combined indices to Bloch indices
//	buildIndeces();
//
//	energies.resize(numEnStates,0.);
//	velocities.resize(numVelStates,complexZero);
//	if ( fullBandStructure.hasEigenvectors ) {
//		hasEigenvectors = true;
//		eigenvectors.resize(numEigStates,complexZero);
//	}
//
//	// now we store the energies and other quantities in the reordered array
//	// loop over new points
//
//	for ( long ik=0; ik<numPoints; ik++ ) {
//
//		auto point = activePoints.getPoint(ik);
//		long ikOld = filter(ik);
//
//		// load the old state
//		auto state = fullBandStructure.getState(ikOld);
//
//		// copy the filtered energies
//		Eigen::VectorXd thisEnergies = state.getEnergies();
//		Eigen::VectorXd eigEns(numBands(ik));
//		long ibAct = 0;
//		for ( long ibFull=filteredBands[ik][0];
//				ibFull<=filteredBands[ik][1]; ibFull++) {
//			eigEns(ibAct) = thisEnergies(ibFull);
//			ibAct++;
//		}
//		setEnergies(point, eigEns);
//
//		if ( fullBandStructure.hasEigenvectors ) {
//			Eigen::MatrixXcd theseEigenvectors;
//			state.getEigenvectors(theseEigenvectors);
//
//			// we are reducing the basis size!
//			// the first index has the size of the Hamiltonian
//			// the second index has the size of the filtered bands
//			Eigen::MatrixXcd theseEigvecs(numFullBands,numBands(ik));
//			long ibAct = 0;
//			for ( long ibFull=filteredBands[ik][0];
//					ibFull<=filteredBands[ik][1]; ibFull++) {
//				theseEigvecs.col(ibAct) = theseEigenvectors.col(ibFull);
//				ibAct++;
//			}
//			setEigenvectors(point, theseEigvecs);
//		}
//
//		// thisVelocity is a tensor of dimensions (ib, ib, 3)
//		auto thisVelocity = state.getVelocities();
//
//		// now we filter it
//		Eigen::Tensor<std::complex<double>,3> thisVels(numBands(ik),
//				numBands(ik),3);
//		long ib1New = 0;
//		for ( long ib1Old = filteredBands[ik][0];
//				ib1Old<filteredBands[ik][1]+1; ib1Old++ ) {
//			long ib2New = 0;
//			for ( long ib2Old = filteredBands[ik][0];
//					ib2Old<filteredBands[ik][1]+1; ib2Old++ ) {
//				for ( long i=0; i<3; i++ ) {
//					thisVels(ib1New,ib2New,i) = thisVelocity(ib1Old,ib2Old,i);
//				}
//				ib2New++;
//			}
//			ib1New++;
//		}
//		setVelocities(point, thisVels);
//	}
//
//	return activePoints_;
//}

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
    return {auxBloch2Comb(is,0),auxBloch2Comb(is,1)};
}

void ActiveBandStructure::buildIndeces() {
    auxBloch2Comb = Eigen::MatrixXi::Zero(numStates, 2);
    cumulativeKbOffset = Eigen::VectorXi::Zero(numPoints);
    cumulativeKbbOffset = Eigen::VectorXi::Zero(numPoints);

    for (long ik = 1; ik < numPoints; ik++) {
        cumulativeKbOffset(ik) = cumulativeKbOffset(ik - 1) + numBands(ik - 1);
        cumulativeKbbOffset(ik) = cumulativeKbbOffset(ik - 1)
                + 3 * numBands(ik - 1) * numBands(ik - 1);
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
        Context &context, HarmonicHamiltonian &h0, FullPoints &fullPoints,
        const bool &withEigenvectors, const bool &withVelocities) {

    Particle particle = h0.getParticle();

    Eigen::VectorXd temperatures = context.getTemperatures();

    ActivePoints bogusPoints(fullPoints, Eigen::VectorXi::Zero(1));
    ActiveBandStructure activeBandStructure(particle, bogusPoints);

    if (particle.isPhonon()) {
        double temperatureMin = temperatures.minCoeff();
        double temperatureMax = temperatures.maxCoeff();

        Window window(context, particle, temperatureMin, temperatureMax);

        activeBandStructure.buildOnTheFly(window, fullPoints, h0,
                withEigenvectors, withVelocities);
        StatisticsSweep statisticsSweep(context);
        return {activeBandStructure, statisticsSweep};
    } else {
        Error e("apply window for electrons not implemented");
    }
}

void ActiveBandStructure::buildOnTheFly(Window &window, FullPoints &fullPoints,
        HarmonicHamiltonian &h0, const bool &withEigenvectors,
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

    numFullBands = 0; // save the unfiltered number of bands
    std::vector<int> myFilteredPoints;
    std::vector<std::vector<long>> myFilteredBands;
    for (long ik : mpi->divideWorkIter(fullPoints.getNumPoints())) {
      //    for (long ik = 0; ik < fullPoints.getNumPoints(); ik++) {
      Point point = fullPoints.getPoint(ik);
      // diagonalize harmonic hamiltonian
      auto [theseEnergies, theseEigenvectors] = h0.diagonalize(point);
      // ens is empty if no "relevant" energy is found.
      // bandsExtrema contains the lower and upper band index of "relevant"
      // bands at this point
      auto [ens, bandsExtrema] = window.apply(theseEnergies);
      if (ens.empty()) {  // nothing to do
        continue;
      } else {  // save point index and "relevant" band indices
        myFilteredPoints.push_back(ik);
        myFilteredBands.push_back(bandsExtrema);
      }
      numFullBands = theseEnergies.size();
    }

    // now, we let each MPI process now how many points each process has found
    int mySize = myFilteredPoints.size();
    int mpiSize = mpi->getSize();
    int *receiveCounts = nullptr;
    receiveCounts = (int *)malloc(mpiSize * sizeof(int));
    mpi->gather(&mySize, receiveCounts);
#ifdef MPI_AVAIL
    MPI_Bcast(receiveCounts, mpiSize, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    // receivecounts now tells how many wavevectors found per MPI process

    // now we count the total number of wavevectors
    numPoints = 0;
    for ( int i=0; i<mpi->getSize(); i++ ) {
      numPoints += *(receiveCounts+i);
    }

    // now we collect the wavevector indeces
    // first we find the offset to compute global indeces from local indices
    std::vector<int> displacements(mpiSize,0);
    for (int i = 1; i < mpiSize; i++) {
      displacements[i] = displacements[i - 1] + *(receiveCounts + i - 1);
    }

    // collect all the indeces in the filteredPoints vector
    Eigen::VectorXi filter(numPoints);
    filter.setZero();
    for ( int i=0; i<mySize; i++ ) {
      int index = i + displacements[mpi->getRank()];
      filter(index) = myFilteredPoints[i];
    }
    mpi->allReduceSum(&filter);

    // unfortunately, a vector<vector> isn't contiguous
    // let's use Eigen matrices
    Eigen::MatrixXi filteredBands(numPoints,2);
    filteredBands.setZero();
    for ( int i=0; i<mySize; i++ ) {
      int index = i + displacements[mpi->getRank()];
      filteredBands(index,0) = myFilteredBands[i][0];
      filteredBands(index,1) = myFilteredBands[i][1];
    }
    mpi->allReduceSum(&filteredBands);

    free(receiveCounts);

    //////////////// Done MPI recollection

    // numBands is a book-keeping of how many bands per kpoint there are
    // this isn't a constant number.
    // on top of that, we look for the size of the arrays containing bandstruc.
    numBands = Eigen::VectorXi::Zero(numPoints);
    long numEnStates = 0;
    long numVelStates = 0;
    long numEigStates = 0;
    for (long ik = 0; ik < numPoints; ik++) {
        numBands(ik) = filteredBands(ik,1) - filteredBands(ik,0) + 1;
        //
        numEnStates += numBands(ik);
        numVelStates += 3 * numBands(ik) * numBands(ik);
        numEigStates += numBands(ik) * numFullBands;
    }
    numStates = numEnStates;

    // initialize the raw data buffers of the activeBandStructure

    ActivePoints activePoints_(fullPoints, filter);
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
    for (long ik : mpi->divideWorkIter(numPoints)) {
        Point point = activePoints.getPoint(ik);
        auto [theseEnergies, theseEigenvectors] = h0.diagonalize(point);
        // eigenvectors(3,numAtoms,numBands)
        auto [ens, bandsExtrema] = window.apply(theseEnergies);

        Eigen::VectorXd eigEns(numBands(ik));
        long ibAct = 0;
        for (long ibFull = filteredBands(ik,0); ibFull <= filteredBands(ik,1);
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
            for (long ibFull = filteredBands(ik,0);
                    ibFull <= filteredBands(ik,1); ibFull++) {
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
            for (long ib1Old = filteredBands(ik,0);
                    ib1Old < filteredBands(ik,1) + 1; ib1Old++) {
                long ib2New = 0;
                for (long ib2Old = filteredBands(ik,0);
                        ib2Old < filteredBands(ik,1) + 1; ib2Old++) {
                    for (long i = 0; i < 3; i++) {
                        thisVels(ib1New, ib2New, i) = thisVelocity(ib1Old,
                                ib2Old, i);
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
