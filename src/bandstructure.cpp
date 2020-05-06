#include "bandstructure.h"
#include "exceptions.h"
#include "window.h"

long compressIndeces(long i, long j, long k, long size1, long size2,
		long size3) {
	return i*size2*size3 + j*size3 + k;
}

//auxiliary functions to work with the velocity
Eigen::VectorXcd packXcd(Eigen::Tensor<std::complex<double>,3> a,
		long size1, long size2, long size3) {
	Eigen::VectorXcd b(size1*size2*size3);
	for ( long i=0; i<size1; i++ ) {
		for ( long j=0; j<size2; j++ ) {
			for ( long k=0; k<size3; k++ ) {
				b(i*size2*size3 + j*size3 + k) = a(i,j,k);
			}
		}
	}
	return b;
}

std::vector<std::complex<double>> packStdXcd(Eigen::Tensor<std::complex<double>,3> a,
		long size1, long size2, long size3) {
	std::vector<std::complex<double>> b(size1*size2*size3);
	for ( long i=0; i<size1; i++ ) {
		for ( long j=0; j<size2; j++ ) {
			for ( long k=0; k<size3; k++ ) {
				b[i*size2*size3 + j*size3 + k] = a(i,j,k);
			}
		}
	}
	return b;
}

FullBandStructure::FullBandStructure(long numBands_, Statistics & statistics_,
		bool withVelocities, bool withEigenvectors, FullPoints & fullPoints_) :
			statistics{statistics_}, fullPoints(fullPoints_) {

	numBands = numBands_;
	numAtoms = numBands_ / 3;

//	if ( fullPoints_ != nullptr ) {
		useIrreducible = false;
//		fullPoints = fullPoints_;
//	} else if ( irreduciblePoints_ != nullptr ) {
//		useIrreducible = true;
//		irreduciblePoints = irreduciblePoints_;
//	} else {
//		Error e("FullBandStructure must provide one Points mesh.", 1);
//	}
//	if ( ( fullPoints != nullptr ) && ( irreduciblePoints != nullptr ) ) {
//		Error e("FullBandStructure must provide only one Points mesh.", 1);
//	}

	if ( withVelocities ) {
		hasVelocities = true;
		Eigen::MatrixXcd velocities_(getNumPoints(), numBands*numBands*3);
		velocities_.setZero();
		velocities = velocities_;
	}

	if ( withEigenvectors ) {
		hasVelocities = true;
		Eigen::MatrixXcd eigenvectors_(getNumPoints(),3*numAtoms*numBands);
		eigenvectors_.setZero();
		eigenvectors = eigenvectors_;
	}

	Eigen::MatrixXd energies_(getNumPoints(), numBands);
	energies.setZero();
	energies = energies_;

	// now, I want to manipulate the Eigen matrices at lower level
	// I create this pointer to data, so I can move it around
	rawEnergies = energies.data();
	if ( hasVelocities ) {
		rawVelocities = velocities.data();
	}
	if ( hasEigenvectors ) {
		rawEigenvectors = eigenvectors.data();
	}

	energiesRows = numBands;
	velocitiesRows = numBands * numBands * 3;
	eigenvectorsRows = numBands * numAtoms * 3;
}

// copy constructor
FullBandStructure::FullBandStructure(const FullBandStructure & that) :
	statistics(that.statistics), fullPoints(that.fullPoints),
	energies(that.energies), velocities(that.velocities),
	eigenvectors(that.eigenvectors), rawEnergies(that.rawEnergies),
	rawVelocities(that.rawVelocities), rawEigenvectors(that.rawEigenvectors),
	energiesRows(that.energiesRows), velocitiesRows(that.velocitiesRows),
	eigenvectorsRows(that.eigenvectorsRows), numBands(that.numBands),
	numAtoms(that.numAtoms), useIrreducible(that.useIrreducible),
	hasEigenvectors(that.hasEigenvectors), hasVelocities(that.hasVelocities) {
}

FullBandStructure & FullBandStructure::operator = ( // copy assignment
		const FullBandStructure & that) {
	if ( this != &that ) {
		statistics = that.statistics;
		fullPoints = that.fullPoints;
		energies = that.energies;
		velocities = that.velocities;
		eigenvectors = that.eigenvectors;
		rawEnergies = that.rawEnergies;
		rawVelocities = that.rawVelocities;
		rawEigenvectors = that.rawEigenvectors;
		energiesRows = that.energiesRows;
		velocitiesRows = that.velocitiesRows;
		eigenvectorsRows = that.eigenvectorsRows;
		numBands = that.numBands;
		numAtoms = that.numAtoms;
		useIrreducible = that.useIrreducible;
		hasEigenvectors = that.hasEigenvectors;
		hasVelocities = that.hasVelocities;
//		irreduciblePoints = that.irreduciblePoints;
	}
	return *this;
}

Statistics FullBandStructure::getStatistics() {
	return statistics;
}

Statistics ActiveBandStructure::getStatistics() {
	return statistics;
}

long FullBandStructure::getNumBands() {
	return numBands;
}

bool FullBandStructure::hasIrreduciblePoints() {
	return useIrreducible;
}

long FullBandStructure::getNumPoints() {
//	if ( useIrreducible ) {
//		return irreduciblePoints->getNumPoints();
//	} else {
	return fullPoints.getNumPoints();
//	}
}

long FullBandStructure::getIndex(Eigen::Vector3d& pointCoords) {
//	if ( useIrreducible ) {
//		return irreduciblePoints->getIndex(pointCoords);
//	} else {
		return fullPoints.getIndex(pointCoords);
//	}
}

Point FullBandStructure::getPoint(const long& pointIndex) {
//	if ( useIrreducible ) {
//		return irreduciblePoints->getPoint(pointIndex);
//	} else {
		return fullPoints.getPoint(pointIndex);
//	}
}

void FullBandStructure::setEnergies(Eigen::Vector3d& coords,
		Eigen::VectorXd& energies_) {
	long ik = getIndex(coords);
	energies.row(ik) = energies_;
}

void FullBandStructure::setEnergies(Point & point,
		Eigen::VectorXd& energies_) {
	long ik = point.getIndex();
	energies.row(ik) = energies_;
}

void FullBandStructure::setVelocities(Point & point,
		Eigen::Tensor<std::complex<double>,3>& velocities_) {
	if ( ! hasVelocities ) {
		Error e("FullBandStructure was initialized without velocities",1);
	}
	Eigen::VectorXcd tmpVelocities_(numBands*numBands*3);
	tmpVelocities_ = packXcd(velocities_, numBands, numBands, 3);
	long ik = point.getIndex();
	velocities.row(ik) = tmpVelocities_;
}

void FullBandStructure::setEigenvectors(Point & point,
		Eigen::Tensor<std::complex<double>,3>& eigenvectors_) {
	if ( ! hasEigenvectors ) {
		Error e("FullBandStructure was initialized without eigvecs",1);
	}
	Eigen::VectorXcd tmp = packXcd(eigenvectors_, 3, numAtoms, numBands);
	long ik = point.getIndex();
	eigenvectors.row(ik) = tmp;
}

State FullBandStructure::getState(Point & point) {
	long pointIndex = point.getIndex();

	// we construct the vector by defining begin() and end()
	double * thisEn;
	thisEn = rawEnergies + pointIndex * energiesRows;

	// note: in some cases these are nullptr
	std::complex<double> * thisVel = nullptr;
	std::complex<double> * thisEig = nullptr;

	if ( hasVelocities ) {
		thisVel = rawVelocities + pointIndex * velocitiesRows;
	}
	if ( hasEigenvectors ) {
		thisEig = rawEigenvectors + pointIndex * eigenvectorsRows;
	}

	State s(point, thisEn, numAtoms, numBands, thisVel, thisEig);
	return s;
}

void FullBandStructure::populate(PhononH0 & h0) {
	for ( long ik=0; ik<getNumPoints(); ik++ ) {
		Point point = getPoint(ik);
		auto [ens, eigvecs] = h0.diagonalize(point);
		setEnergies(point, ens);
		if ( hasEigenvectors ) {
			setEigenvectors(point, eigvecs);
		}
		if ( hasVelocities) {
			auto vels = h0.diagonalizeVelocity(point);
			setVelocities(point, vels);
		}
	}
}

Eigen::VectorXd FullBandStructure::getBandEnergies(long & bandIndex) {
	Eigen::VectorXd bandEnergies = energies.col(bandIndex);
	return bandEnergies;
}

bool ActiveBandStructure::hasPoints() {
	if ( activePoints != nullptr ) {
		return true;
	} else {
		return false;
	}
}

ActiveBandStructure::ActiveBandStructure(Statistics & statistics_) :
		statistics(statistics_) {
}

long ActiveBandStructure::getNumPoints() {
	if ( ! hasPoints() ) {
		Error e("ActiveBandStructure hasn't been populated yet" ,1);
	}
	return activePoints->getNumPoints();
}

Point ActiveBandStructure::getPoint(const long & pointIndex) {
	if ( ! hasPoints() ) {
		Error e("ActiveBandStructure hasn't been populated yet" ,1);
	}
	return activePoints->getPoint(pointIndex);
}

long ActiveBandStructure::getNumStates() {
	if ( ! hasPoints() ) {
		Error e("ActiveBandStructure hasn't been populated yet" ,1);
	}
	return numStates;
}

State ActiveBandStructure::getState(Point & point) {
	if ( ! hasPoints() ) {
		Error e("ActiveBandStructure hasn't been populated yet" ,1);
	}

	long ik = point.getIndex();
	long zero = 0;

	long ind = bloch2Comb(ik,zero);
	double * thisEn = &energies[ind];

	ind = velBloch2Comb(ik,zero,zero,zero);
	std::complex<double> * thisVel = &velocities[ind];

	std::complex<double> * thisEig = nullptr;
	if ( hasEigenvectors ) {
		ind = eigBloch2Comb(ik,zero,zero,zero);
		thisEig = &eigenvectors[ind];
	}

	State s(point, thisEn, numAtoms, numBands(ik), thisVel, thisEig);
	return s;
}

//State getBlochState(const long & stateIndex); //return a single Bloch state

ActivePoints ActiveBandStructure::buildOnTheFly(Window & window,
		FullPoints & fullPoints, HarmonicHamiltonian & h0) {

	// Note: eigenvectors are assumed to be phonon eigenvectors
	// might need adjustments in the future.

	if ( ! window.canOnTheFly ) {
		Error e("Window cannot be applied in this way" ,1);
	}

	numAtoms = fullPoints.getCrystal().getNumAtoms();

	std::vector<long> filteredPoints;
	std::vector<std::vector<long>> filteredBands;
	std::vector<double> filteredEnergies;
	std::vector<std::complex<double>> filteredEigenvectors;

	for ( long ik=0; ik<fullPoints.getNumPoints(); ik++ ) {
		Point point = fullPoints.getPoint(ik);

		auto [theseEnergies, theseEigenvectors] = h0.diagonalize(point);
		// eigenvectors(3,numAtoms,numBands)

		auto [ens, bandsExtrema] = window.apply(theseEnergies);

		if ( ens.empty() ) {
			continue;
		} else {

			filteredPoints.push_back(ik);
			filteredBands.push_back(bandsExtrema);

			for ( long ib=0; ib<ens.size(); ib++ ) {
				filteredEnergies.push_back(ens[ib]);
			}

			if ( h0.hasEigenvectors ) {
				std::complex<double> x;

				for ( int i=0; i<3; i++ ) {
					for ( int iat=0; iat<numAtoms; iat++ ) {
						for ( long ibOld=bandsExtrema[0];
								ibOld<bandsExtrema[1]+1; ibOld++ ) {
							x = theseEigenvectors(i,iat,ibOld);
							filteredEigenvectors.push_back(x);
						}
					}
				}
			}
		}

	}

	numPoints = filteredPoints.size();
	VectorXl filter(numPoints);
	for ( long i=0; i<numPoints; i++ ) {
		filter(i) = filteredPoints[i];
	}

	VectorXl numBands(numPoints);
	for ( long ik=0; ik<numPoints; ik++ ) {
		numBands(ik) = filteredBands[ik][1] - filteredBands[ik][0] + 1;
	}

	numStates = numBands.sum();

	ActivePoints activePoints_(fullPoints, filter);
	activePoints = &activePoints_;

	energies = filteredEnergies;

	if ( h0.hasEigenvectors ) {
		hasEigenvectors = true;
		std::vector<std::complex<double>> eigenvectors(numStates*3*numAtoms);
		eigenvectors = filteredEigenvectors;
	}

	// now we add velocities
	for ( long ik=0; ik<numPoints; ik++ ) {
		Point point = activePoints_.getPoint(ik);

		// thisVelocity is a tensor of dimensions (ib, ib, 3)
		auto thisVelocity = h0.diagonalizeVelocity(point);

		// now we filter it
		for ( long ib1Old = filteredBands[ik][0];
				ib1Old<filteredBands[ik][1]+1; ib1Old++ ) {
			for ( long i=0; i<3; i++ ) {
				groupVelocities.push_back(thisVelocity(ib1Old,ib1Old,i).real());
				for ( long ib2Old = filteredBands[ik][0];
						ib2Old<filteredBands[ik][1]+1; ib2Old++ ) {
					velocities.push_back(thisVelocity(ib1Old,ib2Old,i));
				}
			}
		}
	}
	return activePoints_;
}

ActivePoints ActiveBandStructure::buildAsPostprocessing(Window & window,
		FullBandStructure & fullBandStructure) {

	if ( fullBandStructure.hasEigenvectors ) {
		hasEigenvectors = true;
	}

	numAtoms = fullBandStructure.fullPoints.getCrystal().getNumAtoms();

	std::vector<long> filteredPoints;
	std::vector<std::vector<long>> filteredBands;
	// first we find the states that can be kept
	for ( long ik=0; ik<fullBandStructure.getNumPoints(); ik++ ) {
		Point point = fullBandStructure.getPoint(ik);
		State s = fullBandStructure.getState(point);
		auto theseEnergies = s.getEnergies();
		auto [ens, bandsExtrema] = window.apply(theseEnergies);
		if ( ens.empty() ) {
			continue;
		} else {
			filteredPoints.push_back(ik);
			filteredBands.push_back(bandsExtrema);

			for ( long ib=bandsExtrema[0]; ib<bandsExtrema[1]+1; ib++ ) {
				energies.push_back(ens[ib]);
			}

		}
	}

	// now we store the information we found on the states after the window
	// k/q-points in the active list
	numPoints = filteredPoints.size();
	// filter maps the new to old k index. filter(ik) = ikOld
	VectorXl filter(numPoints);
	for ( long i=0; i<numPoints; i++ ) {
		filter(i) = filteredPoints[i];
	}
	// store the number of bands per kpoint
	VectorXl numBands_(numPoints);
	for ( long ik=0; ik<numPoints; ik++ ) {
		numBands_(ik) = filteredBands[ik][1] - filteredBands[ik][0] + 1;
	}
	numBands = numBands_;
	// total number of active states
	numStates = numBands.sum();
	// initialize the kpoints object
	ActivePoints activePoints_(fullBandStructure.fullPoints, filter);
	activePoints = &activePoints_;
	// construct the mapping from combined indices to Bloch indices
	buildIndeces();

	// now we store the energies and other quantities in the reordered array

	// loop over new points
	long oldNumBands = fullBandStructure.numBands;

	for ( long ik=0; ik<numPoints; ik++ ) {
		long ikOld = filter(ik);
		long ib = 0;
		for ( long ibOld = filteredBands[ik][0]; ibOld<filteredBands[ik][1]+1;
				ibOld++ ) {

			for ( long i=0; i<3; i++) {
				long indOld = compressIndeces(ibOld, ibOld, i,
								oldNumBands, oldNumBands, 3);
				groupVelocities.push_back(
						fullBandStructure.velocities(ik,indOld).real());

				long ib2 = 0;
				for ( long ib2Old = filteredBands[ik][0];
						ib2Old<filteredBands[ik][1]+1; ib2Old++ ) {

					long indOld = compressIndeces(ibOld, ib2Old, i,
									oldNumBands, oldNumBands, 3);
					velocities.push_back(
							fullBandStructure.velocities(ikOld,indOld));
					ib2 += 1;
				}
			}
			ib += 1;
		}

		if ( hasEigenvectors ) {
			for ( long i=0; i<3; i++ ) {
				for ( long iat=0; iat<numAtoms; iat++ ) {
					long ib = 0;
					for ( long ibOld = filteredBands[ik][0];
							ibOld<filteredBands[ik][1]+1; ibOld++ ) {
						long indOld = compressIndeces(i, iat, ibOld, 3,
								numAtoms, oldNumBands);
						eigenvectors.push_back(
								fullBandStructure.eigenvectors(ikOld,indOld));
						ib += 1;
					}
				}
			}
		}
	}
	return activePoints_;
}

long ActiveBandStructure::velBloch2Comb(long & ik, long & ib1, long & ib2,
		long & i) {
	long is;
	long offset = cumulativeKbbOffset(ik);
	is = ik * offset + ib1 * numBands(ik) * 3 + ib2 * 3 + i;
	return is;
}

long ActiveBandStructure::gvelBloch2Comb(long & ik, long & ib, long & i) {
	long is;
	long offset = cumulativeKbOffset(ik)*3;
	is = ik * offset + ib * 3 + i;
	return is;
}

long ActiveBandStructure::eigBloch2Comb(long & ik, long & i, long & iat,
		long & ib) {
	long is;
	long offset = cumulativeKbOffset(ik);
	is = ik * offset + i * numBands(ik) * numAtoms + iat * numBands(ik) + ib;
	return is;
}

long ActiveBandStructure::bloch2Comb(long & ik, long & ib) {
	long is;
	long offset = cumulativeKbOffset(ik);
	is = ik * offset + ib;
	return is;
}

std::tuple<long,long> ActiveBandStructure::comb2Bloch(long & is) {
	long ik, ib;
	ik = auxBloch2Comb(is,0);
	ib = auxBloch2Comb(is,1);
	return {ik,ib};
}

void ActiveBandStructure::buildIndeces() {
	MatrixXl auxBloch2Comb_(numStates,2);
	VectorXl cumulativeKbOffset_(numPoints);
	VectorXl cumulativeKbbOffset_(numPoints);

	cumulativeKbOffset_(0) = 0;
	cumulativeKbbOffset_(0) = 0;
	for ( long ik=1; ik<numPoints; ik++ ) {
		cumulativeKbOffset_(ik) = cumulativeKbOffset_(ik-1) + numBands(ik-1);
		cumulativeKbbOffset_(ik) = cumulativeKbbOffset_(ik-1)
				+ 3 * numBands(ik-1) * numBands(ik-1);
	}

	long is = 0;
	for ( long ik=0; ik<numPoints; ik++ ) {
		for ( long ib=0; ib<numBands(ik); ib++ ) {
			auxBloch2Comb_(is,0) = ik;
			auxBloch2Comb_(is,1) = ib;
			is += 1;
		}
	}
	auxBloch2Comb = auxBloch2Comb_;
	cumulativeKbOffset = cumulativeKbOffset_;
	cumulativeKbbOffset = cumulativeKbbOffset_;
}

double ActiveBandStructure::getEnergy(long & stateIndex) {
	return energies[stateIndex];
}

Eigen::Vector3d ActiveBandStructure::getGroupVelocity(long & stateIndex) {
	Eigen::Vector3d vel;
	vel(0) = groupVelocities[stateIndex*3];
	vel(1) = groupVelocities[stateIndex*3+1];
	vel(2) = groupVelocities[stateIndex*3+2];
	return vel;
}
