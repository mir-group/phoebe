#include "active_bandstructure.h"
#include "bandstructure.h"
#include "exceptions.h"
#include "window.h"

ActiveBandStructure::ActiveBandStructure(Particle & particle_) :
		particle(particle_) {
}

// copy constructor
ActiveBandStructure::ActiveBandStructure(const ActiveBandStructure & that) :
			particle(that.particle),
			energies(that.energies),
			velocities(that.velocities),
			eigenvectors(that.eigenvectors),
			activePoints(that.activePoints),
			hasEigenvectors(that.hasEigenvectors),
			numStates(that.numStates),
			numAtoms(that.numAtoms),
			numPoints(that.numPoints),
			numBands(that.numBands),
			numFullBands(that.numFullBands),
			auxBloch2Comb(that.auxBloch2Comb),
			cumulativeKbOffset(that.cumulativeKbOffset),
			cumulativeKbbOffset(that.cumulativeKbbOffset) {
}

ActiveBandStructure & ActiveBandStructure::operator=(
		const ActiveBandStructure & that) { // assignment operator
	if ( this != & that ) {
		particle = that.particle;
		energies = that.energies;
		velocities = that.velocities;
		eigenvectors = that.eigenvectors;
		activePoints = that.activePoints;
		hasEigenvectors = that.hasEigenvectors;
		numStates = that.numStates;
		numAtoms = that.numAtoms;
		numPoints = that.numPoints;
		numBands = that.numBands;
		numFullBands = that.numFullBands;
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
	if ( activePoints != nullptr ) {
		return true;
	} else {
		return false;
	}
}

Points ActiveBandStructure::getPoints() {
	if ( ! hasPoints() ) {
		Error e("ActiveBandStructure hasn't been populated yet" ,1);
	}
	return *activePoints;
}

Point ActiveBandStructure::getPoint(const long & pointIndex) {
	if ( ! hasPoints() ) {
		Error e("ActiveBandStructure hasn't been populated yet" ,1);
	}
	return activePoints->getPoint(pointIndex);
}

long ActiveBandStructure::getNumPoints() {
	if ( ! hasPoints() ) {
		Error e("ActiveBandStructure hasn't been populated yet" ,1);
	}
	return numPoints;
}

long ActiveBandStructure::getNumBands() {
	Error e("ActiveBandStructure doesn't have constant number of bands");
	return 0;
}

State ActiveBandStructure::getState(Point & point){
	if ( ! hasPoints() ) {
		Error e("ActiveBandStructure hasn't been populated yet" ,1);
	}

	long ik = point.getIndex();

	long ind = bloch2Comb(ik,0);
	double * thisEn = &energies[ind];

	ind = velBloch2Comb(ik,0,0,0);
	std::complex<double> * thisVel = &velocities[ind];

	std::complex<double> * thisEig = nullptr;
	if ( hasEigenvectors ) {
		ind = eigBloch2Comb(ik,0,0);
		thisEig = &eigenvectors[ind];
	}

	State s(point, thisEn, numAtoms, numBands(ik), thisVel,
			thisEig);
	return s;
}

long ActiveBandStructure::getIndex(const WavevectorIndex & ik,
		const BandIndex & ib) {
	if ( ! hasPoints() ) {
		Error e("ActiveBandStructure hasn't been populated yet" ,1);
	}
	return bloch2Comb(ik.get(), ib.get());
}

long ActiveBandStructure::getNumStates() {
	if ( ! hasPoints() ) {
		Error e("ActiveBandStructure hasn't been populated yet" ,1);
	}
	return numStates;
}

const double & ActiveBandStructure::getEnergy(const long & stateIndex) {
	if ( ! hasPoints() ) {
		Error e("ActiveBandStructure hasn't been populated yet" ,1);
	}
	return energies[stateIndex];
}

Eigen::Vector3d ActiveBandStructure::getGroupVelocity(const long & stateIndex){
	if ( ! hasPoints() ) {
		Error e("ActiveBandStructure hasn't been populated yet" ,1);
	}
	auto[ik,ib] = comb2Bloch(stateIndex);
	Eigen::Vector3d vel;
	vel(0) = velocities[velBloch2Comb(ik, ib, ib, 0)].real();
	vel(1) = velocities[velBloch2Comb(ik, ib, ib, 1)].real();
	vel(2) = velocities[velBloch2Comb(ik, ib, ib, 2)].real();
	return vel;
}

Eigen::Vector3d ActiveBandStructure::getWavevector(const long & stateIndex) {
	if ( ! hasPoints() ) {
		Error e("ActiveBandStructure hasn't been populated yet" ,1);
	}
	auto[ik,ib] = comb2Bloch(stateIndex);
	Point p = activePoints->getPoint(ik);
	return p.getCoords(Points::cartesianCoords);
}

void ActiveBandStructure::setEnergies(Point & point,
		Eigen::VectorXd & energies_) {
	long ik = point.getIndex();
	for ( long ib=0; ib<energies_.size(); ib++ ) {
		long index = bloch2Comb(ik,ib);
		energies[index] = energies_(ib);
	}
}

void ActiveBandStructure::setEigenvectors(Point & point,
		Eigen::MatrixXcd & eigenvectors_) {
	long ik = point.getIndex();
	for ( long i=0; i<eigenvectors_.rows(); i++ ) {
		for ( long j=0; j<eigenvectors_.cols(); j++ ) {
			long index = eigBloch2Comb(ik,i,j);
			eigenvectors[index] = eigenvectors_(i,j);
		}
	}
}

void ActiveBandStructure::setVelocities(Point & point,
		Eigen::Tensor<std::complex<double>,3> & velocities_) {
	long ik = point.getIndex();
	for ( long ib1=0; ib1<velocities_.dimension(0); ib1++ ) {
		for ( long ib2=0; ib2<velocities_.dimension(1); ib2++ ) {
			for ( long j : {0,1,2} ) {
				long index = velBloch2Comb(ik,ib1,ib2,j);
				velocities[index] = velocities_(ib1,ib2,j);
			}
		}
	}
}

ActivePoints ActiveBandStructure::buildAsPostprocessing(Window & window,
		FullBandStructure & fullBandStructure) {

	if ( fullBandStructure.hasEigenvectors ) {
		hasEigenvectors = true;
	}

	numAtoms = fullBandStructure.points.getCrystal().getNumAtoms();

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
	auto fullPoints = fullBandStructure.getPoints();
	ActivePoints activePoints_(fullPoints, filter);
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
				long indOld = compress3Indeces(ibOld, ibOld, i,
								oldNumBands, oldNumBands, 3);
				long ib2 = 0;
				for ( long ib2Old = filteredBands[ik][0];
						ib2Old<filteredBands[ik][1]+1; ib2Old++ ) {

					long indOld = compress3Indeces(ibOld, ib2Old, i,
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
						long indOld = compress3Indeces(i, iat, ibOld, 3,
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

long ActiveBandStructure::velBloch2Comb(const long & ik, const long & ib1,
		const long & ib2, const long & i) {
	return ik * cumulativeKbbOffset(ik) + ib1 * numBands(ik) * 3 + ib2 * 3 + i;
}

long ActiveBandStructure::eigBloch2Comb(const long & ik, const long & ib1,
		const long & ib2) {
	return ik * cumulativeKbOffset(ik) + ib1 * numFullBands + ib2;
}

long ActiveBandStructure::bloch2Comb(const long & ik, const long & ib) {
	return ik * cumulativeKbOffset(ik) + ib;
}

std::tuple<long,long> ActiveBandStructure::comb2Bloch(const long & is) {
	return {auxBloch2Comb(is,0),auxBloch2Comb(is,1)};
}

void ActiveBandStructure::buildIndeces() {
	MatrixXl auxBloch2Comb_(numStates,2);
	VectorXl cumulativeKbOffset_(numPoints);
	VectorXl cumulativeKbbOffset_(numPoints);

	cumulativeKbOffset_(0) = 0;
	cumulativeKbbOffset_(0) = 0;
	for ( long ik=1; ik<numPoints; ik++ ) {
		cumulativeKbOffset_(ik) = cumulativeKbOffset_(ik-1)
				+ numBands(ik-1) * numFullBands;
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
