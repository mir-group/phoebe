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

//Eigen::Tensor<std::complex<double>,3> unpackStdXcd(std::vector<std::complex<double>> b,
//		long size1, long size2, long size3) {
//	Eigen::Tensor<std::complex<double>,3> a(size1,size2,size3);
//	long ind;
//	for ( long i=0; i<size1; i++ ) {
//		for ( long j=0; j<size2; j++ ) {
//			for ( long k=0; k<size3; k++ ) {
//				ind = compressIndeces(i, j, k, size1, size2, size3);
//				a(i,j,k) = b[ind];
//			}
//		}
//	}
//	return a;
//}

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

//Eigen::Tensor<std::complex<double>,3> unpackXcd(Eigen::VectorXcd b,
//		long size1, long size2, long size3) {
//	Eigen::Tensor<std::complex<double>,3> a(size1,size2,size3);
//	long ind;
//	for ( long i=0; i<size1; i++ ) {
//		for ( long j=0; j<size2; j++ ) {
//			for ( long k=0; k<size3; k++ ) {
//				ind = compressIndeces(i, j, k, size1, size2, size3);
//				a(i,j,k) = b(ind);
//			}
//		}
//	}
//	return a;
//}

FullBandStructure::FullBandStructure(long numBands_, Statistics & statistics_,
		bool withVelocities, bool withEigenvectors,
		FullPoints * fullPoints_, IrreduciblePoints * irreduciblePoints_) :
			statistics{statistics_} {

	numBands = numBands_;
	numAtoms = numBands_ / 3;

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

	if ( fullPoints != nullptr ) {
		useIrreducible = false;
		fullPoints = fullPoints_;
	} else if ( irreduciblePoints != nullptr ) {
		useIrreducible = true;
		irreduciblePoints = irreduciblePoints_;
	} else {
		Error e("FullBandStructure must provide one Points mesh.", 1);
	}
	if ( ( fullPoints != nullptr ) && ( irreduciblePoints != nullptr ) ) {
		Error e("FullBandStructure must provide only one Points mesh.", 1);
	}

	Eigen::MatrixXd energies_(getNumPoints(), numBands);
	energies.setZero();
	energies = energies_;

//	Eigen::MatrixXd dndt_(getNumPoints(), numBands);
//	Eigen::MatrixXd dnde_(getNumPoints(), numBands);
//	dndt_.setZero();
//	dnde_.setZero();
//	dndt = dndt_;
//	dnde = dnde_;

	// here is some black magic
	// http://dovgalecs.com/blog/eigen-how-to-get-in-and-out-data-from-eigen-matrix/
	// we want to access the raw data of Eigen Matrix without moving it
	// to do so, we need to use the map functionality of Eigen
	Eigen::Map<Eigen::MatrixXd> (rawEnergies, energies.rows(), energies.cols() ) = energies;
	// After the call to map, rawEnergies is a pointer to the beginning
	// of the raw data. Note that by default, Eigen stores in column major
	// so we had to specify that we want row-major order.
	if ( hasVelocities ) {
		Eigen::Map<Eigen::MatrixXcd> (rawVelocities, velocities.rows(), velocities.cols()) = velocities;
	}
	if ( hasEigenvectors ) {
		Eigen::Map<Eigen::MatrixXcd> (rawEigenvectors, eigenvectors.rows(), eigenvectors.cols()) = eigenvectors;
	}

	energiesRows = numBands;
	velocitiesRows = numBands*numBands*3;
	eigenvectorsRows = numBands * numAtoms * 3;
}

long FullBandStructure::getNumBands() {
	return numBands;
}

bool FullBandStructure::hasIrreduciblePoints() {
	return useIrreducible;
}

long FullBandStructure::getNumPoints() {
	if ( useIrreducible ) {
		return irreduciblePoints->getNumPoints();
	} else {
		return fullPoints->getNumPoints();
	}
}

long FullBandStructure::getIndex(Eigen::Vector3d& pointCoords) {
	if ( useIrreducible ) {
		return irreduciblePoints->getIndex(pointCoords);
	} else {
		return fullPoints->getIndex(pointCoords);
	}
}

Point FullBandStructure::getPoint(const long& pointIndex) {
	if ( useIrreducible ) {
		return irreduciblePoints->getPoint(pointIndex);
	} else {
		return fullPoints->getPoint(pointIndex);
	}
}

//void FullBandStructure::populate() {
//	Error e("populate() not implemented in FullBandStructure", 1);
//}

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

////TODO: maybe combine set chemPot with set Temp?
//void FullBandStructure::setChemicalPotential(double chemPot) {
//	chemicalPotential = chemPot;
//	FullBandStructure::setOccupations();
//}
//
//void FullBandStructure::setTemperature(double temp) {
//	temperature = temp;
//	FullBandStructure::setOccupations();
//}

//void FullBandStructure::setNumValenceElectrons(long numElectrons) {
//	numValenceElectrons = numElectrons;
//}
//void FullBandStructure::setHomo(double homo_) {
//	homo = homo_;
//}

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

void FullBandStructure::populate(HarmonicHamiltonian & h0) {
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
	long bandIndexMin = filteredBands[ik][0];
	long bandIndexMax = filteredBands[ik][1];
	long numBands = bandIndexMax - bandIndexMin + 1;

	double * thisEn = &energies[ik][0];

	std::complex<double> * thisVel;
	thisVel = &velocities[ik][0];

	std::complex<double> * thisEig = nullptr;

	if ( hasEigenvectors ) {
		thisEig = &eigenvectors[ik][0];
	}

	State s(point, thisEn, numAtoms, numBands, thisVel, thisEig);
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
	std::vector<std::vector<double>> filteredEnergies;
	std::vector<std::vector<std::complex<double>>> filteredEigenvectors;

	numStates = 0;

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
			filteredEnergies.push_back(ens);

			int numTheseBands = bandsExtrema[1] - bandsExtrema[0] + 1;
			numStates += numTheseBands;

			if ( h0.hasEigenvectors ) {
				Eigen::Tensor<std::complex<double>,3> eigvec(3,numAtoms,
						numTheseBands);

				for ( int i=0; i<3; i++ ) {
					for ( int iat=0; iat<numAtoms; iat++ ) {
						for ( int ib=0; ib<numTheseBands; ib++ ) {
							eigvec(i,iat,ib) = theseEigenvectors(i,iat,
									ib+bandsExtrema[0]);
						}
					}
				}

				filteredEigenvectors.push_back(
						packStdXcd(eigvec, 3, numAtoms, numTheseBands));
			}
		}

	}

	numPoints = filteredPoints.size();
	VectorXl filter(numPoints);
	for ( long i=0; i<numPoints; i++ ) {
		filter(i) = filteredPoints[i];
	}

	ActivePoints activePoints_(fullPoints, filter);
	activePoints = &activePoints_;

	energies = filteredEnergies;

	if ( h0.hasEigenvectors ) {
		hasEigenvectors = true;
		eigenvectors = filteredEigenvectors;
	}

	// now we add velocities
	std::vector<std::vector<std::complex<double>>> velocities;
	for ( long ik=0; ik<numPoints; ik++ ) {
		Point point = activePoints_.getPoint(ik);

		// thisVelocity is a tensor of dimensions (ib, ib, 3)
		auto thisVelocity = h0.diagonalizeVelocity(point);

		// now we filter it

		int numTheseBands = filteredBands[ik][1] - filteredBands[ik][0] + 1;
		Eigen::Tensor<std::complex<double>,3> vel(numTheseBands,numTheseBands,3);
		vel.setZero();

		for ( int ib1 = 0; ib1<numTheseBands; ib1++ ) {
			for ( int ib2 = 0; ib2<numTheseBands; ib2++ ) {
				for ( int i=0; i<3; i++ ) {
					vel(ib1,ib2,i) =
							thisVelocity(ib1+filteredBands[ik][0],
									ib2+filteredBands[ik][0],i);
				}
			}
		}

		// reshape and store it
		velocities.push_back(packStdXcd(vel, numTheseBands, numTheseBands, 3));
	}
	return activePoints_;
}

ActivePoints ActiveBandStructure::buildAsPostprocessing(Window & window,
		FullBandStructure & fullBandStructure) {

	if ( fullBandStructure.hasEigenvectors ) {
		hasEigenvectors = true;
	}

	numAtoms = fullBandStructure.fullPoints->getCrystal().getNumAtoms();

	std::vector<long> filteredPoints;
	std::vector<std::vector<long>> filteredBands;

	numStates = 0;
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
			energies.push_back(ens);

			long numTheseBands = bandsExtrema[1] - bandsExtrema[0] + 1;
			numStates += numTheseBands;

			if ( hasEigenvectors ) {
				auto theseEigenvectors = s.getEigenvectors();

				Eigen::Tensor<std::complex<double>,3> eigvec(3,numAtoms,
						numTheseBands);

				for ( long i=0; i<3; i++ ) {
					for ( long iat=0; iat<numAtoms; iat++ ) {
						for ( long ib=0; ib<numTheseBands; ib++ ) {
							eigvec(i,iat,ib) = theseEigenvectors(i,iat,
									ib+bandsExtrema[0]);
						}
					}
				}
				eigenvectors.push_back(
						packStdXcd(eigvec, 3, numAtoms, numTheseBands));
			}

			auto theseVelocities = s.getVelocities();
			Eigen::Tensor<std::complex<double>,3> vel(numTheseBands,
					numTheseBands, 3);

			for ( int ib1=0; ib1<numTheseBands; ib1++ ) {
				for ( int ib2=0; ib2<numTheseBands; ib2++ ) {
					for ( int i=0; i<3; i++ ) {
						vel(ib1,ib2,i) = theseVelocities(ib1+bandsExtrema[0],
								ib2+bandsExtrema[0],i);
					}
				}
			}
			velocities.push_back(
					packStdXcd(vel, numTheseBands, numTheseBands, 3));
		}
	}

	numPoints = filteredPoints.size();
	VectorXl filter(numPoints);
	for ( long i=0; i<numPoints; i++ ) {
		filter(i) = filteredPoints[i];
	}

	ActivePoints activePoints_(*fullBandStructure.fullPoints, filter);
	activePoints = &activePoints_;
	return activePoints_;
}
