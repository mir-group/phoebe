#include <string>
#include <cmath>
#include "electron_h0_fourier.h"
#include "exceptions.h"
#include "constants.h"
#include "particle.h"

ElectronH0Fourier::ElectronH0Fourier(Crystal & crystal_,
		FullPoints coarsePoints_,
		FullBandStructure coarseBandStructure_, double cutoff_) :
		crystal(crystal_), coarseBandStructure(coarseBandStructure_),
		coarsePoints(coarsePoints_), particle(Particle::electron) {

	numBands = coarseBandStructure.getNumBands();
	cutoff = cutoff_;
//	if ( coarseBandStructure.hasIrreduciblePoints() ) {
//		Error e("input electronic band structure must be specified on grid",1);
//	}
	numDataPoints = coarseBandStructure.getNumPoints();
	refWavevector = coarseBandStructure.getPoint(0).getCoords(
			Points::cartesianCoords);

	// to do the interpolation, we need the lattice vector basis:
	setPositionVectors();
	// now we look for the expansion coefficients that interpolates the bands
	// note that setPositionVectors must stay above this call
	Eigen::MatrixXcd expansionCoefficients_(numBands, numPositionVectors);
	Eigen::VectorXd energies(numDataPoints);
	expansionCoefficients_.setZero();
	energies.setZero();
	std::cout << "Building coefficients" << std::endl;
	for ( long iBand=0; iBand<numBands; iBand++ ) {
		std::cout << iBand+1 << "/" << numBands << std::endl;
		energies = coarseBandStructure.getBandEnergies(iBand);
		expansionCoefficients_.row(iBand) = getCoefficients(energies);
	}
	std::cout << "Done building coefficients" << std::endl;
	expansionCoefficients = expansionCoefficients_;
}

Particle ElectronH0Fourier::getParticle() {
	return particle;
}

double ElectronH0Fourier::getRoughnessFunction(Eigen::Vector3d position) {
	double norm = position.norm();
	return pow(1. - coeff1 * norm / minDistance, 2)
			+ coeff2 * pow(norm/minDistance,6);
}

std::complex<double> ElectronH0Fourier::getStarFunction(
		Eigen::Vector3d & wavevector, long & iR) {
	std::complex<double> phase = complexI
			* wavevector.transpose() * positionVectors.col(iR);
	std::complex<double> starFunction = exp(phase) / positionDegeneracies(iR);
	return starFunction;
}

Eigen::Vector3cd ElectronH0Fourier::getDerivativeStarFunction(
		Eigen::Vector3d & wavevector, long & iR) {
	std::complex<double> phase = complexI
			* wavevector.transpose() * positionVectors.col(iR);
	Eigen::Vector3cd starFunctionDerivative = complexI
			* positionVectors.col(iR) * exp(phase) / positionDegeneracies(iR);
	return starFunctionDerivative;
}

void ElectronH0Fourier::setPositionVectors() {
	// lattice parameters
	Eigen::Matrix3d directUnitCell = crystal.getDirectUnitCell();
	Eigen::Vector3d a1 = directUnitCell.row(0);
	Eigen::Vector3d a2 = directUnitCell.row(1);
	Eigen::Vector3d a3 = directUnitCell.row(2);

	const long N0 = long(cutoff / a1.norm() * 6.);
	const long N1 = long(cutoff / a2.norm() * 6.);
	const long N2 = long(cutoff / a3.norm() * 6.);
	long nMax = N0*N1*N2;

	Eigen::MatrixXd tmpVectors(3,nMax);
	std::vector<bool> isDegenerate(nMax, false);
	Eigen::VectorXd tmpDegeneracies(nMax);
	// the first vector is the vector zero
	tmpVectors.setZero();
	tmpDegeneracies.setZero();
	tmpDegeneracies.array() += 1.;

	// we also need the smallest non zero vector norm
	minDistance = a1.norm(); // this is a first guess

	// build the list of positions within the cutoff sphere
	long numVec_ = 1;
	for ( long i0=-N0; i0<=N0; i0++ ) {
		for ( long i1=-N1; i1<=N1; i1++ ) {
			for ( long i2=-N2; i2<=N2; i2++ ) {
				Eigen::Vector3d thisVec = i0 * a1 + i1 * a2 + i2 * a3;
				double thisNorm = thisVec.norm();
				if ( ( thisNorm > epsilon8 ) && ( thisNorm < cutoff ) ) {
					tmpVectors.col(numVec_) = thisVec;
					numVec_ += 1;
					if ( thisNorm < minDistance ) {
						minDistance = thisNorm;
					}
				}
			}
		}
	}

	if ( numDataPoints >= numVec_ ) {
		Error e("The number of interpolating coefficients is smaller than data"
				" points: increase Fourier cutoff");
	}

	// this is wrong: only one symmetry recognized
//	std::cout << "Symm mats\n";
//	std::cout << crystal.getSymmetryMatrices().size() << "\n";
//	std::cout << crystal.getSymmetryMatrices()[0] << "\n";
//	std::cout << crystal.getSymmetryMatrices()[1] << "\n";
//	std::cout << crystal.getSymmetryMatrices()[2] << "\n";
//	std::cout << crystal.getSymmetryMatrices()[3] << "\n";
//	std::cout << crystal.getSymmetryMatrices()[4] << "\n";

//	 deactivate symmetries for debugging
	Eigen::Matrix3d symmMatrix;
	symmMatrix.setZero();
	symmMatrix(0,0) = 1.;
	symmMatrix(1,1) = 1.;
	symmMatrix(2,2) = 1.;

//	auto crystSymmMatrices = crystal.getSymmetryMatrices();
//	std::vector<Eigen::Matrix3d> symmMatrices;
//	for ( auto s : crystSymmMatrices ) {
//		Eigen::Matrix3d s2;
//		s2 = directUnitCell * s;
//		s2 = s2 * directUnitCell.inverse();
//		symmMatrices.push_back(s2);
//	}
//	numpy.dot(numpy.dot(a,sCryst),aInv)

	// now we build the list of irreducible vectors
	for ( long iR=1; iR<numVec_; iR++ ) { // skip first vec!
		Eigen::Vector3d vec = tmpVectors.col(iR);
		if ( ! isDegenerate[iR] ) { // then we look for reducible points
			for ( long iR2=iR+1; iR2<numVec_; iR2++ ) {
				if ( ! isDegenerate[iR2] ) {
//					for ( auto symmMatrix : crystal.getSymmetryMatrices() ) {
						Eigen::Vector3d rotVec =
								symmMatrix * tmpVectors.col(iR2);
						if ( (rotVec-vec).norm() < 1.0e-6 ) {
							isDegenerate[iR2] = true;
							tmpDegeneracies(iR) += 1;
						}
//					}
				}
			}
		}
	}

	numPositionVectors = 0;
	for ( long iR=0; iR<numVec_; iR++ ) {
		if ( ! isDegenerate[iR] ) {
			numPositionVectors += 1;
		}
	}

	Eigen::VectorXd positionDegeneracies_(numPositionVectors);
	Eigen::MatrixXd positionVectors_(3,numPositionVectors);

	long counter = 0;
	for ( long iR=0; iR<numVec_; iR++ ) {
		if ( ! isDegenerate[iR] ) { // then we look for reducible points
			positionVectors_.col(counter) = tmpVectors.col(iR);
			positionDegeneracies_(counter) = tmpDegeneracies(iR);
			counter += 1;
		}
	}

	positionVectors = positionVectors_;
	positionDegeneracies = positionDegeneracies_;
}

Eigen::VectorXcd ElectronH0Fourier::getLagrangeMultipliers(
		Eigen::VectorXd energies) {
	Eigen::VectorXcd multipliers(numDataPoints-1);
	Eigen::VectorXcd deltaEnergies(numDataPoints-1);
	Eigen::MatrixXcd h(numDataPoints-1,numDataPoints-1);
	h.setZero();

	for ( long i=1; i<numDataPoints; i++ ) {
		deltaEnergies(i-1) = energies(i) - energies(0);
	}

	Eigen::Vector3d iWavevector;
	std::complex<double> smki, smkj;
	Eigen::MatrixXcd smk(numDataPoints,numPositionVectors);
	for ( long iR=0; iR<numPositionVectors; iR++ ) {
		for ( long i=0; i<numDataPoints; i++ ) {
			iWavevector =
					coarseBandStructure.getPoint(i).getCoords(
							Points::cartesianCoords);
			std::complex<double> smki = getStarFunction(iWavevector, iR);
			smk(i,iR) = smki;
		}
	}

	for ( long iR=0; iR<numPositionVectors; iR++ ) {
		Eigen::Vector3d position = positionVectors.col(iR);
		double rho = getRoughnessFunction(position);
		std::complex<double> smk0 = getStarFunction(refWavevector, iR);
		for ( long i=0; i<numDataPoints-1; i++ ) {
			smki = smk(i+1,iR);
			for ( long j=0; j<numDataPoints-1; j++ ) {
				smkj = smk(j+1,iR);
				h(i,j) += (smki-smk0) * std::conj(smkj-smk0) / rho;
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

	for ( long m=1; m<numPositionVectors; m++ ) {
		std::complex<double> smk0 = getStarFunction(refWavevector,m);
		for ( long i=1; i<numDataPoints; i++ ) {
			wavevector = coarseBandStructure.getPoint(i
					).getCoords(Points::cartesianCoords);
			coefficients(m) += multipliers(i-1)
					* std::conj(getStarFunction(wavevector,m) - smk0);
		}
		Eigen::Vector3d position = positionVectors.col(m);
		coefficients(m) /= getRoughnessFunction(position);
	}

	// special case for R=0
	coefficients(0) = energies(0);
	for ( long m=1; m<numPositionVectors; m++ ) {
		coefficients(0) -= coefficients(m) * getStarFunction(refWavevector, m);
	}

	return coefficients;
}

double ElectronH0Fourier::getEnergyFromCoords(Eigen::Vector3d & wavevector,
		long & bandIndex) {
	double energy = 0.;
	for ( long m=0; m<numPositionVectors; m++ ) {
		std::complex<double> c = expansionCoefficients(bandIndex, m)
				* getStarFunction(wavevector, m);
		energy += c.real();
	}
	return energy;
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
		ElectronH0Fourier::diagonalizeFromCoords(Eigen::Vector3d & wavevector){
	Eigen::MatrixXcd eigvecs(numBands,numBands);
	eigvecs.setZero();
	Eigen::VectorXd energies(numBands);
	for ( long ib = 0; ib<numBands; ib++ ) {
		energies(ib) = getEnergyFromCoords(wavevector, ib);
	}
	return {energies,eigvecs};
}

Eigen::Vector3d ElectronH0Fourier::getGroupVelocityFromCoords(
		Eigen::Vector3d & wavevector, long & bandIndex) {
	Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
	for ( long m=0; m<numPositionVectors; m++ ) {
		Eigen::Vector3cd c = expansionCoefficients(bandIndex, m)
				* getDerivativeStarFunction(wavevector, m);
		velocity += c.real();
	}
	return velocity;
}

long ElectronH0Fourier::getNumBands() {
	return numBands;
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
		ElectronH0Fourier::diagonalize(Point & point) {

	Eigen::Vector3d coords = point.getCoords(Points::cartesianCoords);
	auto [energies,eigvecs] = diagonalizeFromCoords(coords);

	// this is to return something aligned with the phonon case
	// One should investigate how to return a null pointer
//	Eigen::MatrixXcdTensor<std::complex<double>,3> eigvecs;
//	eigvecs.setZero();

	return {energies,eigvecs};
}

Eigen::Tensor<std::complex<double>,3> ElectronH0Fourier::diagonalizeVelocity(
			Point & point) {
	Eigen::Tensor<std::complex<double>,3> velocity(numBands,numBands,3);
	velocity.setZero();
	Eigen::Vector3d coords = point.getCoords(Points::cartesianCoords);
	for ( long ib=0; ib<numBands; ib++ ) {
		Eigen::Vector3d v = getGroupVelocityFromCoords(coords,ib);
		for ( long i=0; i<3; i++ ) {
			velocity(ib,ib,i) = v(i);
		}
	}
	return velocity;
}

FullBandStructure ElectronH0Fourier::populate(Points & fullPoints,
		bool & withVelocities, bool & withEigenvectors) {

	FullBandStructure fullBandStructure(numBands, particle,
			withVelocities, withEigenvectors, fullPoints);

	for ( long ik=0; ik<fullBandStructure.getNumPoints(); ik++ ) {
		Point point = fullBandStructure.getPoint(ik);
		auto [ens, eigvecs] = diagonalize(point);
		fullBandStructure.setEnergies(point, ens);
		if ( withVelocities ) {
			auto vels = diagonalizeVelocity(point);
			fullBandStructure.setVelocities(point, vels);
		}
	}
	return fullBandStructure;
}
