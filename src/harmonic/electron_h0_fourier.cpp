#include <string>
#include <cmath>
#include "electron_h0_fourier.h"
#include "exceptions.h"
#include "constants.h"
#include "statistics.h"

ElectronH0Fourier::ElectronH0Fourier(Crystal & crystal_,
		FullPoints coarsePoints_, FullBandStructure coarseBandStructure_,
		double cutoff_) :
		crystal(crystal_), coarseBandStructure(coarseBandStructure_),
		coarsePoints(coarsePoints_) {
	numBands = coarseBandStructure.getNumBands();
	cutoff = cutoff_;
	if ( coarseBandStructure.hasIrreduciblePoints() ) {
		Error e("input electronic band structure must be specified on grid",1);
	}
	numDataPoints = coarseBandStructure.getNumPoints();
	refWavevector = coarseBandStructure.getPoint(0).getCoords("cartesian");

	// to do the interpolation, we need the lattice vector basis:
	setPositionVectors();

	// now we look for the expansion coefficients that interpolates the bands
	// note that setPositionVectors must stay above this call
	Eigen::VectorXcd expansionCoefficients_(numBands, numPositionVectors);
	Eigen::VectorXd energies(numDataPoints);
	for ( long iBand=0; iBand<numBands; iBand++ ) {
		energies = coarseBandStructure.getBandEnergies(iBand);
		expansionCoefficients_.row(iBand) = getCoefficients(energies);
	}
	expansionCoefficients = expansionCoefficients_;
}

double ElectronH0Fourier::getRoughnessFunction(Eigen::Vector3d position) {
	double norm = position.norm();
	return pow(1. - coeff1 * norm/minDistance,2)
			+ coeff2 * pow(norm/minDistance,6);
}

std::complex<double> ElectronH0Fourier::getStarFunction(
		Eigen::Vector3d& wavevector, long & iR) {
	std::complex<double> phase = complexI * twoPi
			* wavevector.transpose() * positionVectors[iR];
	std::complex<double> starFunction = exp(phase)
			/ double(positionDegeneracies[iR]);
	return starFunction;
}

Eigen::Vector3cd ElectronH0Fourier::getDerivativeStarFunction(
		Eigen::Vector3d& wavevector, long & iR) {
	std::complex<double> phase = complexI * twoPi * wavevector.transpose()
			* positionVectors[iR];
	Eigen::Vector3cd starFunctionDerivative = complexI * positionVectors[iR]
			* exp(phase) / double(positionDegeneracies[iR]);
	return starFunctionDerivative;
}

void ElectronH0Fourier::setPositionVectors() {
	// lattice parameters
	Eigen::Matrix3d directUnitCell = crystal.getDirectUnitCell();
	Eigen::Vector3d a1 = directUnitCell.row(0);
	Eigen::Vector3d a2 = directUnitCell.row(1);
	Eigen::Vector3d a3 = directUnitCell.row(2);

	std::vector<Eigen::Vector3d> tmpVectors;
	std::vector<bool> isDegenerate;
	std::vector<long> tmpDegeneracies;
	// the first vector is the vector zero
	tmpVectors.push_back(Eigen::Vector3d::Zero());

	long N0 = long(cutoff / a1.norm() * 2.);
	long N1 = long(cutoff / a2.norm() * 2.);
	long N2 = long(cutoff / a3.norm() * 2.);

	// we also need the smallest non zero vector norm
	minDistance = a1.norm();

	// build the list of positions
	for ( long i0=-N0; i0<=N0; i0++ ) {
		for ( long i1=-N1; i1<=N1; i1++ ) {
			for ( long i2=-N2; i2<=N2; i2++ ) {
				Eigen::Vector3d thisVec = (double)i0 * a1 + (double)i1 * a1
						+ (double)i2 * a2;
				if ( ( i0 != 0 ) && ( i1 != 0 ) && ( i2 != 0 ) ) {
					tmpVectors.push_back(thisVec);
					isDegenerate.push_back(false);
					tmpDegeneracies.push_back(1);

					double thisNorm = thisVec.norm();
					if ( thisNorm < minDistance ) {
						minDistance = thisNorm;
					}

				}
			}
		}
	}

	// now we build the list of irreducible vectors
	for ( long iR=0; iR<tmpVectors.size(); iR++ ) {
		Eigen::Vector3d vec = tmpVectors[iR];
		if ( ! isDegenerate[iR] ) { // then we look for reducible points
			for ( long iR2=iR+1; iR2<tmpVectors.size(); iR2++ ) {
				if ( ! isDegenerate[iR2] ) {
					for ( auto symmMatrix : crystal.getSymmetryMatrices() ) {
						Eigen::Vector3d rotVec = symmMatrix * tmpVectors[iR2];
						if ( (rotVec-vec).norm() < 1.0e-6 ) {
							isDegenerate[iR2] = true;
							tmpDegeneracies[iR] += 1;
						}
					}
				}
			}
		}
	}

	for ( long iR=0; iR<tmpVectors.size(); iR++ ) {
		if ( ! isDegenerate[iR] ) { // then we look for reducible points
			positionVectors.push_back(tmpVectors[iR]);
			positionDegeneracies.push_back(tmpDegeneracies[iR]);
		}
	}
	numPositionVectors = positionVectors.size();
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

	Eigen::Vector3d iWavevector, jWavevector;

	for ( long iR=0; iR<numPositionVectors; iR++ ) {
		Eigen::Vector3d position = positionVectors[iR];
		double rho = getRoughnessFunction(position);
		std::complex<double> smk0 = getStarFunction(refWavevector, iR);
		for ( long i=0; i<numDataPoints-1; i++ ) {
			iWavevector = coarseBandStructure.getPoint(i+1).getCoords("cartesian");
			std::complex<double> smki = getStarFunction(iWavevector, iR);
			for ( long j=0; j<numDataPoints-1; j++ ) {
				jWavevector = coarseBandStructure.getPoint(j+1).getCoords("cartesian");
				std::complex<double> smkj = getStarFunction(jWavevector, iR);
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
	Eigen::Vector3d wavevector;

	for ( long m=1; m<numPositionVectors; m++ ) {
		Eigen::Vector3d position = positionVectors[m];
		std::complex<double> smk0 = getStarFunction(refWavevector,m);
		for ( long i=1; i<numDataPoints; i++ ) {
			wavevector = coarseBandStructure.getPoint(i).getCoords("cartesian");
			coefficients(m) += multipliers(i)
					* std::conj(getStarFunction(wavevector,m) - smk0);
		}
		coefficients(m) /= getRoughnessFunction(position);
	}

	// special case for R=0
	coefficients(0) = energies(0);
	for ( long m=1; m<numPositionVectors; m++ ) {
		coefficients(0) -= coefficients(m)
				* getStarFunction(refWavevector, m);
	}

	return coefficients;
}

double ElectronH0Fourier::getEnergy(Point& point, long& bandIndex) {
	double energy = 0.;
	Eigen::Vector3d wavevector = point.getCoords("cartesian");
	for ( long m=0; m<numPositionVectors; m++ ) {
		std::complex<double> c = expansionCoefficients(bandIndex, m)
				* getStarFunction(wavevector, m);
		energy += c.real();
	}
	return energy;
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

std::tuple<Eigen::VectorXd, Eigen::Tensor<std::complex<double>,3>>
		ElectronH0Fourier::diagonalize(Point & point) {

	Eigen::Vector3d coords = point.getCoords("cartesian");
	auto [energies,x] = diagonalizeFromCoords(coords);

	// this is to return something aligned with the phonon case
	// One should investigate how to return a null pointer
	Eigen::Tensor<std::complex<double>,3> eigvecs;
	eigvecs.setZero();

	return {energies,eigvecs};
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


Eigen::Tensor<std::complex<double>,3> ElectronH0Fourier::diagonalizeVelocity(
			Point & point) {
	Eigen::Tensor<std::complex<double>,3> velocity(numBands,numBands,3);
	velocity.setZero();
	auto v = getGroupVelocities(point);
	for ( long ib=0; ib<numBands; ib++ ) {
		for ( long i=0; i<3; i++ ) {
			velocity(ib,ib,i) = v(ib,i);
		}
	}
}

Eigen::VectorXd ElectronH0Fourier::getEnergies(Point& point) {
	Eigen::VectorXd energies(numBands);
	for ( long bandIndex=0; bandIndex<numBands; bandIndex++ ) {
		energies(bandIndex) = getEnergy(point, bandIndex);
	}
	return energies;
}

Eigen::Vector3d ElectronH0Fourier::getGroupVelocity(Point& point,
		long& bandIndex) {
	Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
	Eigen::Vector3d wavevector = point.getCoords("cartesian");
	for ( long m=0; m<numPositionVectors; m++ ) {
		Eigen::Vector3cd c = expansionCoefficients(bandIndex, m)
				* getDerivativeStarFunction(wavevector, m);
		velocity += c.real();
	}
	return velocity;
}

Eigen::MatrixXd ElectronH0Fourier::getGroupVelocities(Point& point) {
	Eigen::MatrixXd velocities(numBands,3);
	for ( long bandIndex=0; bandIndex<numBands; bandIndex++ ) {
		velocities.row(bandIndex) = getGroupVelocity(point, bandIndex);
	}
	return velocities;
}

FullBandStructure ElectronH0Fourier::populateBandStructure(
		FullPoints* fullPoints, IrreduciblePoints* irreduciblePoints) {

	if ( ( fullPoints == nullptr ) && ( irreduciblePoints == nullptr ) ) {
		Error e("From populateBandStructure: must provide mesh of points", 1);
	}

	Statistics statistics(Statistics::fermi);
	bool withVelocities = false;
	bool withEigenvectors = false;
	FullBandStructure denseBandStructure(numBands, statistics,
			withVelocities, withEigenvectors, fullPoints, irreduciblePoints);

	Points * points;
	if ( fullPoints != nullptr ) {
		points = fullPoints;
	} else {
		points = irreduciblePoints;
	}

	Eigen::VectorXd energies(numBands);
	for ( long ik=0; ik<points->getNumPoints(); ik++ ) {
		Point point = points->getPoint(ik);
		energies = getEnergies(point);
		denseBandStructure.setEnergies(point, energies);
	}
	return denseBandStructure;
}
