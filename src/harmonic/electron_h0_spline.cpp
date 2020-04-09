#include <string>
#include <cmath>
#include <unsupported/Eigen/Splines>
#include "electron_h0_spline.h"
#include "exceptions.h"
#include "constants.h"

ElectronH0Spline::ElectronH0Spline(Crystal& crystal_,
		ElBandStructure coarseBandStructure_, double cutoff_) :
		crystal{crystal_}, coarseBandStructure{coarseBandStructure_} {
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
	for ( int iBand=0; iBand<numBands; iBand++ ) {
		energies = coarseBandStructure.getBandEnergies(iBand);
		expansionCoefficients_.row(iBand) = getCoefficients(energies);
	}
	expansionCoefficients = expansionCoefficients_;
}

double ElectronH0Spline::getRoughnessFunction(Eigen::Vector3d position) {
	double norm = position.norm();
	return pow(1. - coeff1 * norm/minDistance,2)
			+ coeff2 * pow(norm/minDistance,6);
}

std::complex<double> ElectronH0Spline::getStarFunction(
		Eigen::Vector3d& wavevector, Eigen::Vector3d& position) {
	std::complex<double> starFunction, phase;
	starFunction = complexZero;

	for ( auto symmMatrix : crystal.getSymmetryMatrices() ) {
		phase = complexI * wavevector.transpose() * symmMatrix * position;
		starFunction += exp(phase);
	}
	starFunction /= (double)crystal.getNumSymmetries();

	return starFunction;
}

Eigen::Vector3cd ElectronH0Spline::getDerivativeStarFunction(
		Eigen::Vector3d& wavevector, Eigen::Vector3d& position) {
	std::complex<double> phase;
	Eigen::Vector3cd starFunctionDerivative;
	starFunctionDerivative.setZero();

	for ( auto symmMatrix : crystal.getSymmetryMatrices() ) {
		phase = complexI * wavevector.transpose() * symmMatrix * position;
		starFunctionDerivative += complexI * position * exp(phase);
	}
	starFunctionDerivative /= (double)crystal.getNumSymmetries();
	return starFunctionDerivative;
}

void ElectronH0Spline::setPositionVectors() {
	// lattice parameters
	Eigen::Matrix3d directUnitCell = crystal.getDirectUnitCell();
	Eigen::Vector3d a1 = directUnitCell.row(0);
	Eigen::Vector3d a2 = directUnitCell.row(1);
	Eigen::Vector3d a3 = directUnitCell.row(2);

	std::vector<Eigen::Vector3d> tmpVectors;
	// the first vector is the vector zero
	Eigen::Vector3d thisVec;
	thisVec.setZero();
	positionVectors.push_back(thisVec);

	int N0 = int(cutoff / a1.norm() * 2.);
	int N1 = int(cutoff / a2.norm() * 2.);
	int N2 = int(cutoff / a3.norm() * 2.);

	// we also need the smallest non zero vector norm
	double thisNorm;
	minDistance = a1.norm();

	// build the list of positions
	for ( int i0=-N0; i0<=N0; i0++ ) {
		for ( int i1=-N1; i1<=N1; i1++ ) {
			for ( int i2=-N2; i2<=N2; i2++ ) {
				thisVec = (double)i0 * a1 + (double)i1 * a1 + (double)i2 * a2;
				if ( ( i0 != 0 ) && ( i1 != 0 ) && ( i2 != 0 ) ) {
					positionVectors.push_back(thisVec);

					thisNorm = thisVec.norm();
					if ( thisNorm < minDistance ) {
						minDistance = thisNorm;
					}

				}
			}
		}
	}
	numPositionVectors = positionVectors.size();
}

Eigen::VectorXcd ElectronH0Spline::getLagrangeMultipliers(
		Eigen::VectorXd energies) {
	Eigen::VectorXcd multipliers(numDataPoints-1);
	Eigen::VectorXcd deltaEnergies(numDataPoints-1);
	Eigen::MatrixXcd h(numDataPoints-1,numDataPoints-1);
	h.setZero();

	for ( int i=1; i<numDataPoints; i++ ) {
		deltaEnergies(i-1) = energies(i) - energies(0);
	}

	Eigen::Vector3d iWavevector, jWavevector;
	std::complex<double> smki, smkj, smk0;
	double rho;

	for ( auto position : positionVectors ) {
		rho = getRoughnessFunction(position);
		smk0 = getStarFunction(refWavevector, position);
		for ( int i=0; i<numDataPoints-1; i++ ) {
			iWavevector = coarseBandStructure.getPoint(i+1).getCoords("cartesian");
			smki = getStarFunction(iWavevector, position);
			for ( int j=0; j<numDataPoints-1; j++ ) {
				jWavevector = coarseBandStructure.getPoint(j+1).getCoords("cartesian");
				smkj = getStarFunction(jWavevector, position);
				h(i,j) += (smki-smk0) * std::conj(smkj-smk0) / rho;
			}
		}
	}

	// solve h*multipliers = deltaEnergies
	multipliers = h.ldlt().solve(deltaEnergies);
	return multipliers;
}

Eigen::VectorXcd ElectronH0Spline::getCoefficients(Eigen::VectorXd energies) {
	Eigen::VectorXcd multipliers = getLagrangeMultipliers(energies);

	Eigen::VectorXcd coefficients(numPositionVectors);
	Eigen::Vector3d position, wavevector;
	std::complex<double> smk0;

	for ( int m=1; m<numPositionVectors; m++ ) {
		position = positionVectors[m];
		smk0 = getStarFunction(refWavevector,position);
		for ( int i=1; i<numDataPoints; i++ ) {
			wavevector = coarseBandStructure.getPoint(i).getCoords("cartesian");
			coefficients(m) += multipliers(i)
					* std::conj(getStarFunction(wavevector,position) - smk0);
		}
		coefficients(m) /= getRoughnessFunction(position);
	}

	// special case for R=0
	coefficients(0) = energies(0);
	for ( int m=1; m<numPositionVectors; m++ ) {
		coefficients(0) -= coefficients(m)
				* getStarFunction(refWavevector, position);
	}

	return coefficients;
}

double ElectronH0Spline::getEnergy(Point& point, int& bandIndex) {
	double energy = 0.;
	std::complex<double> c;
	Eigen::Vector3d wavevector = point.getCoords("cartesian");
	for ( int m=0; m<numPositionVectors; m++ ) {
		c = expansionCoefficients(bandIndex, m)
				* getStarFunction(wavevector, positionVectors[m]);
		energy += c.real();
	}
	return energy;
}

Eigen::VectorXd ElectronH0Spline::getEnergies(Point& point) {
	Eigen::VectorXd energies(numBands);
	for ( int bandIndex=0; bandIndex<numBands; bandIndex++ ) {
		energies(bandIndex) = getEnergy(point, bandIndex);
	}
	return energies;
}

Eigen::Vector3d ElectronH0Spline::getGroupVelocity(Point& point,
		int& bandIndex) {
	Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
	Eigen::Vector3d wavevector = point.getCoords("cartesian");
	Eigen::Vector3cd c;
	for ( int m=0; m<numPositionVectors; m++ ) {
		c = expansionCoefficients(bandIndex, m)
				* getDerivativeStarFunction(wavevector, positionVectors[m]);
		velocity += c.real();
	}
	return velocity;
}

Eigen::MatrixXd ElectronH0Spline::getGroupVelocities(Point& point) {
	Eigen::MatrixXd velocities(numBands,3);
//	velocities.setZero();
	for ( int bandIndex=0; bandIndex<numBands; bandIndex++ ) {
		velocities.row(bandIndex) = getGroupVelocity(point, bandIndex);
	}
	return velocities;
}

ElBandStructure ElectronH0Spline::populateBandStructure(
		FullPoints* fullPoints, IrreduciblePoints* irreduciblePoints) {

	if ( ( fullPoints == nullptr ) && ( irreduciblePoints == nullptr ) ) {
		Error e("From populateBandStructure: must provide mesh of points", 1);
	}

	Points* points;
	ElBandStructure denseBandStructure(numBands, fullPoints, irreduciblePoints);
	if ( fullPoints != nullptr ) {
		points = fullPoints;
	} else {
		points = irreduciblePoints;
	}

	Eigen::VectorXd energies(numBands);
	for ( int ik=0; ik<points->getNumPoints(); ik++ ) {
		Point point = points->getPoint(ik);
		energies = getEnergies(point);
		denseBandStructure.setEnergies(point, energies);
	}
	return denseBandStructure;
}
