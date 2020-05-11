#include "electron_h0_wannier.h"
#include "constants.h"
#include "exceptions.h"

ElectronH0Wannier::ElectronH0Wannier(const Eigen::Matrix3d & directUnitCell_,
		const Eigen::MatrixXd & crystalVectors_,
		const Eigen::VectorXd & vectorsDegeneracies_,
		const Eigen::Tensor<std::complex<double>,3> & h0R_) :
		statistics(Statistics::electron) {

	h0R = h0R_;
	directUnitCell = directUnitCell_;
	crystalVectors = crystalVectors_;
	vectorsDegeneracies = vectorsDegeneracies_;

	if ( crystalVectors.cols() != 3 ) {
		Error e("WannierH0(): crystalVectors should have dimensions (R,3)", 1);
	}
	if ( h0R.dimension(1) != h0R.dimension(2) ) {
		Error e("WannierH0(): h0R should have dimensions (R,bands,bands)", 1);
	}
	if ( h0R.dimension(0) != crystalVectors.rows() ) {
		Error e("WannierH0(): h0R and crystalVectors not aligned", 1);
	}
	if ( vectorsDegeneracies.size() != crystalVectors.rows() ) {
		Error e("WannierH0(): degeneracies not aligned with vectors", 1);
	}
	numBands = h0R.dimension(1);
	numVectors = vectorsDegeneracies.size();
}

// copy constructor
ElectronH0Wannier::ElectronH0Wannier( const ElectronH0Wannier & that ) :
	statistics(Statistics::electron) {
		h0R = that.h0R;
		directUnitCell = that.directUnitCell;
		numBands = that.numBands;
		crystalVectors = that.crystalVectors;
		numVectors = that.numVectors;
		vectorsDegeneracies = that.vectorsDegeneracies;
}

// copy assignment
ElectronH0Wannier & ElectronH0Wannier::operator = (
		const ElectronH0Wannier & that ) {
	if ( this != & that ) {
	    crystalVectors.resize(0,0);
	    vectorsDegeneracies.resize(0);
		h0R.resize(0,0,0);
		statistics = that.statistics;
		numVectors = that.numVectors;
		numBands = that.numBands;
	    crystalVectors = that.crystalVectors;
	    vectorsDegeneracies = that.vectorsDegeneracies;
		directUnitCell = that.directUnitCell;
	    h0R = that.h0R;
	}
	return *this;
}

// default constructor
ElectronH0Wannier::ElectronH0Wannier() : statistics(Statistics::electron) {
}

long ElectronH0Wannier::getNumBands() {
	return numBands;
}

Statistics ElectronH0Wannier::getStatistics() {
	return statistics;
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
	ElectronH0Wannier::diagonalizeFromCoords(Eigen::Vector3d & k) {

	Eigen::MatrixXcd h0K(numBands,numBands);
	h0K.setZero();

	for ( long iR=0; iR<crystalVectors.rows(); iR++ ) {
		Eigen::Vector3d R = crystalVectors.row(iR);
		double phase = twoPi * k.transpose() * R;
		std::complex<double> phaseFactor = {cos(phase),sin(phase)};
		for ( long m=0; m<numBands; m++ ) {
			for ( long n=0; n<numBands; n++ ) {
				h0K(m,n) += phaseFactor *h0R(iR,m,n) /vectorsDegeneracies(iR);
			}
		}
	}

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(h0K);
	Eigen::VectorXd energies = eigensolver.eigenvalues();
	Eigen::MatrixXcd eigenvectors = eigensolver.eigenvectors();

	return {energies, eigenvectors};
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
		ElectronH0Wannier::diagonalize(Point & point) {
	Eigen::Vector3d k = point.getCoords("cartesian");

	auto [energies,eigenvectors] = diagonalizeFromCoords(k);

	// note: the eigenvector matrix is the unitary transformation matrix U
	// from the Bloch to the Wannier gauge.

	return {energies, eigenvectors};
}

Eigen::Tensor<std::complex<double>,3> ElectronH0Wannier::diagonalizeVelocity(
		Point & point) {
	Eigen::Vector3d coords = point.getCoords("cartesian");
	double delta = 1.0e-8;
	double threshold = 0.000001 / energyRyToEv; // = 1 micro-eV
	auto velocity = HarmonicHamiltonian::internalDiagonalizeVelocity(coords,
			delta, threshold);
	return velocity;
}
