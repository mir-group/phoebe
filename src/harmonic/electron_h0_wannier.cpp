#include "electron_h0_wannier.h"
#include "constants.h"
#include "exceptions.h"

ElectronH0Wannier::ElectronH0Wannier(const Eigen::Matrix3d & directUnitCell_,
		const Eigen::Matrix<double,3,Eigen::Dynamic> & bravaisVectors_,
		const Eigen::VectorXd & vectorsDegeneracies_,
		const Eigen::Tensor<std::complex<double>,3> & h0R_,
		const Eigen::Tensor<std::complex<double>,4> & rMatrix_) :
		statistics(Statistics::electron) {

	h0R = h0R_;
	rMatrix = rMatrix_;
	directUnitCell = directUnitCell_;
	bravaisVectors = bravaisVectors_;
	vectorsDegeneracies = vectorsDegeneracies_;

	if ( h0R.dimension(1) != h0R.dimension(2) ) {
		Error e("WannierH0(): h0R should have dimensions (R,bands,bands)");
	}
	if ( h0R.dimension(0) != bravaisVectors.cols() ) {
		Error e("WannierH0(): h0R and bravaisVectors not aligned");
	}
	if ( vectorsDegeneracies.size() != bravaisVectors.cols() ) {
		Error e("WannierH0(): degeneracies not aligned with vectors");
	}

	if ( ( rMatrix.dimension(1) != h0R.dimension(0) ) ||
			( rMatrix.dimension(2) != h0R.dimension(1) ) ||
			( rMatrix.dimension(3) != h0R.dimension(2) ) ) {
		Error e("WannierH0(): h0R and rMatrix should be aligned");
	}

	if ( rMatrix.dimension(0) != 3 ) {
		Error e("WannierH0(): rMatrix should be a vector");
	}

	numBands = h0R.dimension(1);
	numVectors = vectorsDegeneracies.size();
}

// copy constructor
ElectronH0Wannier::ElectronH0Wannier( const ElectronH0Wannier & that ) :
	statistics(Statistics::electron) {
		h0R = that.h0R;
		rMatrix = that.rMatrix;
		directUnitCell = that.directUnitCell;
		numBands = that.numBands;
		bravaisVectors = that.bravaisVectors;
		numVectors = that.numVectors;
		vectorsDegeneracies = that.vectorsDegeneracies;
}

// copy assignment
ElectronH0Wannier & ElectronH0Wannier::operator = (
		const ElectronH0Wannier & that ) {
	if ( this != & that ) {
	    bravaisVectors.resize(0,0);
	    vectorsDegeneracies.resize(0);
		h0R.resize(0,0,0);
		rMatrix.resize(0,0,0,0);
		statistics = that.statistics;
		numVectors = that.numVectors;
		numBands = that.numBands;
	    bravaisVectors = that.bravaisVectors;
	    vectorsDegeneracies = that.vectorsDegeneracies;
		directUnitCell = that.directUnitCell;
	    h0R = that.h0R;
	    rMatrix = that.rMatrix;
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

	for ( long iR=0; iR<bravaisVectors.cols(); iR++ ) {
		Eigen::Vector3d R = bravaisVectors.col(iR);
		double phase = k.dot(R);
		std::complex<double> phaseFactor = {cos(phase),sin(phase)};
		for ( long m=0; m<numBands; m++ ) {
			for ( long n=0; n<numBands; n++ ) {
				h0K(m,n) += phaseFactor *h0R(iR,m,n) /vectorsDegeneracies(iR);
			}
		}
	}

//	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(numBands);
//	eigensolver.compute(h0K);
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(h0K);
	Eigen::VectorXd energies = eigensolver.eigenvalues();
	Eigen::MatrixXcd eigenvectors = eigensolver.eigenvectors();

	return {energies, eigenvectors};
}
