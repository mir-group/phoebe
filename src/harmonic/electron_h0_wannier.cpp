#include "electron_h0_wannier.h"
#include "constants.h"
#include "exceptions.h"

ElectronH0Wannier::ElectronH0Wannier(Eigen::Matrix3d & directUnitCell_,
		Eigen::MatrixXd & crystalVectors_,
		Eigen::VectorXd & vectorsDegeneracies_,
		Eigen::Tensor<std::complex<double>,3> & h0R_) {

	crystalVectors = crystalVectors_;
	vectorsDegeneracies = vectorsDegeneracies_;
	directUnitCell = directUnitCell_;
	h0R = h0R_;

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

std::tuple<Eigen::VectorXd, Eigen::Tensor<std::complex<double>,3>>
ElectronH0Wannier::diagonalize(Point & point) {
	Eigen::Vector3d k = point.getCoords("cartesian");

	auto [energies,xxx] = diagonalizeFromCoords(k);

	// dummy eigenvector in output
	//TODO:  maybe we need to use these eigenvectors in a smarter way.
	Eigen::Tensor<std::complex<double>,3> eigvecs(1,1,1);
	eigvecs.setZero();

	return {energies, eigvecs};
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
