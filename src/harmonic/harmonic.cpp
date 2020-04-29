#include "harmonic.h"
#include "exceptions.h"
#include "eigen.h"
#include "points.h"

// note: these are dummy functions, that should be overwritten in the

std::tuple<Eigen::VectorXd,Eigen::Tensor<std::complex<double>,3>>
	HarmonicHamiltonian::diagonalize(Point & point) {
  Eigen::VectorXd energies(1);
  Eigen::Tensor<std::complex<double>,3> eigvecs(1,1,1);
  energies.setZero();
  eigvecs.setZero();
  return {energies, eigvecs};
}

Eigen::VectorXd HarmonicHamiltonian::diagonalizeEnergy(Point & point) {
	auto [energies, eigvecs] = diagonalize(point);
	return energies;
}

std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
	HarmonicHamiltonian::diagonalizeFromCoords(Eigen::Vector3d & k) {
	Eigen::VectorXd energies(numBands);
	Eigen::MatrixXcd eigvecs(numBands,numBands);
	return {energies, eigvecs};
}

Eigen::Tensor<std::complex<double>,3> HarmonicHamiltonian::diagonalizeVelocity(
				Point & point) {
  Eigen::Tensor<std::complex<double>,3> c(1,1,1);
  c.setZero();
  return c;
}

Eigen::Tensor<std::complex<double>,3>
		HarmonicHamiltonian::internalDiagonalizeVelocity(
				Eigen::Vector3d & coords, double & delta, double & threshold) {


	Eigen::Tensor<std::complex<double>,3> velocity(numBands,numBands,3);
	velocity.setZero();

	// get the eigenvectors and the energies of the k-point
	auto [energies,eigenvectors] = diagonalizeFromCoords(coords);

	// now we compute the velocity operator, diagonalizing the expectation
	// value of the derivative of the dynamical matrix.
	// This works better than doing finite differences on the frequencies.
	for ( long i=0; i<3; i++ ) {
		// define q+ and q- from finite differences.
		Eigen::Vector3d qPlus = coords;
		Eigen::Vector3d qMins = coords;
		qPlus(i) += delta;
		qMins(i) -= delta;

		// diagonalize the dynamical matrix at q+ and q-
		auto [enPlus,eigPlus] = diagonalizeFromCoords(qPlus);
		auto [enMins,eigMins] = diagonalizeFromCoords(qMins);

		// build diagonal matrices with frequencies
		Eigen::MatrixXd enPlusMat(numBands,numBands);
		Eigen::MatrixXd enMinsMat(numBands,numBands);
		enPlusMat.diagonal() << enPlus;
		enMinsMat.diagonal() << enMins;

		// build the dynamical matrix at the two wavevectors
		// since we diagonalized it before, A = M.U.M*
		Eigen::MatrixXcd sqrtDPlus(numBands,numBands);
		sqrtDPlus = eigPlus * enPlusMat * eigPlus.adjoint();
		Eigen::MatrixXcd sqrtDMins(numBands,numBands);
		sqrtDMins = eigMins * enMinsMat * eigMins.adjoint();

		// now we can build the velocity operator
		Eigen::MatrixXcd der(numBands,numBands);
		der = (sqrtDPlus - sqrtDMins) / ( 2. * delta );

		// and to be safe, we reimpose hermiticity
		der = 0.5 * ( der + der.adjoint() );

		// now we rotate in the basis of the eigenvectors at q.
		der = eigenvectors.adjoint() * der * eigenvectors;

		for ( long ib1=0; ib1<numBands; ib1++ ) {
			for ( long ib2=0; ib2<numBands; ib2++ ) {
				velocity(ib1,ib2,i) = der(ib1,ib2);
			}
		}
	}

	// turns out that the above algorithm has problems with degenerate bands
	// so, we diagonalize the velocity operator in the degenerate subspace,

	for ( long ib=0; ib<numBands; ib++ ) {

		// first, we check if the band is degenerate, and the size of the
		// degenerate subspace
		long sizeSubspace = 1;
		for ( long ib2=ib+1; ib2<numBands; ib2++ ) {
			// I consider bands degenerate if their frequencies are the same
			// within 0.0001 cm^-1
			if ( abs(energies(ib)-energies(ib2)) > threshold ) {
				break;
			}
			sizeSubspace += 1;
		}

		if ( sizeSubspace > 1 ) {
			Eigen::MatrixXcd subMat(sizeSubspace,sizeSubspace);
			// we have to repeat for every direction
			for ( long iCart=0; iCart<3; iCart++ ) {

				// take the velocity matrix of the degenerate subspace
				for ( long i=0; i<sizeSubspace; i++ ) {
					for ( long j=0; j<sizeSubspace; j++ ) {
						subMat(i,j) = velocity(ib+i,ib+j,iCart);
					}
				}

				// reinforce hermiticity
				subMat = 0.5 * ( subMat + subMat.adjoint());

				// diagonalize the subMatrix
				Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>
						eigensolver(subMat);
				Eigen::MatrixXcd newEigvecs = eigensolver.eigenvectors();

				// rotate the original matrix in the new basis
				// that diagonalizes the subspace.
				subMat = newEigvecs.adjoint() * subMat * newEigvecs;

				// reinforce hermiticity
				subMat = 0.5 * ( subMat + subMat.adjoint());

				// substitute back
				for ( long i=0; i<sizeSubspace; i++ ) {
					for ( long j=0; j<sizeSubspace; j++ ) {
						velocity(ib+i,ib+j,iCart) = subMat(i,j);
					}
				}
			}
		}

		// we skip the bands in the subspace, since we corrected them already
		ib += sizeSubspace - 1;
	}
	return velocity;
}
