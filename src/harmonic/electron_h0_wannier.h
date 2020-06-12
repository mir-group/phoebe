#ifndef WANNIERH0_H
#define WANNIERH0_H

#include <math.h>
#include "eigen.h"
#include "harmonic.h"
#include "points.h"
#include "bandstructure.h"
#include "constants.h"

/** Class for diagonalizing electronic energies with the Wannier interpolation
 * The object is built passing the information produced by the file _tb.dat of
 * Wannier90, and can be used to compute the interpolated band structure.
 */
class ElectronH0Wannier : public HarmonicHamiltonian {
public:
	ElectronH0Wannier(const Eigen::Matrix3d & directUnitCell_,
			const Eigen::Matrix<double,3,Eigen::Dynamic> & bravaisVectors_,
			const Eigen::VectorXd & vectorsDegeneracies_,
			const Eigen::Tensor<std::complex<double>,3> & h0R_,
			const Eigen::Tensor<std::complex<double>,4> & rMatrix_);

	/** get the electronic energies (in Ry) at a single k-point.
	 * Energies don't have any reference value, and must be used in connection
	 * with a chemical potential.
	 * @param k: a point object with the wavevector. Must have the cartesian
	 * coordinates of the wavevector.
	 * @return tuple(energies, eigenvectors): the energies are a double vector
	 * of size (numBands). Eigenvectors of size (numBands,numBands) are the
	 * complex unitary transformation matrix U, connecting the Wannier gauge
	 * with the Bloch gauge.
	 */
	template<typename T>
	std::tuple<Eigen::VectorXd,Eigen::MatrixXcd> diagonalize(
			Point<T> & point);

	/** get the electron velocities (in atomic units) at a single k-point.
	 * @param k: a Point object with the wavevector coordinates.
	 * @return velocity(numBands,numBands,3): values of the velocity operator
	 * for this state, in atomic units.
	 */
	template<typename T>
	Eigen::Tensor<std::complex<double>,3> diagonalizeVelocity(Point<T> &point);

	// checks whether this object can compute eigenvectors.
    const bool hasEigenvectors = true;

    /** Method to return that the underlying is that of an electronic Fermion.
     */
    Particle getParticle();

    /** get the total number of bands.
     * This is a constant for all wavevectors.
     */
    long getNumBands();

    /** Copy constructor
     */
    ElectronH0Wannier( const ElectronH0Wannier & that );

    /** Copy assignment
     */
    ElectronH0Wannier & operator = ( const ElectronH0Wannier & that );

	/** This method constructs an electron bandstructure.
	 * @param points: the object with the list/mesh of wavevectors
	 * @param withVelocities: if true, compute the electron velocity operator.
	 * @param withEigenvectors: if true, stores the Wannier eigenvectors.
	 * @return FullBandStructure: the bandstructure object containing the
	 * complete electronic band structure.
	 */
    template<typename T>
    FullBandStructure<T> populate(T & fullPoints, bool & withVelocities,
    		bool &withEigenvectors);

    /** compute the Berry connection <u_mk| nabla_k |u_nk> at arb. wavevectors.
     * @param point: the Point coordinates of the wavevector.
     * @return Berry connection: a generalized Berry connection in the form of
     * a matrix <u_mk| nabla_k |u_nk> for a fixed wavevector. The Berry
     * connection is actually just the diagonal matrix elements.
     */
    template<typename T>
    std::vector<Eigen::MatrixXcd> getBerryConnection(Point<T> & point);
protected:
    Particle particle;
    virtual std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
    	diagonalizeFromCoords(Eigen::Vector3d & k);

    // list of lattice vectors, used for the Fourier transform from real
    // to reciprocal space
    Eigen::Matrix<double,3,Eigen::Dynamic> bravaisVectors;
    // cound the vector degeneracies, to use symmetries
    Eigen::VectorXd vectorsDegeneracies;
    Eigen::Matrix3d directUnitCell;
    // hamiltonian matrix in real space <0m|H|nR>
    Eigen::Tensor<std::complex<double>,3> h0R;
    // position matrix elements <0m|r|nR>
    Eigen::Tensor<std::complex<double>,4> rMatrix;

    long numBands;
    long numVectors;
};

template<typename T>
FullBandStructure<T> ElectronH0Wannier::populate(T & fullPoints,
		bool & withVelocities, bool & withEigenvectors) {

	FullBandStructure<T> fullBandStructure(numBands, particle,
			withVelocities, withEigenvectors, fullPoints);

	for ( long ik=0; ik<fullBandStructure.getNumPoints(); ik++ ) {
		Point<T> point = fullBandStructure.getPoint(ik);
		auto [ens, eigvecs] = diagonalize(point);
		fullBandStructure.setEnergies(point, ens);
		if ( withVelocities ) {
			auto vels = diagonalizeVelocity(point);
			fullBandStructure.setVelocities(point, vels);
		}
		//TODO: I must fix the different shape of eigenvectors w.r.t. phonons
//		if ( withEigenvectors ) {
//			fullBandStructure.setEigenvectors(point, eigvecs);
//		}
	}
	return fullBandStructure;
}

template<typename T>
std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
		ElectronH0Wannier::diagonalize(Point<T> & point) {
	Eigen::Vector3d k = point.getCoords(Points::cartesianCoords);

	auto [energies,eigenvectors] = diagonalizeFromCoords(k);

	// note: the eigenvector matrix is the unitary transformation matrix U
	// from the Bloch to the Wannier gauge.

	return {energies, eigenvectors};
}

template<typename T>
Eigen::Tensor<std::complex<double>,3> ElectronH0Wannier::diagonalizeVelocity(
		Point<T> & point) {
	Eigen::Vector3d coords = point.getCoords(Points::cartesianCoords);
	double delta = 1.0e-8;
	double threshold = 0.000001 / energyRyToEv; // = 1 micro-eV
	auto velocity = HarmonicHamiltonian::internalDiagonalizeVelocity(coords,
			delta, threshold);
	return velocity;
}

template<typename T>
std::vector<Eigen::MatrixXcd> ElectronH0Wannier::getBerryConnection(
		Point<T> & point) {
	Eigen::Vector3d k = point.getCoords(Points::cartesianCoords);

	// first we diagonalize the hamiltonian
	auto [ens, eigvecs] = diagonalize(point);

	// note: the eigenvector matrix is the unitary transformation matrix U
	// from the Bloch to the Wannier gauge.

	std::vector<Eigen::MatrixXcd> bc;

	for ( long i=0; i<3; i++ ) {

		// now construct the berryConnection in reciprocal space and Wannier gauge
		Eigen::MatrixXcd berryConnectionW(numBands,numBands);
		berryConnectionW.setZero();

		for ( long iR=0; iR<bravaisVectors.cols(); iR++ ) {
			Eigen::Vector3d R = bravaisVectors.col(iR);
			double phase = k.dot(R);
			std::complex<double> phaseFactor = {cos(phase),sin(phase)};
			for ( long m=0; m<numBands; m++ ) {
				for ( long n=0; n<numBands; n++ ) {
					berryConnectionW(m,n) +=
							phaseFactor * rMatrix(i,iR,m,n)
							/ vectorsDegeneracies(iR);
				}
			}
		}

		Eigen::MatrixXcd berryConnection(numBands,numBands);
		berryConnection = eigvecs.adjoint() * berryConnectionW * eigvecs;
		bc.push_back(berryConnection);
	}
	return bc;
}

#endif
