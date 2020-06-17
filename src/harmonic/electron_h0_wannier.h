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
class ElectronH0Wannier: public HarmonicHamiltonian {
public:
    ElectronH0Wannier(const Eigen::Matrix3d &directUnitCell_,
            const Eigen::Matrix<double, 3, Eigen::Dynamic> &bravaisVectors_,
            const Eigen::VectorXd &vectorsDegeneracies_,
            const Eigen::Tensor<std::complex<double>, 3> &h0R_,
            const Eigen::Tensor<std::complex<double>, 4> &rMatrix_);

    /** Copy constructor
     */
    ElectronH0Wannier(const ElectronH0Wannier &that);

    /** Copy assignment
     */
    ElectronH0Wannier& operator =(const ElectronH0Wannier &that);

    /** Method to return that the underlying is that of an electronic Fermion.
     */
    Particle getParticle();

    /** get the total number of bands.
     * This is a constant for all wavevectors.
     */
    long getNumBands();

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
    std::tuple<Eigen::VectorXd, Eigen::MatrixXcd> diagonalize(Point &point);

    virtual std::tuple<Eigen::VectorXd, Eigen::MatrixXcd>
    diagonalizeFromCoords(Eigen::Vector3d &k);

    /** get the electron velocities (in atomic units) at a single k-point.
     * @param k: a Point object with the wavevector coordinates.
     * @return velocity(numBands,numBands,3): values of the velocity operator
     * for this state, in atomic units.
     */
    Eigen::Tensor<std::complex<double>, 3> diagonalizeVelocity(Point &point);

    /** This method constructs an electron bandstructure.
     * @param points: the object with the list/mesh of wavevectors
     * @param withVelocities: if true, compute the electron velocity operator.
     * @param withEigenvectors: if true, stores the Wannier eigenvectors.
     * @return FullBandStructure: the bandstructure object containing the
     * complete electronic band structure.
     */
    FullBandStructure populate(Points &fullPoints, bool &withVelocities,
            bool &withEigenvectors);

    /** compute the Berry connection <u_mk| nabla_k |u_nk> at arb. wavevectors.
     * @param point: the Point coordinates of the wavevector.
     * @return Berry connection: a generalized Berry connection in the form of
     * a matrix <u_mk| nabla_k |u_nk> for a fixed wavevector. The Berry
     * connection is actually just the diagonal matrix elements.
     */
    std::vector<Eigen::MatrixXcd> getBerryConnection(Point &point);
protected:
    Particle particle;

    // list of lattice vectors, used for the Fourier transform from real
    // to reciprocal space
    Eigen::Matrix<double, 3, Eigen::Dynamic> bravaisVectors;
    // cound the vector degeneracies, to use symmetries
    Eigen::VectorXd vectorsDegeneracies;
    Eigen::Matrix3d directUnitCell;
    // hamiltonian matrix in real space <0m|H|nR>
    Eigen::Tensor<std::complex<double>, 3> h0R;
    // position matrix elements <0m|r|nR>
    Eigen::Tensor<std::complex<double>, 4> rMatrix;

    long numBands;
    long numVectors;
};

#endif
