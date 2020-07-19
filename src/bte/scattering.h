#ifndef SCATTERING_H
#define SCATTERING_H

#include "context.h"
#include "vector_bte.h"
#include "delta_function.h"
#include "Matrix.h"

/** Base class of the scattering matrix.
 * Note: this is an abstract class, which can only work if builder() is defined
 */
class ScatteringMatrix {
public:
    /** Scattering matrix constructor
     * @param context: object containing user input configuration.
     * @param statisticsSweep: object controlling the loops over temperature
     * and chemical potential.
     * @param innerBandStructure: bandStructure object. This is the mesh used
     * to integrate the anharmonic properties for each state of outerBandStruc.
     * For transport calculation, this object is typically equal to
     * outerBandStructure. Might differ when outerBS is on a path of points.
     * @param outerBandStructure: bandStructure object. The anharmonic
     * properties are computed on the grid of points specified by this object,
     * e.g. this could be a path of points or a uniform mesh of points
     */
    ScatteringMatrix(Context &context_, StatisticsSweep &statisticsSweep_,
            BaseBandStructure &innerBandStructure_,
            BaseBandStructure &outerBandStructure_);

    /** Copy constructor
     */
    ScatteringMatrix(const ScatteringMatrix &that);

    /** Copy assignment operator
     */
    ScatteringMatrix& operator=(const ScatteringMatrix &that);

    /** Destructor
     */
    ~ScatteringMatrix();

    /** This method needs to be called right after the constructor.
     * It will setup variables of the scattering matrix, and compute at least
     * the linewidths.
     * If the user input ask to store the matrix in memory, the scattering
     * matrix gets allocated and built here, through a call to a virtual member
     * builder(), which is defined in subclasses.
     */
    void setup();

    /** Returns the diagonal matrix elements.
     * @return diagonal: a VectorBTE containing the linewidths.
     */
    VectorBTE diagonal();

    /** Computes the product A*f - diag(A)*f
     * where A is the scattering matrix, f is the vector of quasiparticle
     * populations, and diag(A) is the diagonal of the matrix.
     */
    VectorBTE offDiagonalDot(VectorBTE &inPopulation);

    /** Computes the product A*f
     * where A is the scattering matrix and f is the vector of quasiparticle
     * populations.
     */
    VectorBTE dot(VectorBTE &inPopulation);

    /** Computes the product A*B, where A is the scattering matrix, and
     * B is an Eigen::MatrixXd. This can be used to compute products of the
     * scattering matrix with other vectors.
     */
    ParallelMatrix<double> dot(const ParallelMatrix<double> &otherMatrix);

    /** Call to obtain the single-particle relaxation times of the system.s
     * @return tau: a VectorBTE object storing the relaxation times
     */
    VectorBTE getSingleModeTimes();

    /** Converts the scttering matrix from the form A to the symmetrized Omega.
     * A acts on the canonical phonon population f, while Omega acts on the
     * symmetrized phonon population \tilde{n}.
     * Note: for phonons, n = bose(bose+1)f , with bose being the
     * bose--einstein distribution, while n = sqrt(bose(bose+1)) tilde(n).
     * Only works if the matrix is kept in memory and after setup() has been
     * called.
     */
    void a2Omega();

    /** The inverse of a2Omega, converts the matrix Omega to A
     */
    void omega2A();

    /** Diagonalize the scattering matrix
     * @return eigenvalues: a VectorBTE object with the eigenvalues
     * @return eigenvectors: a Eigen::MatrixXd with the eigenvectors
     * Eigenvectors are aligned on rows: eigenvectors(qpState,eigenIndex)
     */
    std::tuple<VectorBTE, ParallelMatrix<double>> diagonalize();

protected:
    Context &context;
    StatisticsSweep &statisticsSweep;

    // Smearing is a pointer created in constructor with a smearing factory
    // Used in the construction of the scattering matrix to approximate the
    // Dirac-delta function in transition rates.
    DeltaFunction *smearing;

    BaseBandStructure &innerBandStructure;
    BaseBandStructure &outerBandStructure;

    // constant relaxation time approximation -> the matrix is just a scalar
    // and there are simplified evaluations taking place
    bool constantRTA = false;
    bool highMemory = true; // whether the matrix is kept in memory
    bool isMatrixOmega = false; // whether the matrix is Omega or A
    // A acts on the canonical population, omega on the population
    // A acts on f, omega on n, with n = bose(bose+1)f e.g. for phonons

    // we save the diagonal matrix element in a dedicated vector
    VectorBTE internalDiagonal;
    // the scattering matrix, initialized if highMemory==true
    ParallelMatrix<double> theMatrix;

    long numStates; // number of Bloch states (i.e. the size of theMatrix)
    long numPoints; // number of wavevectors
    long numCalcs; // number of "Calculations", i.e. number of temperatures and
    // chemical potentials on which we compute scattering.
    int dimensionality_;

    // this is to exclude problematic Bloch states, e.g. acoustic phonon modes
    // at gamma, which have zero frequencies and thus have several non-analytic
    // behaviors (e.g. Bose--Einstein population=infinity). We set to zero
    // terms related to these states.
    std::vector<long> excludeIndeces;

    /** Method that actually computes the scattering matrix.
     * Pure virtual function: needs an implementation in every subclass.
     * Builder has three behaviors:
     * 1) if matrix.size()==0 and linewidth is passed, builder computes the
     * quasiparticle linewidths.
     * 2) if matrix.size > 0 and linewidth is passed, builder computes the
     * quasiparticle linewidths and the scattering matrix. Memory intensive!
     * 3) if matrix.size()==0, linewidth is not passed, but we pass in+out
     * populations, we compute outPopulation = scattMatrix * inPopulation.
     * This doesn't require to store the matrix in memory.
     */
    virtual void builder(ParallelMatrix<double> &matrix, VectorBTE *linewidth,
            VectorBTE *inPopulation, VectorBTE *outPopulation) = 0;

    /** Returns a vector of pairs of wavevector indices to iterate over during
     * the construction of the scattering matrix.
     * @param switchCase: if 0, returns the pairs of wavevectors to loop for
     * the case where the scattering matrix is built in memory.
     * If != 0, returns the pairs of wavevectors to loop for the case where
     * only the action of the scattering matrix is computed.
     * @return vector<tuple<iq1,iq2>>: a tuple of wavevector indices to loop
     * over in the construction of the scattering matrix.
     */
    std::vector<std::tuple<std::vector<long>,long>> getIteratorWavevectorPairs(
            const int & switchCase);
};

#endif
