#ifndef SCATTERING_H
#define SCATTERING_H

#include "context.h"
#include "vector_bte.h"
#include "delta_function.h"

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
	ScatteringMatrix(Context & context_, StatisticsSweep & statisticsSweep_,
			FullBandStructure<FullPoints> & innerBandStructure_,
			FullBandStructure<FullPoints> & outerBandStructure_);

	/** Copy constructor
	 *
	 */
	ScatteringMatrix(const ScatteringMatrix & that);

	/** Assignment operator
	 *
	 */
	ScatteringMatrix & operator=(const ScatteringMatrix & that);

	/** Destructor
	 *
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
	VectorBTE offDiagonalDot(VectorBTE & inPopulation);

	/** Computes the product A*f
	 * where A is the scattering matrix and f is the vector of quasiparticle
	 * populations.
	 */
	VectorBTE dot(VectorBTE & inPopulation);

	/** Computes the product A*B, where A is the scattering matrix, and
	 * B is an Eigen::MatrixXd. This can be used to compute products of the
	 * scattering matrix with other vectors.
	 */
	Eigen::MatrixXd dot(const Eigen::MatrixXd & otherMatrix);

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
	std::tuple<VectorBTE,Eigen::MatrixXd> diagonalize();
protected:
	Context & context;
	StatisticsSweep & statisticsSweep;
	DeltaFunction * smearing;

	FullBandStructure<FullPoints> & innerBandStructure;
	FullBandStructure<FullPoints> & outerBandStructure;

	 // constant relaxation time approximation -> the matrix is just a scalar
	bool constantRTA = false;
	bool highMemory = true;
	bool hasCGScaling = false;
	bool isMatrixOmega = false;

	VectorBTE internalDiagonal;
	Eigen::MatrixXd theMatrix;
	long numStates;
	long numPoints;
	long numCalcs;

	std::vector<long> excludeIndeces;

	// pure virtual function
	// needs an implementation in every subclass
	virtual void builder(Eigen::MatrixXd & matrix, VectorBTE * linewidth,
			VectorBTE * inPopulation, VectorBTE * outPopulation) = 0;
};

#endif
