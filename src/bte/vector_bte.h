#ifndef VECTORBTE_H
#define VECTORBTE_H

#include "eigen.h"
#include "context.h"
#include "active_bandstructure.h"

/** Class used to store the "vector" of out-of-equilibrium populations.
 * The vector indices are over the Bloch state. Additionally, there is one of
 * these vector for each pair of (temperatures,chemicalPotentials).
 * Can be used to store both scalar quantities (linewidths) or vectors
 * (like velocities).
 */
class VectorBTE {
public:
	/** Constructor method, initializes raw buffer data and saves helper
	 * variables.
	 * @param statisticsSweep: saves the info on how many temperatures/chemical
	 * potentials we are evaluating.
	 * @param bandStructure: saves the underlying bandstructure for which we
	 * compute the out-of-equilibrium populations.
	 * @param dimensionality: determines the size of the vector on cartesian
	 * indices. 1 for scalar quantities like linewidths Gamma(BlochIndeces), 3
	 * for vector quantities like phonon populations f(blochIndeces,cartesian).
	 */
	VectorBTE(StatisticsSweep & statisticsSweep_,
			FullBandStructure & bandStructure_,
			const long & dimensionality_=3);
	// copy constructor
	VectorBTE(const VectorBTE & that);
	// copy assignment
	VectorBTE & operator = (const VectorBTE & that);
	// scalar product
	Eigen::VectorXd dot(const VectorBTE & that);
	// product operator overload
	VectorBTE operator * (VectorBTE & that);
	// product operator overload
	VectorBTE operator * (const double & scalar);
	// product operator overload
	VectorBTE operator * (const Eigen::VectorXd & vector);
	// sum operator overload
	VectorBTE operator + (VectorBTE & that);
	// diff operator overload
	VectorBTE operator - (VectorBTE & that);
	// diff operator overload
	VectorBTE operator - ();
	// division operator overload
	VectorBTE operator / (VectorBTE & that);

	VectorBTE sqrt();
	VectorBTE reciprocal();

	void setConst(const double & constant);

	/** Convert an out-of-equilibrium population from the canonical form f to
	 * the absolute value n, such that n = bose(bose+1)f or n=fermi(1-fermi)f.
	 */
	void canonical2Population();

	/** Convert an out-of-equilibrium population from the absolute value n to
	 * the canonical value n, such that n = bose(bose+1)f or n=fermi(1-fermi)f.
	 */
	void population2Canonical();

//protected:

	// we store auxiliary objects and parameters
	StatisticsSweep & statisticsSweep;
	FullBandStructure & bandStructure;
	long numCalcs;
	long numStates;
	long numChemPots;
	long numTemps;
	long dimensionality;

	/** glob2Loc and loc2Glob compress/decompress the indices on temperature,
	 * chemical potential, and cartesian direction into/from a single index.
	 * TODO: these indeces, and how they are used elsewhere, is rather messy
	 * That's because we have to work both with quantities such as linewidths,
	 * which are a scalar over the Bloch states, and phonon populations, which
	 * are cartesian vectors over the Bloch states.
	 * I should probably create two different classes for these.
	 */
	long glob2Loc(const ChemPotIndex & imu, const TempIndex & it,
			const DimIndex & idim);
	std::tuple<ChemPotIndex,TempIndex,DimIndex> loc2Glob(const long & i);

	/** raw buffer containing the values of the vector
	 *  The matrix has size (numCalcs, numStates), where numCalcs is the number
	 *  of pairs of temperature and chemical potentials, and numStates is the
	 *  number of Bloch states used in the Boltzmann equation.
	 */
	Eigen::MatrixXd data;

	friend class ScatteringMatrix; // this is also to remember that
	// if we change the index order of this class, we should check the
	// ScatteringMatrix implementations: they are high efficiency methods
	// so they need low-level access to the raw buffer

	/** base class to implement +, -, / and * operations.
	 * It's split separately so that subclasses can create the correct output
	 * class, and also because operations are rather similar.
	 */
	VectorBTE baseOperator(VectorBTE & that, const int & operatorType);
	const int operatorSums = 0;
	const int operatorDivs = 1;
	const int operatorProd = 2;
	const int operatorDiff = 3;

	/** List of Bloch states to be excluded from the calculation (i.e. for
	 * which vectorBTE values are 0), for example, the acoustic modes at the
	 * gamma point, whose zero frequencies may cause problems.
	 */
	std::vector<long> excludeIndeces;
};

#endif
