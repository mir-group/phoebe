#ifndef VECTORBTE_H
#define VECTORBTE_H

#include "eigen.h"
#include "context.h"
#include "active_bandstructure.h"
#include "Matrix.h"
#include "PMatrix.h"
#include "basevector_bte.h"

/** Class used to store the "vector" of out-of-equilibrium populations.
 * The vector indices are over the Bloch state. Additionally, there is one of
 * these vector for each pair of (temperatures,chemicalPotentials).
 * Can be used to store both scalar quantities (linewidths) or vectors
 * (like velocities).
 */
class VectorBTE : public BaseVectorBTE {
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
    VectorBTE(StatisticsSweep &statisticsSweep_,
            BaseBandStructure &bandStructure_, const long &dimensionality_ = 3);

    /** Copy constructor
     */
    VectorBTE(const VectorBTE &that);

    /** Copy assignment operator
     */
    VectorBTE& operator =(const VectorBTE &that);

    /** Computes the scalar product between two VectorBTE objects.
     * The scalar product of x and y, is defined such as
     * z(iCalc) = sum_i x(iCalc,i) y(iCalc,i), where i is an index over Bloch
     * states, and iCalc is an index over temperatures and chemical potentials.
     * @param that: the second vector used in the scalar product
     */
    Eigen::VectorXd dot(const VectorBTE &that);

    /** element wise product between two VectorBTE objects x and y.
     * If the dimensionality of the two objects is the same, we compute
     * element-wise result = x*y.
     * If y has dimensionality 1, we compute x(every dim)*y(0), and the result
     * has the dimensionality of x.
     * @param that: the second VectorBTE object y, such that result = *this*y
     */
    VectorBTE operator *(VectorBTE &that);

    /** Computes the product of a VectorBTE with a scalar, i.e. all elements
     * of vectorBTE x -> x * scalar.
     * @param scalar: a double with the constant factor to be used in the
     * element-wise multiplication.
     */
    VectorBTE operator *(const double &scalar);

    /** Computes the product of a VectorBTE with a vector. The vector has
     * size equal to the number of calculations (i.e. number of temperatures
     * times the number of chemical potentials) used in the run. Given a
     * calculation index iCalc, the result is an element-wise x(it)*vector(it).
     * @param vector: a double vector to be used in the product, of size
     * equal to numCalcs.
     */
    VectorBTE operator *(const Eigen::VectorXd &vector);

    /** Computes the product of a VectorBTE with a parallel matrix. Only works
     * if the number of temperatures/chemical potentials (numCalcs) is equal
     * to one. At fixed calculation index iCalc, the result is an matrix-vector
     * multiplication x(it,i)*pMatrix(i,j).
     * @param pMatrix: a parallel distributed double matrix to be used in the
     * product, of size equal to numStates x numStates.
     */
    VectorBTE operator *(ParallelMatrix<double> &matrix);

    /** element wise sum between two VectorBTE objects x and y.
     * If the dimensionality of the two objects is the same, we compute
     * element-wise result = x+y.
     * If y has dimensionality 1, we compute x(every dim)+y(0), and the result
     * has the dimensionality of x.
     * @param that: the second VectorBTE object y, such that result = *this+y
     */
    VectorBTE operator +(VectorBTE &that);

    /** element wise difference between two VectorBTE objects x and y.
     * If the dimensionality of the two objects is the same, we compute
     * element-wise result = x-y.
     * If y has dimensionality 1, we compute x(every dim)-y(0), and the result
     * has the dimensionality of x.
     * @param that: the second VectorBTE object y, such that result = *this-y
     */
    VectorBTE operator -(VectorBTE &that);

    /** Invert the sign of the VectorBTE content i.e. x -> -x
     */
    VectorBTE operator -();

    /** element wise division between two VectorBTE objects x and y.
     * If the dimensionality of the two objects is the same, we compute
     * element-wise result = x/y.
     * If y has dimensionality 1, we compute x(every dim)/y(0), and the result
     * has the dimensionality of x.
     * @param that: the second VectorBTE object y, such that result = *this/y
     */
    VectorBTE operator /(VectorBTE &that);

    /** Replace the content of VectorBTE with its square root
     * (element-wise x -> sqrt(x) ).
     */
    VectorBTE sqrt();
    /** Replace the content of VectorBTE with its reciprocal
     * (element-wise x -> 1/x).
     */
    VectorBTE reciprocal();

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
    BaseBandStructure &bandStructure;

    friend class ScatteringMatrix; // this is also to remember that
    // if we change the index order of this class, we should check the
    // ScatteringMatrix implementations: they are high efficiency methods
    // so they need low-level access to the raw buffer

    /** base class to implement +, -, / and * operations.
     * It's split separately so that subclasses can create the correct output
     * class, and also because operations are rather similar.
     */
    VectorBTE baseOperator(VectorBTE &that, const int &operatorType);
    const int operatorSums = 0;
    const int operatorDivs = 1;
    const int operatorProd = 2;
    const int operatorDiff = 3;
};

#endif
