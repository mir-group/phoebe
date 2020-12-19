#ifndef BASE_VECTOR_BTE_H
#define BASE_VECTOR_BTE_H

#include "Matrix.h"
#include "PMatrix.h"
#include "active_bandstructure.h"
#include "context.h"
#include "eigen.h"

/** Class used to store a "matrix" of data.
 * The class member "data" has size data(numCalcs,numStates), where numCalcs
 * is the number of temperatures/chemical potentials, and numStates is
 * specified in input. Used for the EPA electron transport calculation.
 * It is subclassed to a VectorBTE case when numStates is aligned with the
 * Bloch states of the band structure class.
 */
class BaseVectorBTE {
 public:
  /** Constructor method, initializes raw buffer data and saves helper
   * variables.
   * @param statisticsSweep: saves the info on how many temperatures/chemical
   * potentials we are evaluating.
   * @param numStates: saves the number of states on which we compute the
   * vector.
   * @param dimensionality: determines the size of the vector on cartesian
   * indices. 1 for scalar quantities like line-widths Gamma(BlochIndices), 3
   * for vector quantities like phonon populations f(blochIndices,cartesian).
   */
  BaseVectorBTE(StatisticsSweep &statisticsSweep_, const long &numStates_,
                const long &dimensionality_ = 3);

  /** Copy constructor
   */
  BaseVectorBTE(const BaseVectorBTE &that);

  /** Copy assignment operator
   */
  BaseVectorBTE &operator=(const BaseVectorBTE &that);

  /** Computes the scalar product between two BaseVectorBTE objects.
   * The scalar product of x and y, is defined such as
   * z(iCalc) = sum_i x(iCalc,i) y(iCalc,i), where i is an index over Bloch
   * states, and iCalc is an index over temperatures and chemical potentials.
   * @param that: the second vector used in the scalar product
   */
  virtual Eigen::MatrixXd dot(const BaseVectorBTE &that);

  /** element wise product between two BaseVectorBTE objects x and y.
   * If the dimensionality of the two objects is the same, we compute
   * element-wise result = x*y.
   * If y has dimensionality 1, we compute x(every dim)*y(0), and the result
   * has the dimensionality of x.
   * @param that: the second BaseVectorBTE object y, such that result = *this*y
   */
  BaseVectorBTE operator*(BaseVectorBTE &that);

  /** Computes the product of a BaseVectorBTE with a scalar, i.e. all elements
   * of vectorBTE x -> x * scalar.
   * @param scalar: a double with the constant factor to be used in the
   * element-wise multiplication.
   */
  BaseVectorBTE operator*(const double &scalar);

  /** Computes the product of a BaseVectorBTE with a vector. The vector has
   * size equal to the number of calculations (i.e. number of temperatures
   * times the number of chemical potentials) used in the run. Given a
   * calculation index iCalc, the result is an element-wise x(it)*vector(it).
   * @param vector: a double vector to be used in the product, of size
   * equal to numCalcs.
   */
  BaseVectorBTE operator*(const Eigen::VectorXd &vector);

  /** Computes the product of a BaseVectorBTE with a parallel matrix. Only works
   * if the number of temperatures/chemical potentials (numCalcs) is equal
   * to one. At fixed calculation index iCalc, the result is an matrix-vector
   * multiplication x(it,i)*pMatrix(i,j).
   * @param pMatrix: a parallel distributed double matrix to be used in the
   * product, of size equal to numStates x numStates.
   */
  BaseVectorBTE operator*(ParallelMatrix<double> &matrix);

  /** element wise sum between two BaseVectorBTE objects x and y.
   * If the dimensionality of the two objects is the same, we compute
   * element-wise result = x+y.
   * If y has dimensionality 1, we compute x(every dim)+y(0), and the result
   * has the dimensionality of x.
   * @param that: the second BaseVectorBTE object y, such that result = *this+y
   */
  BaseVectorBTE operator+(BaseVectorBTE &that);

  /** element wise difference between two BaseVectorBTE objects x and y.
   * If the dimensionality of the two objects is the same, we compute
   * element-wise result = x-y.
   * If y has dimensionality 1, we compute x(every dim)-y(0), and the result
   * has the dimensionality of x.
   * @param that: the second BaseVectorBTE object y, such that result = *this-y
   */
  BaseVectorBTE operator-(BaseVectorBTE &that);

  /** Invert the sign of the BaseVectorBTE content i.e. x -> -x
   */
  BaseVectorBTE operator-();

  /** element wise division between two BaseVectorBTE objects x and y.
   * If the dimensionality of the two objects is the same, we compute
   * element-wise result = x/y.
   * If y has dimensionality 1, we compute x(every dim)/y(0), and the result
   * has the dimensionality of x.
   * @param that: the second BaseVectorBTE object y, such that result = *this/y
   */
  BaseVectorBTE operator/(BaseVectorBTE &that);

  /** Replace the content of BaseVectorBTE with its square root
   * (element-wise x -> sqrt(x) ).
   */
  BaseVectorBTE sqrt();
  /** Replace the content of BaseVectorBTE with its reciprocal
   * (element-wise x -> 1/x).
   */
  BaseVectorBTE reciprocal();

  /** Set the whole content (raw buffer) of BaseVectorBTE to a scalar value.
   * @param constant: the value to be used in the set.
   */
  void setConst(const double &constant);

  /** Get and set operator
   */
  double &operator()(const int iCalc, const int iDim, const int iState);

  /** Const get and set operator
   */
  const double &operator()(const int iCalc, const int iDim,
                           const int iState) const;

  /** raw buffer containing the values of the vector
   *  The matrix has size (numCalcs, numStates), where numCalcs is the number
   *  of pairs of temperature and chemical potentials, and numStates is the
   *  number of Bloch states used in the Boltzmann equation.
   */
  Eigen::MatrixXd data;

  // we store auxiliary objects and parameters
  StatisticsSweep &statisticsSweep;
  long numCalcs;
  long numStates;
  long numChemPots;
  long numTemps;
  long dimensionality;

  /** glob2Loc and loc2Glob compress/decompress the indices on temperature,
   * chemical potential, and cartesian direction into/from a single index.
   * TODO: these indices, and how they are used elsewhere, is rather messy
   * That's because we have to work both with quantities such as line-widths,
   * which are a scalar over the Bloch states, and phonon populations, which
   * are cartesian vectors over the Bloch states.
   * I should probably create two different classes for these.
   */
  long glob2Loc(const ChemPotIndex &imu, const TempIndex &it,
                const CartIndex &idim);
  std::tuple<ChemPotIndex, TempIndex, CartIndex> loc2Glob(const long &i);

  /** List of Bloch states to be excluded from the calculation (i.e. for
   * which vectorBTE values are 0), for example, the acoustic modes at the
   * gamma point, whose zero frequencies may cause problems.
   */
  std::vector<long> excludeIndices;

 protected:
  /** base class to implement +, -, / and * operations.
   * It's split separately so that subclasses can create the correct output
   * class, and also because operations are rather similar.
   */
  BaseVectorBTE baseOperator(BaseVectorBTE &that, const int &operatorType);
  const int operatorSums = 0;
  const int operatorDivs = 1;
  const int operatorProd = 2;
  const int operatorDiff = 3;
};

#endif
