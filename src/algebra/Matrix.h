#include "PMatrix.h"
// note: this goes outside MPI_AVAIL
// in the rest of the code we can just load this header, it will be subsitituted
// by the Matrix class defined in PMatrix.h, which relies on scalapack

#ifndef MPI_AVAIL

#ifndef MATRIX_H
#define MATRIX_H

// include statements
#include <type_traits>
#include <complex>
#include <assert.h>
#include <vector>
#include <iostream>
#include <cmath>
#include "Blas.h" 

// Trick to allow type promotion below
template <typename T>
struct identity_t { typedef T type; };

// Allow promotion of complex types (for example, complex<double> * int)
#define COMPLEX_OPS(OP)                                                         \
	template <typename _Tp> std::complex<_Tp>                               \
	operator OP(std::complex<_Tp> lhs, const typename identity_t<_Tp>::type & rhs) \
	{                                                                       \
		return lhs OP rhs;                                              \
	}                                                                       \
	template <typename _Tp> std::complex<_Tp>                               \
	operator OP(const typename identity_t<_Tp>::type & lhs, const std::complex<_Tp> & rhs) \
	{                                                                       \
		return lhs OP rhs;                                              \
	}
COMPLEX_OPS(+)
COMPLEX_OPS(-)
COMPLEX_OPS(*)
COMPLEX_OPS(/)
#undef COMPLEX_OPS

//! Matrix parent class, which can be used to define matrix classes of different types
/** \file matrix.h
 * \brief General templated matrix class, with explicit specialization for double and complex<double> types. */
template <typename T>
class Matrix
{
	/// Class variables
	int nRows;
	int nCols;
	
	public:

		T* mat; // pointer to the internal array structure. NOTE: can we find a way to make this private?

		/// Matrix class constructors -----------------------------------
		Matrix(const int rows, const int cols); // Construct using row, col numbers
		Matrix(); // default constructor
		// NOTE: currently relying on the default destructor instead.
		// if I add this explicit delete, I get a double free error when the program exits 
		//~Matrix(){delete[] mat;}  
		/// Copy constructor
		Matrix(const Matrix<T>& m1); 
		
		/// Matrix class methods -----------------------------------
		int rows() const;
		int cols() const;
		int size() const;
                /// A reshape function. Note, this is not the same as pythons -- python wraps around the rows, 
                /// this function instead wraps down columns. This is because the matrix is set up in col major order. 
		void reshape(const int rows, const int cols);
		/// Print the matrix
		void print() const;
		void copy(const Matrix<T>& toCopy);         

		/// Index from a 1D array to a position in a 2D array (matrix)
		int index(int row,int col) const;

		// Get and set operators
		/// Get an element of the matrix
		T& operator()(const int row, const int col) const;
		/// Set an element of the matrix
		void operator()(const int row, const int col, const T value);
 
		/// Get a block of the matrix
		Matrix<T> operator()(const int rowStart, const int rowStop, const int colStart, const int colStop) const;
		/// Take a whole row or column slice
		Matrix<T> getRow(const int row) const;
		Matrix<T> getCol(const int col) const;
	
		/// Set a block of the matrix
		void operator()(const int rowStart, const int rowStop, const int colStart, const int colStop, const Matrix<T>& value);
		/// Set a whole row or column slice
		void setRow(const int row, const Matrix<T>& value);
		void setCol(const int col, const Matrix<T>& value);
	
		/// Addition and subtraction member functions
		// Generic matrix addition.  
		template<typename U>
		Matrix<T>& operator+=(const Matrix<U>& m1) {
			assert( (m1.rows() == nRows) && (m1.cols() == nCols) );
			static_assert( std::is_same<T, typename std::common_type<U,T>::type >::value, "You cannot += against a type which would return a type other than the lhs value." );
			for(int s = 0; s<getSize(); s++) mat[s] += m1.mat[s];
			return *this;
		} 
		// Generic matrix subtraction. 
		template<typename U>
		Matrix<T>& operator-=(const Matrix<U>& m1) {
			assert( (m1.rows() == nRows) && (m1.cols() == nCols) );
			static_assert( std::is_same<T, typename std::common_type<U,T>::type >::value, "You cannot -= against a type which would return a type other than the lhs value." );
			for(int s = 0; s<getSize(); s++) mat[s] -= m1.mat[s];
			return *this;
		} 
		// Generic matrix multiplication. 
		template<typename U> Matrix<T>& operator*=(const Matrix<U>& m1) {
			assert(nRows == m1.cols());
			static_assert( std::is_same<T, typename std::common_type<U,T>::type >::value, "You cannot *= against a type which would return a type other than the lhs value." );

			Matrix<T> ret(nRows,m1.cols());
			// Loop over the rows of this matrix
			for(int i =0; i<nRows; i++) {
				// loop over the cols of this matrix and rows of m1
				for(int j =0; j<nCols; j++) {
						// loop over the cols of m1
						for(int k =0; k<m1.cols(); k++) ret(i,k) += (*this)(i,j) * m1(j,k);
				}
			} 
			this->mat = ret.mat; 
			return *this;
		}

		/// Declaration of nonmember friend functions for operators like +, -, *. auto type here specifies the use of a "trailing return type"
		// Generic matrix addition. 
		template<typename U, typename V>
		friend auto operator+(const Matrix<U>& m1, const Matrix<V>& m2) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> ;
		// Generic matrix subtraction. 
		template<typename U, typename V>
		friend auto operator-(const Matrix<U>& m1, const Matrix<V>& m2) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type>; 
		/// Multiplication-based operations
		/// General matrix multiplication
		template<typename U, typename V>
		friend auto operator*(const Matrix<U>& m1, const Matrix<V>& m2) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type>; 

		// Matrix times a scalar. -- TODO: write specific versions with  z/dscal
		// Generic matrix * scalar product
		template<typename U, typename V>
		friend auto operator*( const Matrix<U>& m1,  const V scalar) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type>; 

		// TODO: replace with z or dscal
		// Generic matrix times a scalar, have to define again with argument order flipped.
		template<typename U, typename V>
		friend auto operator*(const V scalar, const Matrix<U>& m1) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type>; 
		// Generic matrix / a scalar -- TODO: replace with z or dscal
		template<typename U, typename V>
		friend auto operator/( const Matrix<U>& m1,  const V scalar) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type>; 
		
		/// Unary negation
		Matrix<T> operator-() const;
		/// Equivalence operators
		bool operator==(const Matrix<T>& m2) const;
		bool operator!=(const Matrix<T>& m2) const;

		/// Non-member functions specfic to complex matrices. 
		template<typename U> friend Matrix<U> real(const Matrix<std::complex<U>>& m);
		template<typename U> friend Matrix<U> imag(const Matrix<std::complex<U>>& m);

		/// Useful matrix functions
		T trace() const;        ///< Returns the trace of a matrix
                T det() const;                        ///< Calculate the determinant of a diagonal matrix
                double norm() const;         ///< Returns the Frobenius norm of a matrix 
		void transposeIP(Matrix<T>& m);     ///< A function to transpose in place
		Matrix<T> transpose();  			///< Diagonalize a matrix
		void conjugateIP(Matrix<T>& m);     ///< A function to transpose in place
		Matrix<T> conjugate();  ///< Returns the complex conjugate of the matrix
		void eye();             ///< Sets this matrix as the identity
	
		/// Matrix algebra functions
		void diagonalize(Matrix<std::complex<double>>& eigvecs, Matrix<std::complex<double>>& eigvals) const; ///< Diagonalize a real matrix
		Matrix<T> dagger(); ///<  Provide  A(dag)

	private:

	
}; // end of class Matrix

//  forward declarations to satisfy compiler when explicitly specializing
template <> void Matrix<std::complex<double>>::conjugateIP(Matrix<std::complex<double>>& m);
template <> void Matrix<double>::diagonalize(Matrix<std::complex<double> >& eigvecs, Matrix<std::complex<double> >& eigvals) const;

// Implementation of non-member matrix functions
// Generic matrix addition
template<typename U, typename V>
auto operator+(const Matrix<U>& m1, const Matrix<V>& m2) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type>{
	assert( (m1.rows() == m2.rows()) || m1.cols() == m2.cols() );
	typedef typename std::common_type<decltype(U{}),decltype(V{})>::type W; 
	Matrix<W> c(m1.rows(), m1.cols());
	for(int s = 0; s<m1.getSize(); s++) c.mat[s] = m1.mat[s] + m2.mat[s];  
	return c;
}
// Generic matrix subtraction. 
template<typename U, typename V>
auto operator-(const Matrix<U>& m1, const Matrix<V>& m2) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> {
	assert( (m1.rows() == m2.rows()) || m1.cols() == m2.cols() );
	typedef typename std::common_type<decltype(U{}),decltype(V{})>::type W; 
	Matrix<W> c(m1.rows(), m1.cols());
	for(int s = 0; s<m1.getSize(); s++) c.mat[s] = m1.mat[s] - m2.mat[s];
	return c;
}
// Generic matrix multiplication.  
template<typename U, typename V>
auto operator*(const Matrix<U>& m1, const Matrix<V>& m2) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> {
	assert(m1.cols() == m2.rows());
	typedef typename std::common_type<decltype(U{}),decltype(V{})>::type W;
	Matrix<W> ret(m1.rows(),m1.cols());
	// Loop over the rows of this matrix
	for(int i =0; i<m1.rows(); i++) {
		// loop over the cols of this matrix and rows of m1
		for(int j =0; j<m1.cols(); j++) {
				// loop over the cols of m1
				for(int k =0; k<m2.cols(); k++) ret(i,k) += m1(i,j) * m2(j,k);
		}
	}
	return ret;
}
// Explict specialization of BLAS matrix-matrix mult for Matrix<complex<double>>
template <>
inline Matrix<std::complex<double>> operator*(const Matrix<std::complex<double>>& m1, const Matrix<std::complex<double>>& m2) {
	assert(m1.cols() == m2.rows());
	Matrix<std::complex<double>> ret(m1.rows(),m2.cols()); // newly sized matrix
	// throw away variables
	char transa = 'n'; char transb = 'n';
	std::complex<double> alpha(1.0,0.0);
	std::complex<double> beta(0.0,0.0);
	zgemm_(transa,transb,m1.rows(),m2.cols(),m1.cols(),alpha,m1.mat,m1.rows(), m2.mat,m2.rows(),beta,ret.mat,m1.rows());
	return ret;
}
// Explicit specializiation of BLAS matrix-matrix mult for Matrix<double>
template <>
inline Matrix<double> operator*(const Matrix<double>& m1, const Matrix<double>& m2) {
	assert(m1.cols() == m2.rows());
	Matrix<double> ret(m1.rows(),m2.cols()); // newly sized matrix
	// throw away variables
	char transa = 'n'; char transb = 'n'; // compute only reigs
	double alpha = 1.0;
	double beta = 0.0;
	dgemm_(transa,transb,m1.rows(),m2.cols(),m1.cols(),alpha,m1.mat,m1.rows(), m2.mat,m2.rows(),beta,ret.mat,m1.rows());
	return ret;
}
// Generic matrix * scalar product
template<typename U, typename V>
auto operator*( const Matrix<U>& m1,  const V scalar) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> {
	typedef typename std::common_type<decltype(U{}),decltype(V{})>::type W;
	Matrix<W> c(m1.rows(), m1.cols());
	for(int s = 0; s<m1.getSize(); s++) c.mat[s] = m1.mat[s] * scalar;
	return c;
}
// TODO: replace with z/dscal spcialization
// Generic matrix times a scalar, have to define again with argument order flipped.
template<typename U, typename V>
auto operator*(const V scalar, const Matrix<U>& m1) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> {
	typedef typename std::common_type<decltype(U{}),decltype(V{})>::type W;
	Matrix<W> c(m1.rows(), m1.cols());
	for(int s = 0; s<m1.getSize(); s++) c.mat[s] = m1.mat[s] * scalar;
	return c;
}
// Generic matrix / a scalar -- TODO: replace with z/dscal spcialization
template<typename U, typename V>
auto operator/( const Matrix<U>& m1,  const V scalar) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> {
	typedef typename std::common_type<decltype(U{}),decltype(V{})>::type W;
	Matrix<W> c(m1.rows(), m1.cols());
	for(int s = 0; s<m1.getSize(); s++) c.mat[s] = m1.mat[s] / scalar;
	return c;
} 
/// Non-member functions specfic to complex matrices. TODO: is there a faster way to do this?
/// Get the real part of a complex matrix
template<typename U> Matrix<U> real(const Matrix<std::complex<U>>& m) {
	Matrix<U> ret(m.rows(),m.cols());
	for(int s = 0; s<m.getSize(); s++) ret.mat[s] = m.mat[s].real();
	return ret;   
}
/// Get the imag part of a complex matrix
template<typename U> Matrix<U> imag(const Matrix<std::complex<U>>& m) {
	Matrix<U> ret(m.rows(),m.cols());
	for(int s = 0; s<m.getSize(); s++) ret.mat[s] = m.mat[s].imag();
	return ret;
}
/// Set the real part of a complex matrix
template<typename U> void setReal(Matrix<std::complex<U>>& m, Matrix<U>& realMat) {
	for(int s = 0; s<m.getSize(); s++) m.mat[s].real(realMat.mat[s]);
}
/// Set the imag part of a complex matrix
template<typename U> void setImag(Matrix<std::complex<U>>& m, Matrix<U>& imagMat) {
	for(int s = 0; s<m.getSize(); s++) m.mat[s].imag(imagMat.mat[s]);
}

/* This is here to fix linker errors which otherwise require the
 entire template class to be defined in this header file */
#include "Matrix.cpp"

#endif // MATRIX_H

#endif // mpi_AVAIL
