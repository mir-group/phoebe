
#ifndef MATRIX_H
#define MATRIX_H

// include statements
#include <Eigen/Dense>

//! Matrix parent class, which can be used to define matrix classes of different types
/** \file matrix.h
 * \brief Header file for a basic complex matrix parent class, as well as several specialized matrix types */
template <typename T>
class Matrix
{
    /// Class variables
    int nRows;
    int nCols;
    T* mat; // pointer to the internal array structure
    
    public:

        /// Matrix class constructors -----------------------------------
        Matrix(const int rows, const int cols); // Construct using row, col numbers
        ~Matrix(){delete[] mat};   // currently relying on the default destructor instead
        /// Copy constructor
        Matrix(const Matrix& m1); 
        
        /// Matrix class methods -----------------------------------
        int numRows() const;
        int numCols() const;
        int getSize() const;
        void reshape(const int rows, const int cols) { nRows = rows; nCols = cols; } 
        /// Print the matrix
        void print() const;
        
        /// Index from a 1D array to a position in a 2D array (matrix)
        int index(int row,int col) const;

        // Get and set operators
        /// Get an element of the matrix
        T operator()(const int row, const int col) const;
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
        // Brute force generic matrix addition.  
        template<typename U>
        Matrix<T>& operator+=(const Matrix<U>& m1) {
            assert( (m1.numRows() == nRows) && (m1.numCols() == nCols) );
            static_assert( std::is_same<T, typename std::common_type<U,T>::type >::value, "You cannot += against a type which would return a type other than the lhs value." );
            for(int s = 0; s<getSize(); s++) mat[s] += m1.mat[s];
            return *this;
        } 
        // Brute force generic matrix subtraction. 
        template<typename U>
        Matrix<T>& operator-=(const Matrix<U>& m1) {
            assert( (m1.numRows() == nRows) && (m1.numCols() == nCols) );
            static_assert( std::is_same<T, typename std::common_type<U,T>::type >::value, "You cannot -= against a type which would return a type other than the lhs value." );
            for(int s = 0; s<getSize(); s++) mat[s] -= m1.mat[s];
            return *this;
        } // Generic brute force matrix multiplication. 
        template<typename U> Matrix<T>& operator*=(const Matrix<U>& m1) {
            assert(nRows == m1.numCols());
            static_assert( std::is_same<T, typename std::common_type<U,T>::type >::value, "You cannot *= against a type which would return a type other than the lhs value." );

            Matrix<T> ret(nRows,m1.numCols()); // newly sized matrix 
            // Loop over the rows of this matrix
            for(int i =0; i<nRows; i++) {
                // loop over the cols of this matrix and rows of m1
                for(int j =0; j<nCols; j++) {
                        // loop over the cols of m1
                        for(int k =0; k<m1.numCols(); k++) ret(i,k) += mat(i,j) * m1(j,k); 
                }
            } 
            this.mat = ret.mat; 
            return *this;
        }

        /// Nonmember friend functions for operators like +, -, *. auto type here specifies the use of a "trailing return type" after the
        /// function delcaration
        // Generic brute force matrix addition. 
        template<typename U, typename V>
        friend auto operator+(const Matrix<U>& m1, const Matrix<V>& m2) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> {
            assert( (m1.numRows() == m2.numRows()) || m1.numCols() == m2.numCols() );
            typedef W = typename std::common_type<decltype(U{}),decltype(V{})>::type; 
            Matrix<W> c(m1.numRows(), m1.numCols());
            for(int s = 0; s<m1.getSize(); s++) c.mat[s] = m1[s] + m2[s];  
            return c;
        }
        // Generic brute force matrix subtraction. 
        template<typename U, typename V>
        friend auto operator-(const Matrix<U>& m1, const Matrix<V>& m2) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> {
            assert( (m1.numRows() == m2.numRows()) || m1.numCols() == m2.numCols() );
            typedef W = typename std::common_type<decltype(U{}),decltype(V{})>::type; 
            Matrix<W> c(m1.numRows(), m1.numCols());
            for(int s = 0; s<m1.getSize(); s++) c.mat[s] = m1[s] - m2[s];
            return c;
        } 
        /// Multiplication-based operations
        /// General matrix multiplication
        template<typename U, typename V>
        friend auto operator*(const Matrix<U>& m1, const Matrix<V>& m2) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type>  {
            assert(m1.numCols() == m2.numRows());
            typedef W = typename std::common_type<decltype(U{}),decltype(V{})>::type;
            Matrix<T> ret(m1.numRows,m1.numCols()); // newly sized matrix 

            // Loop over the rows of this matrix
            for(int i =0; i<m1.numRows(); i++) {
                // loop over the cols of this matrix and rows of m1
                for(int j =0; j<m1.numCols(); j++) {
                        // loop over the cols of m1
                        for(int k =0; k<m2.numCols(); k++) ret(i,k) += m1(i,j) * m2(j,k);
                }
            }
            return ret;
        }
        // Matrix times a scalar. -- TODO: write specific versions with  z/dscal
        // Generic matrix * scalar product
        template<typename U, typename V>
        friend auto operator*( const Matrix<U>& m1,  const V scalar) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> {
            typedef W = typename std::common_type<decltype(U{}),decltype(V{})>::type;
            Matrix<W> c(m1.numRows(), m2.numCols());
            for(int s = 0; s<m1.getSize(); s++) c.mat[s] = m1[s] * scalar;
            return c;
        }
        // TODO: replace with z or dscal
        // Generic matrix times a scalar, have to define again with argument order flipped.
        template<typename U, typename V>
        friend auto operator*(const V scalar, const Matrix<U>& m1) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> {
            typedef W = typename std::common_type<decltype(U{}),decltype(V{})>::type;
            Matrix<W> c(m1.numRows(), m1.numCols());
            for(int s = 0; s<m1.getSize(); s++) c.mat[s] = m1[s] * scalar;
            return c;
        }
        // Generic matrix / a scalar -- TODO: replace with z or dscal
        template<typename U, typename V>
        friend auto operator/( const Matrix<U>& m1,  const V scalar) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> {
            typedef W = typename std::common_type<decltype(U{}),decltype(V{})>::type;
            Matrix<W> c(m1.numRows(), m1.numCols());
            for(int s = 0; s<m1.getSize(); s++) c.mat[s] = m1[s] / scalar;
            return c;
        }
        
        /// Unary negation
        Matrix<T> operator-() const;
        /// Equivalence operators
        bool operator==(const Matrix<T>& m2) const;
        bool operator!=(const Matrix<T>& m2) const;

        /// Non-member functions specfic to complex matrices. 
        /// TODO: Get the real part of a complex matrix

        /// TODO: Get the imag part of a complex matrix

        /// Useful matrix functions
        T det() const;          ///< Calculate the determinant of a diagonal matrix
        T trace() const;        ///< Returns the trace of a matrix
        void transposeIP(Matrix<T>& m);     ///< A function to transpose in place
        Matrix<T> transpose();  ///< Diagonalize a matrix
        void conjugateIP(Matrix<T>& m);     ///< A function to transpose in place
        Matrix<T> conjugate();  ///< Returns the complex conjugate of the matrix
        void eye();             ///< Sets this matrix as the identity
    
        /// Matrix algebra functions
        void diagonalize(Matrix<T>& eigvecs, Matrix<T>& eigvals) const; ///< Diagonalize a real matrix
        Matrix<T> dagger(); ///<  Provide  A(dag)

    private:
        void resize(int rows, int cols);

    
}; // end of class Matrix

/* there may be a more glamorous/efficient way of doing this.
this is here to fix linker errors which otherwise require the
 entire template class to be defined in this header file */
#include "Matrix.cpp"

#endif // MATRIX_H
