
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
    
    public:
        // NOTE: I don't like leaving this as public. However, don't know another
        // way to handle diagonalization with passed pointers
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> mat;
    
        /// Matrix class constructors -----------------------------------
        Matrix(); // Default constructor -- creates a size 0 matrix
        Matrix(const int rows, const int cols); // Construct using row, col numbers
        //~Matrix();   // currently relying on the default destructor instead
        /// Copy constructor
        Matrix(const Matrix& m1) {
            mat = m1.mat;
            resize(m1.numRows(), m1.numCols());
        }
        
        /// Matrix class methods -----------------------------------
        int numRows() const;
        int numCols() const;
        /// Print the matrix
        void print() const;
        
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
        template<typename U>
        Matrix<T>& operator+=(const Matrix<U>& m1) {
            assert( (m1.numRows() == nRows) && (m1.numCols() == nCols) );
            static_assert( std::is_same<T, typename std::common_type<U,T>::type >::value, "You cannot += against a type which would return a type other than the lhs value." );
            mat = mat + m1.mat.template cast<T>();
            return *this;
        }
        template<typename U>
        Matrix<T>& operator-=(const Matrix<U>& m1) {
            assert( (m1.numRows() == nRows) && (m1.numCols() == nCols) );
            static_assert( std::is_same<T, typename std::common_type<U,T>::type >::value, "You cannot -= against a type which would return a type other than the lhs value." );
            mat = mat - m1.mat.template cast<T>();
            return *this;
        }
        template<typename U>
        Matrix<T>& operator*=(const Matrix<U>& m1) {
            assert(nRows == m1.numCols());
            static_assert( std::is_same<T, typename std::common_type<U,T>::type >::value, "You cannot *= against a type which would return a type other than the lhs value." );
            mat = mat * m1.mat.template cast<T>();  resize(nRows,m1.numCols());
            return *this;
        }
    
        /// Nonmember friend functions for operators like +, -, *. auto type here specifies the use of a "trailing return type" after the
        /// function delcaration
        template<typename U, typename V>
        friend auto operator+(const Matrix<U>& m1, const Matrix<V>& m2) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> {
            assert( (m1.numRows() == m2.numRows()) || m1.numCols() == m2.numCols() );
            Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> c;
            c.mat = m1.mat.template cast<typename std::common_type<decltype(U{}),decltype(V{})>::type>() + m2.mat.template cast<typename std::common_type<decltype(U{}),decltype(V{})>::type>();
            return c;
        }
        template<typename U, typename V>
        friend auto operator-(const Matrix<U>& m1, const Matrix<V>& m2) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> {
            assert( (m1.numRows() == m2.numRows()) || m1.numCols() == m2.numCols() );
            Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> c;
            c.mat = m1.mat.template cast<typename std::common_type<decltype(U{}),decltype(V{})>::type>() - m2.mat.template cast<typename std::common_type<decltype(U{}),decltype(V{})>::type>();
            return c;
        }
        /// Multiplication-based operations
        template<typename U, typename V>
        friend auto operator*(const Matrix<U>& m1, const Matrix<V>& m2) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type>  {
            assert(m1.numCols() == m2.numRows());
            Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> c;
            c.mat = m1.mat.template cast<typename std::common_type<decltype(U{}),decltype(V{})>::type>() * m2.mat.template cast<typename std::common_type<decltype(U{}),decltype(V{})>::type>(); c.resize(m2.numRows(),m1.numCols());
            return c;
        }
        // Matrix times a scalar.
        template<typename U, typename V>
        friend auto operator*( const Matrix<U>& m1,  const V scalar) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> {
            Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> c;
            c.mat = m1.mat.template cast<typename std::common_type<decltype(U{}),decltype(V{})>::type>() * static_cast<typename std::common_type<decltype(U{}),decltype(V{})>::type>(scalar);
            c.resize(m1.numRows(),m1.numCols());
            return c;
        }
        // Matrix times a scalar, have to define again with argument order flipped.
        template<typename U, typename V>
        friend auto operator*(const V scalar, const Matrix<U>& m1) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> {
            Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> c;
            c.mat = m1.mat.template cast<typename std::common_type<decltype(U{}),decltype(V{})>::type>() * static_cast<typename std::common_type<decltype(U{}),decltype(V{})>::type>(scalar);
            c.resize(m1.numRows(),m1.numCols());
            return c;
        }
        // Matrix / a scalar.
        template<typename U, typename V>
        friend auto operator/( const Matrix<U>& m1,  const V scalar) -> Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> {
            Matrix<typename std::common_type<decltype(U{}),decltype(V{})>::type> c;
            c.mat = m1.mat.template cast<typename std::common_type<decltype(U{}),decltype(V{})>::type>() / static_cast<typename std::common_type<decltype(U{}),decltype(V{})>::type>(scalar);
            c.resize(m1.numRows(),m1.numCols());
            return c;
        }
        
        /// Unary negation
        Matrix<T> operator-() const;
        /// Equivalence operators
        bool operator==(const Matrix<T>& m2) const { return mat==m2.mat; }
        bool operator!=(const Matrix<T>& m2) const { return mat!=m2.mat;  }

        /// Useful matrix functions
        T det() const;          ///< Calculate the determinant of a diagonal matrix
        T trace() const;        ///< Returns the trace of a matrix
        Matrix<T> transpose();  ///< Diagonalize a matrix
        Matrix<T> conjugate();  ///< Returns the complex conjugate of the matrix
        void eye();             ///< Sets this matrix as the identity
    
        /// Matrix algebra functions
        // NOTE: There are two overloaded definitions of "matrix::diagonalize". This is because
        // if you declare a matrix with a complex type, it causes an error otherwise.
        // I think for this reason, we will want either a complexMatrix or Hermitian type.
        void diagonalize(Matrix<T>& eigvecs, Matrix<T>& eigvals) const; ///< Diagonalize a real matrix
        void diagonalize(Matrix<std::complex<T> >& eigvecs, Matrix<std::complex<T> >& eigvals) const; ///< Diagonalize a complex matrix

        Matrix<T> dagger(); ///<  Provide  A(dag)

    private:
        void resize(int rows, int cols);
    
}; // end of class Matrix

/* there may be a more glamorous/efficient way of doing this.
this is here to fix linker errors which otherwise require the
 entire template class to be defined in this header file */
#include "Matrix.cpp"

#endif // MATRIX_H
