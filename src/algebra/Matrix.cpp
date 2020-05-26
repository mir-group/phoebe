#include "Matrix.h"
#include <iostream>

// TODO: determine a better way to handle errors
// TODO: what should I do about resize functionality/ 

/* ------------- Constructors -------------- */
// A default constructor to build a dense matrix of zeros to be filled
template <typename T> Matrix<T>::Matrix(int nRows, int nCols){
    this->nRows = nRows;
    this->nCols = nCols;
    mat = 0; // intialize internal pointer
        
    // TODO: allocate memory for internal array. This is a potential point where we could 
    // put instructions regarding cpu/gpu memory
    mat = new T[nRows*nCols];
    for(int i=0; i<nRows*nCols; i++) mat[i] = 0; // fill with zeroes
    assert(mat == nullptr); // Memory could not be allocated, end program
}
// default constructor
template <typename T> Matrix<T>::Matrix() { 
   mat = 0; 
   nRows = 0; 
   nCols = 0; 
}
// copy constructor 
template <typename T> Matrix<T>::Matrix(const Matrix<T>& toCopy) { 
   mat = new T[toCopy.getSize()];
   copy(toCopy); 
} 

/* ------------- Very basic operations -------------- */
template <typename T> int Matrix<T>::numRows() const { return nRows; }
template <typename T> int Matrix<T>::numCols() const { return nCols; }
template <typename T> int Matrix<T>::getSize() const { return nCols*nRows; }
template <typename T> void Matrix<T>::print() const { 
        for(int row = 0; row<nRows; row++){
                for(int col = 0; col<nCols; col++){
                        std::cout << mat[index(row,col)] << "\t" << std::endl;
                }  
                std::cout << "\n" << std::endl;
        } 
        std::cout << mat << "\n" << std::endl;
}
// copy the elements of another matrix into this one
template <typename T> void Matrix<T>::copy(const Matrix<T>& toCopy)  { 
   T* p = mat + nRows*nCols; 
   T* q = toCopy.mat + toCopy.getSize(); 
   while( p > mat) { 
       *--p = *--q; // TODO: is it better to use std::swap here? 
   } 
}
// indexing puts the matrix in column major format
template <typename T> int Matrix<T>::index(int row,int col) const { return nRows * col + row;} 

/* ------------- Some operators -------------- */
// Get element
template <typename T> T Matrix<T>::operator()(const int row, const int col) const { return mat[index(row,col)]; }
// Set element
template <typename T> void Matrix<T>::operator()(const int row, const int col, const T value){ mat[index(row,col)] = value; }

// Get a block of the matrix
template <typename T>
Matrix<T> Matrix<T>::operator()(const int rowStart, const int rowStop, const int colStart, const int colStop) const{
    assert( (rowStop <= nRows) && (colStop <= nCols)); //, "Block is out of bounds." );
    Matrix<T> c(rowStop-rowStart,colStop-colStart);
    for(int row = rowStart; row<=rowStop; row++) {
        for(int col = colStart; col<=colStop; col++) {
                c(row-rowStart,col-colStart) = this(row,col);
        }
    } 
    return c;
}
// Get a row of the matrix
template <typename T> Matrix<T> Matrix<T>::getRow(const int row) const {
    assert( (row < nRows)); //, "Row is out of bounds." );
    Matrix<T> c(1,nCols);
    for(int i =0; i<nCols; i++) c(0,i) = this(row,i);
    return c;
}
// Get a col of the matrix
template <typename T> Matrix<T> Matrix<T>::getCol(const int col) const {
    assert( (col < nCols)); //, "Col is out of bounds." );
    Matrix<T> c(nRows,1);
    for(int i =0; i<nRows; i++) c(i,0) = this(i,col);
    return c;
}
// Set a block of the matrix
template <typename T>
void Matrix<T>::operator()(const int rowStart, const int rowStop, const int colStart, const int colStop, const Matrix<T>& value){
    assert( (rowStop - rowStart == value.numRows()) && (colStop - colStart == value.numCols()) ); //"Matrix to set is larger or smaller than the block requested." );
    assert( (rowStop <= nRows) && (colStop <= nCols)); // "Block is out of bounds." );

    for(int row = rowStart; row<=rowStop; row++) {
        for(int col = colStart; col<=colStop; col++) {
                this(row,col) = value(row-rowStart,col-colStart);
        }
    }
}
// Set a row of the matrix
template <typename T> void Matrix<T>::setRow(const int row, const Matrix<T>& value) {
    assert( (row < nRows)); //, "Row is out of bounds." );
    assert( (value.numRows() == 1) && (value.numCols() == nCols) );  // given matrix must match dims
    for(int i = 0; i<nCols; i++) this(row,i) = value(0,i);
}
// Set a col of the matrix
template <typename T> void Matrix<T>::setCol(const int col, const Matrix<T>& value) {
    assert( (col < nCols)); //, "Col is out of bounds." );
    assert( (value.numCols() == 1) && (value.numRows() == nRows) );  // given matrix must match dims
    for(int i = 0; i<nRows; i++) this(i,col) = value(i,0);
    mat.col(col) = value.mat;
}

// General unary negation TODO: dscal/zscal versions for double/complex double
template <typename T> Matrix<T> Matrix<T>::operator-() const{
    Matrix<T> c(nRows,nCols);
    for(int row = 0; row<nRows; row++) {
        for(int col = 0; col<nCols; col++) c(row,col) = -this(row,col);
    }
    return c;
}
/// Equivalence operators
template <typename T> bool Matrix<T>::operator==(const Matrix<T>& m1) const {
    if ( (nRows != m1.numRows() ) || ( nCols != m1.numCols() )) return false;
    for(int s =0; s< getSize(); s++) { if (mat[s] != m1[s]) return false; }
    return true; 
}
template <typename T> bool Matrix<T>::operator!=(const Matrix<T>& m1) const { 
    return !(this == m1);
}
// Explict specialization of BLAS matrix-matrix mult for Matrix<complex<double>>
/**template <> friend
Matrix<std::complex<double>> operator*(const Matrix<std::complex<double>>& m1, const Matrix<std::complex<double>>& m2) { 
    assert(m1.numCols() == m2.numRows());
    Matrix<std::complex<double>> ret(m1.numRows(),m2.numCols()); // newly sized matrix 
    // throw away variables
    Matrix<std::complex<double>> temp(*this); // copy 
    char transa = 'n'; char transb = 'n';
    std::complex<double> alpha(1.0,1.0);
    std::complex<double> beta(0.0,0.0); 
    zgemm_(transa,transb,m1.numRows(),m2.numCols(),m1.numCols(),alpha,m1.mat,m1.numRows(), m2.mat,m2.numRows(),beta,ret.mat,m1.numRows());  
}

// Explicit specializiation of BLAS matrix-matrix mult for Matrix<double>
template <> friend 
Matrix<double> operator*(const Matrix<double>& m1, const Matrix<double>& m2) {
    assert(m1.numCols() == m2.numRows());
    Matrix<double> ret(m1.numRows(),m2.numCols()); // newly sized matrix 
    // throw away variables
    Matrix<double> temp(*this); // copy 
    char transa = 'n'; char transb = 'n'; // compute only reigs
    double alpha = 1.0;
    double beta = 0.0;
    dgemm_(transa,transb,m1.numRows(),m2.numCols(),m1.numCols(),alpha,m1.mat,m1.numRows(), m2.mat,m2.numRows(),beta,ret.mat,m1.numRows());
} 
**/

/* ------------- Basic matrix functions -------------- */

template <typename T> T Matrix<T>::trace() const { 
    assert(nRows == nCols); // no trace for a non-square matrix
    T trace = 0;
    for(int s = 0; s<nRows; s++)  trace += this(s,s);
    return trace;
 }
// Determinant using LU decomp (see pg 52 in Numerical Recipes)
// NOTE: This isn't const because the throw away variables defined below 
// are automatically const if a member function is const. 
template <> double Matrix<double>::det() { 
    assert(nRows == nCols); // no det if non-square matrix
    int info = 0;  
    Matrix<double> lu(*this); // need to make a copy of our data array
    std::vector<int> pivot(nRows); // throw away argument

    // get the LU decomp of the matrix, LAPACK routine for complex val array
    dgetrf_(&nRows, &nCols, lu.mat, &nRows, pivot.data(), &info);
    assert(info>0); // if this not true, LU decomp failed somehow  

    // Product of the diagonal elements = det
    double det = 0; 
    for (int i = 0; i<nRows; i++) det *= lu(i,i);

    //TODO: if det is zero (info > 0) do we want to throw a warning?
    return det;
}
// Explicit specialization of determinant for a complex matrix. 
// are automatically const if a member function is const. 
template <> std::complex<double> Matrix<std::complex<double>>::det() {
    assert(nRows == nCols); // no det if non-square matrix
    int info = 0;
    Matrix<std::complex<double>> lu(*this); // need to make a copy of our data array
    std::vector<int> pivot(nRows); // throw away argument

    // get the LU decomp of the matrix, LAPACK routine for complex val array
    zgetrf_(&nRows, &nCols, lu.mat, &nRows, pivot.data(), &info);
    assert(info>0); // if this not true, LU decomp failed somehow  

    // Product of the diagonal elements = det
    std::complex<double> det(0.0,0.0);
    for (int i = 0; i<nRows; i++) det *= lu(i,i);

    //TODO: if det is zero (info > 0) do we want to throw a warning?
    return det;
}

// transpose in place
template <typename T> void Matrix<T>::transposeIP(Matrix<T>& m){
    for(int row = 0; row<nRows; row++) {
        for(int col = 0; col<row; col++) std::swap(m(row, col), m(col, row));
    }
}
/// Tranpose returning a new matrix 
template <typename T> Matrix<T> Matrix<T>::transpose(){
     Matrix<T> ret(*this); // make a copy to return
     transposeIP(*ret);
     return ret; 
}
// Conjugate in place for complex types
template <> void Matrix<std::complex<double>>::conjugateIP(Matrix<std::complex<double>>& m){
    // scale imag and real parts separately 
    // TODO: complex array is stored interleaved in memory, skip by 2s
    //dscal_(m.getSize(), -1.0, m.mat+1, 2); //imag parts
    for(int s = 0; s<m.getSize(); s++) m.mat[s].imag(-1.0*m.mat[s].imag());
} 
/// Conjugate in place for any real matrix (just returns)
template <typename T> void Matrix<T>::conjugateIP(Matrix<T>& m){ return; }

/// Conjugate returning a new matrix 
template <typename T> Matrix<T> Matrix<T>::conjugate(){
     Matrix<T> ret(*this); // make a copy to return
     conjugateIP(*ret); 
     return ret; 
}

//TODO: check that this produces a matrix of type T and not int 
template <typename T> void Matrix<T>::eye(){
    assert(nRows == nCols); 
    for(int row = 0; row<nRows; row++) this(row,row) = 1;
}
template <typename T> Matrix<T> Matrix<T>::dagger(){
    Matrix<T> ret(*this); // make a copy to operate on
    transposeIP(*ret);
    conjugateIP(*ret);
    return ret;    
}

/* ------------- Linear algebra functions -------------- */
// this one is only called if the matrix is a <Scalar> type.
// TODO: should the arguments be already sized matrices, or just 
//  pointers to empty, uninitialized ones? 

// dsyev -- real symmetric matrix eigenvalues and vectors
// *geev -- computes eigenvals, and left + right eigenvecs
// *heev -- computes eigenvals + eigenvecs of Hermitian matrix 

/// Diagonalize a complex double matrix 
// NOTE: This isn't const because the throw away variables defined below 
// are automatically const if a member function is const. 
template <> 
void Matrix<std::complex<double>>::diagonalize(Matrix<std::complex<double> >& eigvecs, Matrix<std::complex<double> >& eigvals){
    assert( nRows == nCols );   // needs to be square
    assert( (this->det().imag() != 0) && (this->det().real() != 0));   // should not be singular

    // throw away variables
    Matrix<std::complex<double>> temp(*this); // copy 
    char leig = 'N'; char reig = 'V'; // compute only reigs
    int lwork = (65*nRows); // Supposedly a good choice, I've also seen 20*N.  
    std::complex<double>* work = new std::complex<double>[lwork];
    double* rwork = new double[2*nRows]; 
    int info = 0; 
    // TODO: set this up so it calls zheev if the matrix is hermitian. 
    zgeev_(&leig, &reig, &nRows, temp.mat, &nRows, eigvals.mat, eigvecs.mat, &nRows, eigvecs.mat, &nRows, work, &lwork, rwork, &info);

    assert(info==0); // if it doesn't =0, there was an error. Different errors for info< or > 0.   
}
// Diagonalize for real double matrix 
// NOTE: This isn't const because the throw away variables defined below 
// are automatically const if a member function is const. 

// TODO: we need to change this dgeev call to match the actual call definition. 
// We also need to add functionality for setting the real and imag parts of a complex matrix, 
// so that we can fill the return eigvs 
template <>
void Matrix<double>::diagonalize(Matrix<std::complex<double> >& eigvecs, Matrix<std::complex<double> >& eigvals) {
    assert( nRows == nCols );   // needs to be square
    assert( this->det() != 0);   // should not be singular

    // throw away variables
    Matrix<double> temp(*this); // copy of input matrix which will be overwritten 
    char leig = 'N'; char reig = 'V'; // compute only reigs
    double* WR = new double[nRows]; // real part of eigvals
    double* WI = new double[nRows]; // imaginary part of eigvals 
    double* VL = new double[nRows*nRows]; // throw away for left eigenvectors
    double* VR = new double[nRows*nRows]; // holds eigenvectors
    int lwork = (65*nRows); // should be at least 4*nRows
    double* work = new double[lwork];
    int info = 0;
    // TODO: set this up so it calls dheev if the matrix is hermitian. 
    dgeev_(&leig, &reig, &nRows, temp.mat, &nRows, WR, WI, VL, &nRows, VR, &nRows, work, &lwork, &info);

    assert(info==0); // if it doesn't =0, there was an error. Different errors for info< or > 0.   
}
/// Calls the version for Matrix<double> after casting to double. 
// I don't think this will ever be used, but it's here to be more complete
// NOTE: This isn't const because the throw away variables defined below 
// are automatically const if a member function is const. 
template <typename T>
void Matrix<T>::diagonalize(Matrix<std::complex<double> >& eigvecs, Matrix<std::complex<double> >& eigvals){
    Matrix<double> temp(nRows,nCols); 
    // TODO: can we throw a warning about this type cast? 
    for(int s = 0; s<getSize(); s++) { temp.mat[s] = (T)*this.mat[s]; } 
    // call the regular diag for doubles 
    temp.diagonalize(eigvecs, eigvals); 
}

