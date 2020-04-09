#include "Matrix.h"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

/* ------------- Constructors -------------- */
// A default constructor to build a dense matrix of zeros to be filled
template <typename T>
Matrix<T>::Matrix(int nRows, int nCols){
    this->nRows = nRows;
    this->nCols = nCols;
    mat = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::Zero(nRows,nCols);
}

template <typename T>
Matrix<T>::Matrix(){
    this->nRows = 1;
    this->nCols = 1;
    mat = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>();
    //mat = Eigen::Matrix<T,1,1>::Zero(1,1); // if the above causes a problem, swap with this line. Leaving here to see if we encounter a problem in use. 
}

/* ------------- Very basic operations -------------- */
template <typename T> int Matrix<T>::numRows() const { return mat.rows(); }
template <typename T> int Matrix<T>::numCols() const { return mat.cols(); }

template <typename T> void Matrix<T>::resize(int rows, int cols){
    this->nRows = rows;
    this->nCols = cols;
    mat.resize(rows,cols);
}
template <typename T> void Matrix<T>::print() const { std::cout << mat << "\n" << std::endl;
}

/* ------------- Some operators -------------- */
// Get element
template <typename T> T Matrix<T>::operator()(const int row, const int col) const { return mat(row,col); }
// Set element
template <typename T> void Matrix<T>::operator()(const int row, const int col, const T value){ mat(row,col) = value; }

// Get a block of the matrix
template <typename T>
Matrix<T> Matrix<T>::operator()(const int rowStart, const int rowStop, const int colStart, const int colStop) const{
    assert( (rowStop <= nRows) && (colStop <= nCols)); //, "Block is out of bounds." );
    Matrix<T> c;
    c.mat = this->mat.block(rowStart,colStart,rowStop,colStop);
    return c;
}
// Get a row of the matrix
template <typename T> Matrix<T> Matrix<T>::getRow(const int row) const {
    assert( (row < nRows)); //, "Row is out of bounds." );
    Matrix<T> c;
    c.mat = this->mat.row(row);
    return c;
}
// Get a col of the matrix
template <typename T> Matrix<T> Matrix<T>::getCol(const int col) const {
    assert( (col < nCols)); //, "Col is out of bounds." );
    Matrix<T> c;
    c.mat = this->mat.col(col);
    return c;
}
// Set a block of the matrix
template <typename T>
void Matrix<T>::operator()(const int rowStart, const int rowStop, const int colStart, const int colStop, const Matrix<T>& value){
    assert( (rowStop - rowStart == value.numRows()) && (colStop - colStart == value.numCols()) ); //"Matrix to set is larger or smaller than the block requested." );
    assert( (rowStop <= nRows) && (colStop <= nCols)); // "Block is out of bounds." );
    mat.block(rowStart,colStart,rowStop,colStop) = value.mat;
}
// Set a row of the matrix
template <typename T> void Matrix<T>::setRow(const int row, const Matrix<T>& value) {
    assert( (row < nRows)); //, "Row is out of bounds." );
    mat.row(row) = value.mat;
}
// Set a col of the matrix
template <typename T> void Matrix<T>::setCol(const int col, const Matrix<T>& value) {
    assert( (col < nCols)); //, "Col is out of bounds." );
    mat.col(col) = value.mat;
}

// Unary negation
template <typename T> Matrix<T> Matrix<T>::operator-() const{
    Matrix<T> c;
    c.mat = -1*this->mat;
    return c;
}
/* ------------- Basic  matrix functions -------------- */

template <typename T> T Matrix<T>::trace() const{ return mat.trace(); }
template <typename T> T Matrix<T>::det() const{ return mat.determinant(); }

template <typename T> Matrix<T> Matrix<T>::transpose(){
    Matrix<T> temp;
    temp.resize(this->nRows,this->nCols);
    temp.mat = this->mat;
    temp.mat.transposeInPlace(); // TODO: might be better to use ".transpose" but, then I can't return..
    return temp;
}
template <typename T> Matrix<T> Matrix<T>::conjugate(){
    Matrix<T> temp;
    temp.resize(this->nRows,this->nCols);
    temp.mat = this->mat.conjugate();
    //temp.mat.conjugate();
    return temp;
}
template <typename T> void Matrix<T>::eye(){
    mat = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::Identity(nRows,nCols);
; //(nRows, nCols);
}
template <typename T> Matrix<T> Matrix<T>::dagger(){
    Matrix<T> temp;
    temp.resize(this->nRows,this->nCols);
    temp.mat = this->mat.conjugate();
    temp.mat.transposeInPlace();
    return temp;
}

/* ------------- Linear algebra functions -------------- */
// this one is only called if the matrix is a <Scalar> type.
template <typename T>
void Matrix<T>::diagonalize(Matrix<std::complex<T> >& eigvecs, Matrix<std::complex<T> >& eigvals) const{
    Eigen::EigenSolver< Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > solver(this->mat.size());
    assert( nRows == nCols );   // needs to be square
    assert( mat.determinant() != 0);  // should not be singular
    solver.compute(this->mat);
    eigvecs.mat = solver.eigenvectors();
    eigvals.mat = solver.eigenvalues();
}
// this one is only called if the matrix is a complex<Scalar> type.
template <typename T>
void Matrix<T>::diagonalize(Matrix<T>& eigvecs, Matrix<T>& eigvals) const{
    Eigen::ComplexEigenSolver< Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > solver(this->mat.size());
    assert( nRows == nCols );   // needs to be square
    assert( mat.determinant() != 0);  // should not be singular
    solver.compute(this->mat);
    eigvecs.mat = solver.eigenvectors();
    eigvals.mat = solver.eigenvalues();
}
