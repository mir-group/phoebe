#include "gtest/gtest.h"
#include "PMatrix.h" 
#include <cmath>

TEST (PMatrixTest, diagonalize) { 
   
  int numRows = 2; 
  int numCols = 2; 

  // we do not specify numBlocksRows/Cols, 
  // as this matrix will have to be distributed using 
  // a square process grid. Not specifying them 
  // activates the default use case for a square process grid
  ParallelMatrix<double> pMat = ParallelMatrix<double>(numRows,numCols);

  // construct a pauli matrix in a very brute force way
  // let's first do sigma_x
  if(pMat.indicesAreLocal(0,1)) pMat(0,1) = 1.0;
  if(pMat.indicesAreLocal(1,0)) pMat(1,0) = 1.0;

  auto tup = pMat.diagonalize();
  auto eigenvalues = std::get<0>(tup); 
  auto eigenvectors = std::get<1>(tup);

  // check the eigenvalues, which are the same 
  // on each process
  double sumEigValues = 0.0;
  for(int i = 0; i < 2; i++) {
    sumEigValues += eigenvalues[i];
  }
  EXPECT_EQ(sumEigValues, 0);

  // check the eigenvectors, 
  // which themselves are distributed
  double sumEigVectors = 0.0;
  for(int i = 0; i < 2; i++) {
    for(int j = 0; j < 2; j++) {
      if(eigenvectors.indicesAreLocal(i,j)) {
        sumEigVectors += eigenvectors(i,j);
      }
    }  
  }
  mpi->allReduceSum(&sumEigVectors);
  EXPECT_EQ(sumEigVectors, 2*(1./sqrt(2)));

  // check that the original matrix is still intact
  if(pMat.indicesAreLocal(0,0)) { EXPECT_EQ(pMat(0,0), 0.0); }
  if(pMat.indicesAreLocal(1,0)) { EXPECT_EQ(pMat(1,0), 1.0); }
  if(pMat.indicesAreLocal(0,1)) { EXPECT_EQ(pMat(0,1), 1.0); }
  if(pMat.indicesAreLocal(1,1)) { EXPECT_EQ(pMat(1,1), 0.0); }

}
