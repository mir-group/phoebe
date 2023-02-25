#ifndef TEST_UTILS_H
#define TEST_UTILS_H
#include <Eigen/Core>
#include <complex>

/** Utility function for comparing eigenvectors without worrying about
 * phases, subspace choice etc.
 *
 * @param A: nxn matrix
 * @param x: n-dim vector
 * @return: A*diag(x)*A^H
*/
template<class VecType, class MatType>
Eigen::MatrixXcd mat_vec_mat_adj(MatType A, VecType x, int n){
  Eigen::MatrixXcd out = Eigen::MatrixXcd::Zero(n, n);

  for(int j = 0; j < n; j++){
    for(int i = 0; i < n; i++){
      for(int k = 0; k < n; k++){
        out(i,j) += A(i,k) * x(k) * std::conj(A(j,k));
      }
    }
  }

  return out;
}

template<class TensorType>
Eigen::MatrixXcd ev3Dto2D(TensorType A){
  int n = A.dimension(2);
  Eigen::MatrixXcd B(n,n);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      B(i,j) = A(i/3, i%3, j);
    }
  }

  return B;
}
#endif
