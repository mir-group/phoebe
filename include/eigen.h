#ifndef THISEIGEN_H
#define THISEIGEN_H

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>

typedef Eigen::Matrix<long, Eigen::Dynamic, 1> VectorXl;
typedef Eigen::Matrix<long, Eigen::Dynamic, Eigen::Dynamic> MatrixXl;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
  Eigen::RowMajor> MatrixXdRowMajor;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic,
  Eigen::RowMajor> MatrixXcdRowMajor;

#endif
