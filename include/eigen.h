#pragma once
#ifndef TCC_INCLUDE_EIGEN_H
#define TCC_INCLUDE_EIGEN_H

#pragma GCC diagnostic push
#pragma GCC system_header
#include <Eigen/Dense>
#pragma GCC diagnostic pop

template <typename T>
using RowMatrix
    = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

using RowMatrixXd = RowMatrix<double>;

template <typename T>
using eMap = Eigen::Map<T>;


#endif // TCC_INCLUDE_EIGEN_H
