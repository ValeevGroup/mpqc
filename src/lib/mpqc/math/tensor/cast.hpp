#ifndef MPQC_MATH_TENSOR_CAST_HPP
#define MPQC_MATH_TENSOR_CAST_HPP

#include "mpqc/math/tensor/ref.hpp"
#include "mpqc/math/tensor/order.hpp"
#include <Eigen/Dense>

namespace mpqc {

    template<int N, typename T>
    Eigen::Map< Eigen::Matrix<T, Eigen::Dynamic, 1> >
    vector_cast(TensorRef<T,N,TensorColumnMajor> &r) {
        int n = 1;
        for (int j = 0; j < N; ++j) {
            n *= r.dims()[j];
        }
        return Eigen::Matrix<T, Eigen::Dynamic, 1>::Map(r.data(), n);
    }

    template<int M, int N, typename T>
    Eigen::Map< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >
    matrix_cast(TensorRef<T,M+N,TensorRowMajor> &r) {
        int m = 1;
        int n = 1;
        for (int i = 0; i < M; ++i) {
            m *= r.dims()[i];
        }
        for (int j = 0; j < N; ++j) {
            n *= r.dims()[j+M];
        }
        return Eigen::Matrix<
            T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor
            >::Map(r.data(), m, n);
    }

    template<int M, int N, typename T>
    Eigen::Map< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> >
    matrix_cast(TensorRef<T,M+N,TensorColumnMajor> &r) {
        int m = 1;
        int n = 1;
        for (int i = 0; i < M; ++i) {
            m *= r.dims()[i];
        }
        for (int j = 0; j < N; ++j) {
            n *= r.dims()[j+M];
        }
        return Eigen::Matrix<
            T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor
            >::Map(r.data(), m, n);
    }

} // mpqc


#endif /* MPQC_MATH_TENSOR_CAST_HPP */
