#include "mpqc/math/tensor/base.hpp"
#include "mpqc/math/tensor.hpp"
#include "mpqc/math/tensor/permute.hpp"
#include "mpqc/math/tensor/cast.hpp"
#include "mpqc/math/blas.hpp"
#include "mpqc/math/matrix.hpp"
#include <iostream>

using namespace mpqc;

int main() {

    static const int N = 3;

    double data[N*N*N] = { 0 };
    size_t dims[] = { N, N, N };

    typedef mpqc::TensorRowMajor Order;
    //typedef mpqc::TensorColumnMajor Order;

    mpqc::TensorBase<const double, N, Order> t(data, dims);
    //mpqc::TensorBase<double, N, Order> u(data, dims);
    mpqc::Tensor<double, N, Order> u(dims);

    auto s = u(range(0,2), range(1,3), range(0,2));
    //u = (t);

    for (int k : range(0,2)) {
        for (int j : range(0,2)) {
            for (int i : range(0,2)) {
                u(i,j,k) = i+j*2+k*2*2;
                std::cout << boost::tie(i,j,k) << ":"
                          << u(i,j,k) << " "
                          << i+j*2+k*2*2 << std::endl;
            }
        }
    }

    //permute_in_place<0,2,1>(u);

    for (int k : range(0,2)) {
        for (int j : range(0,2)) {
            for (int i : range(0,2)) {
                std::cout << boost::tie(i,j,k) << ":" << u(i,j,k) << std::endl;
                MPQC_ASSERT(u(i,j,k) == i+j*2+k*2*2);
            }
        }
    }

    auto ut = matrix_cast<2,1>(u)*matrix_cast<1,2>(u);
    Matrix b(N,N);
    blas::gemm(1, matrix_cast<1,2>(u), matrix_cast<2,1>(u), 0, b);

    // for (int k : range(0,N)) {
    //     for (int j : range(0,N)) {
    //         for (int i : range(0,N)) {
    //             std::cout << boost::tie(i,j,k) << ":" << u(i,j,k) << " " << i+j*N << std::endl;
    //             MPQC_ASSERT(u(i,j,k) == i+j*N);
    //         }
    //     }
    // }

    // u += u;
    // u /= 1.5;
    // u *= 2.2;

    // s = u;

    // for (int i : range(0,N)) {
    //     for (int j : range(0,N)) {
    //         std::cout << boost::tie(i,j) << ": s = " << s(i,j) << std::endl;
    //     }
    // }
    
    return 0;

}
