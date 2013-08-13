#include "mpqc/tensor/base.hpp"
#include "mpqc/tensor.hpp"
#include "mpqc/tensor/permute.hpp"
#include <iostream>

using namespace mpqc;

int main() {

    static const int N = 3;

    double data[N*N*N] = { 0 };
    size_t dims[] = { N, N, N };

    //typedef mpqc::TensorRowMajor Order;
    typedef mpqc::TensorColumnMajor Order;

    mpqc::TensorBase<const double, N, Order> t(data, dims);
    //mpqc::TensorBase<double, N, Order> u(data, dims);
    mpqc::Tensor<double, N, Order> u(dims);

    //u(range(0,1), range(0,1)) = t(range(1,2), range(1,2));
    //u = (t);

    for (int k : range(0,N)) {
        for (int j : range(0,N)) {
            for (int i : range(0,N)) {
                u(i,j,k) = i+j*N;
                std::cout << boost::tie(i,j,k) << ":"
                          << u(i,j,k) << " "
                          << i+j*N << std::endl;
            }
        }
    }

    permute_in_place<0,2,1>(u);

    for (int k : range(0,N)) {
        for (int j : range(0,N)) {
            for (int i : range(0,N)) {
                std::cout << boost::tie(i,j,k) << ":" << u(i,j,k) << std::endl;
                assert(u(i,k,j) == i+j*N);
            }
        }
    }


    // for (int k : range(0,N)) {
    //     for (int j : range(0,N)) {
    //         for (int i : range(0,N)) {
    //             std::cout << boost::tie(i,j,k) << ":" << u(i,j,k) << " " << i+j*N << std::endl;
    //             assert(u(i,j,k) == i+j*N);
    //         }
    //     }
    // }

    // u += u;
    // u /= 1.5;
    // u *= 2.3;

    // s = u;

    // for (int i : range(0,N)) {
    //     for (int j : range(0,N)) {
    //         std::cout << boost::tie(i,j) << ": s = " << s(i,j) << std::endl;
    //     }
    // }
    
    return 0;

}
