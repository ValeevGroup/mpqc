#include "mpqc/tensor.hpp"
#include <iostream>

#include <boost/mpl/vector.hpp>

using namespace mpqc;

int main() {

    static const int N = 2;

    double data[N*N] = { 0 };
    size_t dims[] = { N, N, N };

    typedef mpqc::TensorColumnMajor<3> Order;

    mpqc::TensorBase<const double, 3, Order> t(data, dims);
    mpqc::TensorRef<double, 3, Order> u(data, dims);
    mpqc::Tensor<double, 3, Order> s(dims);

    //u(range(0,1), range(0,1), range(0,1)) = t;
    //u = (t);

    for (int k : range(0,N)) {
        for (int j : range(0,N)) {
            for (int i : range(0,N)) {
                u(i,j,k) = i+j*N;
            }
        }
    }

    for (int k : range(0,N)) {
        for (int j : range(0,N)) {
            for (int i : range(0,N)) {
                std::cout << boost::tie(i,j,k) << ":" << u(i,j,k) << " " << i+j*N << std::endl;
                assert(u(i,j,k) == i+j*N);
            }
        }
    }

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
