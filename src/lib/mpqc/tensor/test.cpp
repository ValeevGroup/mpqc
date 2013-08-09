#include "mpqc/tensor.hpp"
#include <iostream>

#include <boost/mpl/vector.hpp>

using namespace mpqc;

int main() {

    static const int N = 5;

    double data[N*N] = { 0 };
    size_t dims[] = { N, N, N };

    mpqc::TensorBase<const double,3> t(data, dims);
    mpqc::TensorRef<double,3> u(data, dims);
    mpqc::Tensor<double,3> s(dims);

    u(range(0,1), range(0,1), range(0,1)) = t;
    //u = (t);

    // for (int j : range(0,4)) {
    //     for (int i : range(0,2)) {
    //         b(i,j) = i+j*10;
    //     }
    // }

    // for (int j : range(0,3)) {
    //     for (int i : range(0,2)) {
    //         //std::cout << boost::tie(i,j) << ":" << t(i,j) << std::endl;
    //         assert(t(i+2,j+1) == i+j*10);
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
    

}
