#include "mpqc/tensor/base.hpp"
#include "mpqc/tensor/ref.hpp"
#include <iostream>

#include <boost/mpl/vector.hpp>

using namespace mpqc;

int main() {

    static const int N = 5;

    double data[N*N] = { 0 };
    size_t dims[] = { N, N };

    mpqc::TensorBase<double,2> t(data, dims);
    mpqc::TensorBase<double,2> b = t(range(2,4), range(1,8));
    mpqc::TensorRef<double,2> u(data, dims);

    for (int i : range(0,4)) {
        for (int j : range(0,3)) {
            b(i,j) = i+j*10;
        }
    }

    for (int j : range(0,3)) {
        for (int i : range(0,6)) {
            //std::cout << boost::tie(i,j) << ":" << t(i,j) << std::endl;
            assert(t(i+2,j+1) == i+j*10);
        }
    }

}
