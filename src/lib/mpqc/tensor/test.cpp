#include "mpqc/tensor/base.hpp"
#include <iostream>

#include <boost/mpl/vector.hpp>

using namespace mpqc;

int main() {

    static const int N = 5;

    double data[N*N] = { 0 };
    size_t dims[] = { N, N };
    mpqc::TensorBase<double,2> t(data, dims);

    mpqc::TensorBase<double,2> b = t(range(2,4), range(1,3));

    for (int i : range(0,3)) {
        for (int j : range(0,3)) {
            b(i,j) = i+j*10;
        }
    }

    for (int j : range(0,3)) {
        for (int i : range(0,3)) {
            //std::cout << boost::tie(i,j) << ":" << t(i,j) << std::endl;
            assert(t(i+2,j+1) == i+j*10);
        }
    }

}
