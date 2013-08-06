#include "mpqc/tensor/base.hpp"
#include <iostream>

using namespace mpqc;

int main() {
    double data[100];
    size_t dims[] = { 10, 10 };
    size_t strides[] = { 1, 10 };
    mpqc::TensorBase<double,2> t(data, dims, strides);

    t(0,0);
    std::cout << mpqc::is_index<int, int>::value << std::endl;

}
