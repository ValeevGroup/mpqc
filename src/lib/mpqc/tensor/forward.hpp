#ifndef MPQC_TENSOR_FORWARD_HPP
#define MPQC_TENSOR_FORWARD_HPP

#include "mpqc/tensor/order.hpp"

namespace mpqc {

    template <typename T, size_t N, class Order = TensorColumnMajor<N> >
    struct TensorBase;

}

#endif /* MPQC_TENSOR_FORWARD_HPP */
