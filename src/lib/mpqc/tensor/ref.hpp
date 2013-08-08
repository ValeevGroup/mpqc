#ifndef MPQC_TENSOR_REF_HPP
#define MPQC_TENSOR_REF_HPP

#include "mpqc/tensor/base.hpp"

namespace mpqc {

    /// @addtogroup Tensor
    /// @{

    /// Tensor reference class.
    /// All data is assumed to be contiguous.
    template<typename T, size_t N, class Order = TensorColumnMajor<N> >
    class TensorRef : TensorBase<T,N,Order> {

    public:

        TensorRef(T *data, const size_t (&dims)[N])
            : TensorBase<T,N,Order>(data, dims) {}

        T* data() { return this->data_; }
        const T* data() const { return this->data_; }

        /// fast assignment operator
        template<typename U>
        TensorRef& operator=(const TensorRef<U,N,Order> &u) {
            size_t size = this->size();
            std::copy(u.data(), u.data()+size, this->data());
        }
 
    };

    /// @}

}


#endif /* MPQC_TENSOR_REF_HPP */
