#ifndef MPQC_TENSOR_REF_HPP
#define MPQC_TENSOR_REF_HPP

#include "mpqc/tensor/base.hpp"
#include "mpqc/tensor/functional.hpp"

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

    public:

        /// fast assignment operator
        template<typename U>
        TensorRef& operator=(const TensorRef<U,N,Order> &u) {
            size_t size = this->size();
            std::copy(u.data(), u.data()+size, this->data());
        }

        /// fast += operator
        template<typename U>
        TensorRef& operator+=(const TensorRef<U,N,Order> &u) {
            this->apply(detail::Tensor::plus_assign(), u.data());
        }

        /// fast -= operator
        template<typename U>
        TensorRef& operator-=(const TensorRef<U,N,Order> &u) {
            this->apply(detail::Tensor::minus_assign(), u.data());
        }

        /// fast *= operator
        template<typename U>
        TensorRef& operator*=(const U &u) {
            this->apply(detail::Tensor::multiply_assign<U>(u));
        }

        /// fast /= operator
        template<typename U>
        TensorRef& operator/=(const U &u) {
            this->apply(detail::Tensor::divide_assign<U>(u));
        }

    protected:

        template<class F>
        F apply(F f) {
            size_t size = this->size();
            for (size_t i = 0; i < size; ++i) {
                f(this->data_[i]);
            }
            return f;
        }

        template<class F, typename Iterator>
        F apply(F f, Iterator it) {
            size_t size = this->size();
            for (size_t i = 0; i < size; ++i) {
                f(this->data_[i], *it++);
            }
            return f;
        }
 
    };

    /// @}

}


#endif /* MPQC_TENSOR_REF_HPP */
