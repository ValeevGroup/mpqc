#ifndef MPQC_TENSOR_HPP
#define MPQC_TENSOR_HPP

/// @defgroup Tensor mpqc.Tensor
/// Tensor objects and functions

#include "mpqc/tensor/ref.hpp"

namespace mpqc {

    /// Tensor reference class.
    /// All data is assumed to be contiguous.
    template<typename T, size_t N, class Order = TensorColumnMajor<N> >
    class Tensor : TensorRef<T,N,Order> {

    public:

        explicit Tensor(const size_t (&dims)[N])
            : TensorRef<T,N,Order>(allocate(dims), dims)
        { }

        Tensor(const Tensor &u)
            : TensorRef<T,N,Order>(allocate(u.dims.elems), u.dims.elems)
        {
            TensorRef<T,N,Order>::operator=(u);
        }

        // template<typename U>
        // Tensor(const TensorRef<U,N,Order> &u)
        //     : TensorRef<T,N,Order>(allocate(u.dims.elems), u.dims.elems)
        // {
        //     TensorRef<T,N,Order>::operator=(u);
        // }

        ~Tensor() {
            delete[] TensorRef<T,N,Order>::data();
        }

        Tensor& operator=(const Tensor &u) {
            TensorRef<T,N,Order>::operator=(u);
        }

    protected:

        static T* allocate(const size_t (&dims)[N]) {
            size_t size = 1;
            for (size_t i = 0; i < N; ++i) {
                size *= dims[i];
            }
            return new T[size];
        }
 
    };

    /// @}

}

#endif /* MPQC_TENSOR_HPP */
