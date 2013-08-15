#ifndef MPQC_MATH_TENSOR_PERMUTE_HPP
#define MPQC_MATH_TENSOR_PERMUTE_HPP

#include "mpqc/math/tensor/forward.hpp"
#include "mpqc/math/tensor/exception.hpp"

#include <boost/mpl/int.hpp>
#include <boost/mpl/vector_c.hpp>

namespace mpqc {
 
    /// @addtogroup MathTensor
    /// @{

    template<typename T, class O>
    void permute_in_place(TensorBase<T,3,O> &t, boost::mpl::vector_c<int,1,0,2>) {
        const typename TensorBase<T,3,O>::Dims &dims = t.dims();
        if (dims[0] != dims[1]) {
            throw TensorDimensionsException(0, 1, dims[0], dims[1]);
        }
        for (int k = 0; k < dims[2]; ++k) {
            for (int j = 0; j < dims[1]; ++j) {
                for (int i = 0; i < j; ++i) {
                    std::swap(t(i,j,k), t(j,i,k));
                }
            }
        }
    }

    template<typename T, class O>
    void permute_in_place(TensorBase<T,3,O> &t, boost::mpl::vector_c<int,0,2,1>) {
        const typename TensorBase<T,3,O>::Dims &dims = t.dims();
        if (dims[1] != dims[2])
            throw TensorDimensionsException(1, 2, dims[1], dims[2]);
        for (int k = 0; k < dims[2]; ++k) {
            for (int j = 0; j < k; ++j) {
                for (int i = 0; i < dims[0]; ++i) {
                    std::swap(t(i,j,k), t(i,k,j));
                }
            }
        }
    }

    template<int I, int J, int K, typename T, class O>
    void permute_in_place(TensorBase<T,3,O> &t) {
        permute_in_place(t, boost::mpl::vector_c<int,I,J,K>());
    }

    /// @}

}

#endif /* MPQC_MATH_TENSOR_PERMUTE_HPP */
