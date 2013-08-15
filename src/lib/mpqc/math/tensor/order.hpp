#ifndef MPQC_MATH_TENSOR_ORDER_HPP
#define MPQC_MATH_TENSOR_ORDER_HPP

#include <boost/tuple/tuple.hpp>
#include <boost/fusion/include/boost_tuple.hpp>
#include <boost/fusion/include/boost_array.hpp>
#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/include/pop_front.hpp>
#include <boost/fusion/include/front.hpp>
#include <boost/fusion/include/pop_back.hpp>
#include <boost/fusion/include/back.hpp>
#include <boost/fusion/include/accumulate.hpp>
#include <boost/mpl/int.hpp>

namespace mpqc {

    /// @addtogroup MathTensor
    /// @{

    /// Tensor row major (i.e. last dimension is contiguous) storage order
    struct TensorRowMajor {
        /// Constructs strides array
        /// @param ld Leading dimensions
        template<size_t N>
        static boost::array<size_t,N> strides(const size_t *ld) {
            boost::array<size_t,N> strides;
            size_t stride = 1;
            for (int i = 0; i < N; ++i) {
                strides[N-(i+1)] = stride;
                stride *= ld[N-(i+1)];
                //printf("strides_[%i]=%i\n", N-(i+1), strides_[N-(i+1)]);
            }
            return strides;
        }
        /// given N-d index tuple and strides, compute 1-d index
        template<class Index, class Strides>
        static ptrdiff_t index(const Index &idx, const Strides &strides) {
            namespace fusion = boost::fusion;
            const typename fusion::result_of::as_vector<Index>::type &v =
                fusion::as_vector(idx);
            ptrdiff_t diff = fusion::accumulate(fusion::pop_back(v),
                                                fusion::back(v),
                                                make_index(&strides[0]));
            return diff;
        }
    protected:
        struct make_index {
            typedef ptrdiff_t result_type;
            const size_t *strides;
            explicit make_index(const size_t *strides) : strides(strides) {}
            template<typename T>
            ptrdiff_t operator()(ptrdiff_t idx, const T& t) {
                return idx + t*(*strides++);
            }
        };
    };

    /// Tensor column major (i.e. first dimension is contiguous) storage order
    struct TensorColumnMajor {
        /// Constructs strides array
        /// @param ld Leading dimensions
        template<size_t N>
        static boost::array<size_t,N> strides(const size_t *ld) {
            boost::array<size_t,N> strides;
            size_t stride = 1;
            for (int i = 0; i < N; ++i) {
                strides[i] = stride;
                stride *= ld[i];
                //printf("strides_[%i]=%i\n", N-(i+1), strides_[N-(i+1)]);
            }
            return strides;
        }
        /// given N-d index tuple and strides, compute 1-d index
        template<class Index, class Strides>
        static ptrdiff_t index(const Index &idx, const Strides &strides) {
            namespace fusion = boost::fusion;
            const typename fusion::result_of::as_vector<Index>::type &v =
                fusion::as_vector(idx);
            ptrdiff_t diff = fusion::accumulate(fusion::pop_front(v),
                                                fusion::front(v),
                                                make_index(&strides[1]));
            return diff;
        }
    protected:
        struct make_index {
            typedef ptrdiff_t result_type;
            const size_t *strides;
            explicit make_index(const size_t *strides) : strides(strides) {}
            template<typename T>
            ptrdiff_t operator()(ptrdiff_t idx, const T& t) {
                return idx + t*(*strides++);
            }
        };
    };

    /// @}

}

#endif /* MPQC_MATH_TENSOR_ORDER_HPP */
