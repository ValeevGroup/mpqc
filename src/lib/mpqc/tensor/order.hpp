#ifndef MPQC_TENSOR_ORDER_HPP
#define MPQC_TENSOR_ORDER_HPP

#include <boost/tuple/tuple.hpp>
#include <boost/fusion/include/boost_tuple.hpp>
#include <boost/fusion/include/boost_array.hpp>
#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/mpl/int.hpp>

namespace mpqc {

    /// @addtogroup Tensor
    /// @{

    /// Tensor column major (i.e. first dimension is contiguous) storage order
    template<size_t N>
    struct TensorColumnMajor {
        /// Constructs order object
        /// @param ld Leading dimensions
        TensorColumnMajor(const size_t *ld) {
            size_t stride = 1;
            for (int i = 0; i < N; ++i) {
                strides_[i] = stride;
                stride *= ld[i];
                //std::cout << "strides_[i]=" << strides_[i] << std::endl;
            }
        }
        /// given N-d index tuple, compute 1-d index
        template<class Index>
        ptrdiff_t index(const Index &idx) const {
            return index(boost::fusion::as_vector(idx), boost::mpl::int_<N-1>());
        }
    protected:
        /// recursively add index for each rank K,...,0
        template<class Index, int K>
        ptrdiff_t index(const Index &idx, boost::mpl::int_<K>) const {
            using boost::fusion::at_c;
            // std::cout << "strides_[K]=" << this->strides_[K] << std::endl;
            // std::cout << "base[K]=" << base[K] << std::endl;
            return (index(idx, boost::mpl::int_<K-1>()) +
                    (at_c<K>(idx)*this->strides_[K]));
        }
        /// end recursion and return rank-0 index (first Index element)
        template<class Index>
        ptrdiff_t index(const Index &idx, boost::mpl::int_<0>) const {
            using boost::fusion::at_c;
            //std::cout << base[0] + get<0>(idx) << std::endl;
            return at_c<0>(idx);
        }
    private:
        size_t strides_[N];
    };

    /// @}

}

#endif /* MPQC_TENSOR_ORDER_HPP */
