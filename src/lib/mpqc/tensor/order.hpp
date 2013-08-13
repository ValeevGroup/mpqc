#ifndef MPQC_TENSOR_ORDER_HPP
#define MPQC_TENSOR_ORDER_HPP

#include <boost/tuple/tuple.hpp>
#include <boost/fusion/include/boost_tuple.hpp>
#include <boost/fusion/include/boost_array.hpp>
#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/include/pop_front.hpp>
#include <boost/fusion/include/front.hpp>
#include <boost/fusion/include/pop_back.hpp>
#include <boost/fusion/include/back.hpp>
#include <boost/mpl/int.hpp>

namespace mpqc {

    /// @addtogroup Tensor
    /// @{

    /// Tensor column major (i.e. first dimension is contiguous) storage order
    template<size_t N>
    struct TensorRowMajor {
        /// Constructs order object
        /// @param ld Leading dimensions
        TensorRowMajor(const size_t *ld) {
            size_t stride = 1;
            for (int i = 0; i < N; ++i) {
                strides_[N-(i+1)] = stride;
                stride *= ld[N-(i+1)];
                //printf("strides_[%i]=%i\n", N-(i+1), strides_[N-(i+1)]);
            }
        }
        /// given N-d index tuple, compute 1-d index
        template<class Index>
        ptrdiff_t index(const Index &idx) const {
            namespace fusion = boost::fusion;
            ptrdiff_t diff = index(fusion::as_vector(idx), this->strides_);
            //std::cout << idx << " diff=" << diff << std::endl;
            return diff;
        }
    private:
        size_t strides_[N];
    protected:
        template<class Index>
        static typename boost::enable_if_c<
            (boost::fusion::result_of::size<Index>::value > 1),
            ptrdiff_t
            >::type
        index(const Index &idx, const size_t *strides) {
            namespace fusion = boost::fusion;
            ptrdiff_t diff = ((index(fusion::pop_front(idx), strides+1)) +
                    (fusion::front(idx)*(*strides)));
            return diff;
        }
        template<class Index>
        static typename boost::enable_if_c<
            (boost::fusion::result_of::size<Index>::value == 1),
            ptrdiff_t
            >::type
        index(const Index &idx, const size_t *strides) {
            namespace fusion = boost::fusion;
            return fusion::front(idx);
        }
    };

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
            return index(boost::fusion::as_vector(idx), this->strides_, boost::mpl::int_<N-1>());
        }
    private:
        size_t strides_[N];
    protected:
        /// recursively add index for each rank K,...,0
        template<class Index, class Strides, int K>
        static ptrdiff_t index(const Index &idx, const Strides &strides, boost::mpl::int_<K>) {
            using boost::fusion::at_c;
            // std::cout << "strides_[K]=" << this->strides_[K] << std::endl;
            // std::cout << "base[K]=" << base[K] << std::endl;
            return (index(idx, strides, boost::mpl::int_<K-1>()) + (at_c<K>(idx)*strides[K]));
        }
        /// end recursion and return rank-0 index (first Index element)
        template<class Index, class Strides>
        static ptrdiff_t index(const Index &idx, const Strides &strides, boost::mpl::int_<0>) {
            using boost::fusion::at_c;
            //std::cout << base[0] + get<0>(idx) << std::endl;
            return at_c<0>(idx);
        }
    };

    /// @}

}

#endif /* MPQC_TENSOR_ORDER_HPP */
