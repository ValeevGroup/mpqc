#ifndef MPQC_TENSOR_BASE_HPP
#define MPQC_TENSOR_BASE_HPP

#include <algorithm>
#include <assert.h>

#include <boost/tuple/tuple.hpp>
#include <boost/fusion/include/boost_tuple.hpp>
#include <boost/mpl/accumulate.hpp>
#include <boost/mpl/print.hpp>

#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/remove_reference.hpp>

#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>

#include "mpqc/range.hpp"
#include "mpqc/range/tie.hpp"
#include "mpqc/utility/string.hpp"

#include "mpqc/tensor/order.hpp"

namespace mpqc {
namespace detail {
namespace Tensor {

    template<class T>
    struct is_integral {
        typedef typename boost::is_integral<
            typename boost::remove_reference<T>::type
            >::type type;
    };

    /// "returns" true if every element of tuple T is an integral type
    template<class T>
    struct is_integral_tuple {
        static const bool value =
            boost::mpl::accumulate<
            T,
            boost::mpl::bool_<true>,
            boost::mpl::and_< is_integral<boost::mpl::_2>, boost::mpl::_1 >
            >::type::value;
    };

    /// index tie wrapper
    template<class Tie>
    struct integral_tie : Tie {
        integral_tie(const Tie &t) : Tie(t) {}
    };

    /// range tie wrapper
    template<class Tie>
    struct range_tie : Tie {
        range_tie(const Tie &t) : Tie(t) {}
    };

    /// returns integral_tie if every type in Tie is an integral type,
    /// disabled otherwise
    template<class Tie>
    typename boost::enable_if<
        is_integral_tuple<Tie>,
        integral_tie<Tie>
        >::type
    tie(const Tie &t) {
        return integral_tie<Tie>(t);
    }

    /// returns range_tie if NOT every type in Tie is an integral type,
    /// disabled otherwise
    template<class Tie>
    typename boost::disable_if<
        is_integral_tuple<Tie>,
        range_tie<Tie>
        >::type
    tie(const Tie &t) {
        return range_tie<Tie>(t);
    }

} // Tensor
} // detail
} // mpqc


namespace mpqc {

    /// @addtogroup Tensor
    /// @{

    struct TensorIndexException : std::runtime_error {
        template<typename Rank, typename Index, typename Begin, typename End>
        TensorIndexException(Rank rank, Index index, Begin begin, End end)
            : std::runtime_error("rank<" + string_cast(rank) + "> " +
                                 "index=" + string_cast(index) + " " +
                                 "outside the dimensions [" +
                                 string_cast(begin) + ":" + string_cast(end) + ")")
        {}
    };

    struct TensorRangeException : std::runtime_error {
        template<typename Rank, typename Range, typename Begin, typename End>
        TensorRangeException(Rank rank, Range range, Begin begin, End end)
            : std::runtime_error("rank<" + string_cast(rank) + "> " +
                                 "range=" + string_cast(range) + " " +
                                 "outside the dimensions [" +
                                 string_cast(begin) + ":" + string_cast(end) + ")")
        {}
    };

    /// Tensor base class.
    /// For performance reasons the storage order needs to be known at compile time,
    /// for time being I assume col-major as default (Eigen default).
    /// Only the first major dimenstion is assumed to be contiguous.
    /// Only minimum functionality is provided, namely element and range operators
    template<typename T, size_t N, class Order = TensorColumnMajor<N> >
    class TensorBase {

    public:
        static const size_t RANK = N;
        typedef size_t Dims[N];

    protected:
        T *data_;
        Dims dims_;
        Order order_;

    public:

        TensorBase(T *data,
                   const size_t *dims,
                   const size_t *ld = NULL)
            : order_(ld ? ld : dims),
              data_(data)
        {
            std::copy(dims, dims+N, this->dims_);
        }

        size_t size() const {
            size_t size = 1;
            for (int i = 0; i < N; ++i)
                size *= this->dims_[i];
            return size;
        }

        const Dims& dims() const {
            return dims_;
        }

        // generate index operator of arity N
        // CV may be empty or const
#define MPQC_TENSOR_INDEX_OPERATOR(Z, N, CV)                                    \
        template< BOOST_PP_ENUM_PARAMS(N, class T) >                            \
        typename boost::enable_if                                               \
        < detail::Tensor::is_integral_tuple                                     \
          < boost::tuple<BOOST_PP_ENUM_PARAMS(N,T)> >,                          \
          CV T& >::type                                                         \
        operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, const T, &i)) CV {            \
            using detail::Tensor::tie;                                          \
            return this->operator()(tie(boost::tie(BOOST_PP_ENUM_PARAMS(N,i)))); \
        }                                                                       \

        // generate range operator of arity N
        // CV may be empty or const
#define MPQC_TENSOR_RANGE_OPERATOR(Z, N, CV)                                    \
        template< BOOST_PP_ENUM_PARAMS(N, class T) >                            \
        typename boost::disable_if                                              \
        < detail::Tensor::is_integral_tuple                                     \
          < boost::tuple<BOOST_PP_ENUM_PARAMS(N,T)> >,                          \
          TensorBase<CV T, N> >::type                                           \
        operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, const T, &i)) CV {            \
            using detail::Tensor::tie;                                          \
            return this->operator()(tie(boost::tie(BOOST_PP_ENUM_PARAMS(N,i)))); \
        }                                                                       \
        
        //BOOST_PP_REPEAT_FROM_TO(1, 5, MPQC_TENSOR_INDEX_OPERATOR, nil)
        MPQC_TENSOR_RANGE_OPERATOR(nil, 2,      )
        MPQC_TENSOR_RANGE_OPERATOR(nil, 2, const)

        MPQC_TENSOR_INDEX_OPERATOR(nil, 2,      )
        MPQC_TENSOR_INDEX_OPERATOR(nil, 2, const)

    protected:
        
        /// element-access operator
        template<class Tie>
        T& operator()(const detail::Tensor::integral_tie<Tie> &idx) {
            return this->data_[this->index(idx)];
        }

        /// element-access operator
        template<class Tie>
        const T& operator()(const detail::Tensor::integral_tie<Tie> &idx) const  {
            return this->data_[this->index(idx)];
        }

        template<class Tie>
        TensorBase<T, N>
        operator()(const detail::Tensor::range_tie<Tie> &tie) {
            return block< TensorBase<T, N> >(*this, tie);
        }

        template<class Tie>
        TensorBase<const T, N>
        operator()(const detail::Tensor::range_tie<Tie> &tie) const {
            return block< TensorBase<const T, N> >(*this, tie);
        }

    private:

        template<class Tie, int K>
        void check_index(const detail::Tensor::integral_tie<Tie> &tie,
                         boost::mpl::int_<K>) const {
            if ((tie.template get<K>() < 0) ||
                (tie.template get<K>() > this->dims_[K])) {
                throw TensorIndexException(K, tie.template get<K>(),
                                           0, this->dims_[K]);
            }
            check_index(tie, boost::mpl::int_<K+1>());
        }

        template<class Tie>
        void check_index(const detail::Tensor::integral_tie<Tie> &tie,
                         boost::mpl::int_<N>) const {
        }

        template<class Tie>
        ptrdiff_t index(const detail::Tensor::integral_tie<Tie> &idx) const {
            static_assert(boost::tuples::length<Tie>::value == N,
                          "Invalid TensorBase::operator() arity");
#ifndef NDEBUG
            check_index(idx, boost::mpl::int_<0>());
#endif
            ptrdiff_t index = order_.template index<Tie>(idx);
            //std::cout << idx << ":" << index << "->" << data_+index << std::endl;
            return index;
        }

    private:

        template<class U, class This, class Tie>
        static U block(This &t, detail::Tensor::range_tie<Tie> tie) {
            static_assert(boost::tuples::length<Tie>::value == N,
                          "Invalid TensorBase::operator() arity");
            boost::array<range,N> r = range::tie<Tie>(tie);
            boost::array<ptrdiff_t,N> begin;
            Dims dims;
            for (int i = 0; i < N; ++i) {
                begin[i] = *r[i].begin();
                dims[i] = r[i].size();
#ifndef  NDEBUG
                if ((*r[i].begin() < 0) || (*r[i].end() > t.dims_[i])) {
                    throw TensorRangeException(i, r[i], 0, t.dims_[i]);
                }
#endif
            }
            ptrdiff_t offset = t.order_.index(begin);
            //std::cout << offset << std::endl;
            return U(t.data_+offset, dims, t.order_);
        }

    private:

        friend class TensorBase< typename boost::remove_const<T>::type, N, Order>;
        TensorBase(T *data, const Dims &dims,
                   const Order &order)
            : data_(data), order_(order)
        {
            std::copy(dims, dims+N, dims_);
        }

        TensorBase& operator=(const TensorBase&);        

    };

    /// @}

}


#endif /* MPQC_TENSOR_BASE_HPP */
