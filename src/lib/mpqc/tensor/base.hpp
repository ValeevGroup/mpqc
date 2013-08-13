#ifndef MPQC_TENSOR_BASE_HPP
#define MPQC_TENSOR_BASE_HPP

#include <algorithm>
#include <assert.h>

#include <boost/array.hpp>
#include <boost/tuple/tuple.hpp>

#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>

#include "mpqc/range.hpp"
#include "mpqc/range/tie.hpp"
#include "mpqc/utility/string.hpp"

#include "mpqc/tensor/forward.hpp"
#include "mpqc/tensor/functional.hpp"
#include "mpqc/tensor/exception.hpp"

namespace mpqc {

    /// @addtogroup Tensor
    /// @{

    /// Tensor base class.
    /// For performance reasons the storage order needs to be known at compile time,
    /// for time being I assume col-major as default (Eigen default).
    /// Only the first major dimenstion is assumed to be contiguous.
    /// Only minimum functionality is provided, namely element and range operators
    template<typename T, size_t N, class Order>
    class TensorBase {

    public:
        static const size_t RANK = N;
        typedef boost::array<size_t,N> Dims;
        typedef boost::array<size_t,N> Strides;

    protected:
        T *data_;
        Dims dims_;
        Strides strides_;

    public:

        TensorBase(T *data,
                   const size_t *dims,
                   const size_t *ld = NULL)
        {
            this->data_ = data;
            std::copy(dims, dims+N, this->dims_.begin());
            strides_ = Order::template strides<N>(ld ? ld : dims);
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

    public:
        
        template<typename U, class O>
        void operator=(const TensorBase<const U,N,O> &u) {
            detail::Tensor::apply(detail::Tensor::assign(), *this, u);
        }

        void operator=(const TensorBase& o) {
            this->operator=<T>(o);
        }

    public:

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
          TensorBase<CV T, N, Order> >::type                                    \
        operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, const T, &i)) CV {            \
            using detail::Tensor::tie;                                          \
            return this->operator()(tie(boost::tie(BOOST_PP_ENUM_PARAMS(N,i)))); \
        }                                                                       \
        
        BOOST_PP_REPEAT_FROM_TO(1, 5, MPQC_TENSOR_INDEX_OPERATOR,      )
        BOOST_PP_REPEAT_FROM_TO(1, 5, MPQC_TENSOR_INDEX_OPERATOR, const)
        BOOST_PP_REPEAT_FROM_TO(1, 5, MPQC_TENSOR_RANGE_OPERATOR,      )
        BOOST_PP_REPEAT_FROM_TO(1, 5, MPQC_TENSOR_RANGE_OPERATOR, const)

    public:
        
        /// element-access operator
        template<class Seq>
        T& operator()(const detail::Tensor::integral_tie<Seq> &idx) {
            return this->data_[this->index(idx)];
        }

        /// element-access operator
        template<class Seq>
        const T& operator()(const detail::Tensor::integral_tie<Seq> &idx) const  {
            return this->data_[this->index(idx)];
        }

        template<class Seq>
        TensorBase<T, N, Order>
        operator()(const detail::Tensor::range_tie<Seq> &tie) {
            return block< TensorBase<T, N, Order> >(*this, tie);
        }

        template<class Seq>
        TensorBase<const T, N, Order>
        operator()(const detail::Tensor::range_tie<Seq> &tie) const {
            return block< TensorBase<const T, N, Order> >(*this, tie);
        }

    private:

        template<class Seq, int K>
        void check_index(const detail::Tensor::integral_tie<Seq> &tie,
                         boost::mpl::int_<K>) const {
            using boost::fusion::at_c;
            if ((at_c<K>(tie) < 0) || (at_c<K>(tie) > this->dims_[K])) {
                throw TensorIndexException(K, at_c<K>(tie), 0, this->dims_[K]);
            }
            check_index(tie, boost::mpl::int_<K+1>());
        }

        template<class Seq>
        void check_index(const detail::Tensor::integral_tie<Seq> &tie,
                         boost::mpl::int_<N>) const {}

        template<class Seq>
        ptrdiff_t index(const detail::Tensor::integral_tie<Seq> &idx) const {
            static_assert(boost::fusion::result_of::size<Seq>::value == N,
                          "Invalid TensorBase::operator() arity");
// #ifndef NDEBUG
//             check_index(idx, boost::mpl::int_<0>());
// #endif
            ptrdiff_t index = Order::index((const Seq&)idx, this->strides_);
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
            ptrdiff_t offset = Order::index(begin, t.strides_);
            //std::cout << offset << std::endl;
            return U(t.data_+offset, dims, t.strides_);
        }

    private:

        friend class TensorBase< typename boost::remove_const<T>::type, N, Order>;
        TensorBase(T *data, const Dims &dims, const Strides &strides) 
            : data_(data), dims_(dims), strides_(strides)
        {
        }

    };

    /// @}

}


#endif /* MPQC_TENSOR_BASE_HPP */
