#ifndef MPQC_TENSOR_BASE_HPP
#define MPQC_TENSOR_BASE_HPP

#include <algorithm>

#include <boost/tuple/tuple.hpp>
#include <boost/fusion/tuple.hpp>
#include <boost/mpl/accumulate.hpp>
#include <boost/mpl/print.hpp>

#include <boost/preprocessor/repetition/enum_params_with_a_default.hpp>

#include "mpqc/range/operator.hpp"

namespace mpqc {

    template<class T>
    struct all_integral {
        static const bool value =
            boost::mpl::accumulate<
            T,
            boost::mpl::bool_<true>,
            boost::mpl::and_< boost::is_integral<boost::mpl::_2>, boost::mpl::_1 >
            >::type::value;
    };

    template<class T, BOOST_PP_ENUM_PARAMS_WITH_A_DEFAULT(5, class T, int)>
    struct is_index {
        typedef boost::mpl::print<boost::tuple<T, BOOST_PP_ENUM_PARAMS(5,T)> > print;
        static const bool value = 
            all_integral< boost::tuple<T, BOOST_PP_ENUM_PARAMS(5,T)> >::value;
    };

    template<class Tie>
    struct index_tie : Tie {
        index_tie(const Tie &t) : Tie(t) {}
    };

    template<class Tie>
    typename boost::enable_if<
        all_integral<Tie>, index_tie<Tie>
        >::type
    tie(const Tie &t) {
        return index_tie<Tie>(t);
    }


    template<typename T, size_t N>
    class TensorBase {

    public:

        static const size_t RANK = N;

        TensorBase(T *data,
                   const size_t *dims,
                   const size_t *strides) {
            //this->order_ = TENSOR_STORAGE_UNDEFINED;
            this->data_ = data;
            std::fill(this->base_, this->base_+N, 0);
            std::copy(dims, dims+N, this->dims_);
            std::copy(strides, strides+N, this->strides_);
        }

        MPQC_RANGE_OPERATORS(5, TensorBase, this->operator())

#define MPQC_TENSOR_INDEX_OPERATOR(Z, N, DATA)                                  \
        template< BOOST_PP_ENUM_PARAMS(N, class T) >                            \
        typename boost::enable_if                                               \
        < is_index<BOOST_PP_ENUM_PARAMS(N,T)>, T&>::type                        \
        operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, const T, &i)) {               \
            return this->operator()(tie(boost::tie(BOOST_PP_ENUM_PARAMS(N,i)))); \
        }                                                                       \
        
        //BOOST_PP_REPEAT_FROM_TO(1, 5, MPQC_TENSOR_INDEX_OPERATOR, nil)
        MPQC_TENSOR_INDEX_OPERATOR((), 2, ())

    protected:
 
        T *data_;
        size_t base_[N];
        size_t dims_[N];
        size_t strides_[N];

    protected:
        
        /// element-access operator
        template<class Tie>
        T& operator()(const index_tie<Tie> &idx);//  {
        //     return this->data_[this->index(idx)];
        // }

    private:

        template<class Index>
        size_t index(const Index &idx, boost::mpl::int_<0>) const {
            return 0;
        }

        template<class Index, int K>
        size_t index(const Index &idx, boost::mpl::int_<K>) const {
            static const int J = K-1;
            return ((this->base_[J] + boost::fusion::at_c<J>(idx))*this->strides_[J] +
                    (index(idx, boost::mpl::int_<J>())));
        }

        // template<class Index>
        // size_t index(const Index &idx) const {
        //     static_assert(boost::fusion::result_of::size<Index>::value == N, "");
        //     size_t i = index(boost::fusion::as_vector(idx), boost::mpl::int_<N>());
        //     //std::cout << idx << "->" << i << std::endl;
        //     return i;
        // }

        // template<class U, class This, class Tie>
        // static U generate(This &t, Operator::RangeTie<Tie> r) {
        //     static_assert(Operator::RangeTie<Tie>::N == RANK,
        //                   "Tensor rank mismatch in range operator()");
        //     U u(t);
        //     for (int i = 0; i < N; ++i) {
        //         const range &r = r[i];
        //         u.dims_[i] = r[i].size();
        //         u.base_[i] += *(range_cast(r).begin();
        //         u.strides_[i] = u.strides_[i];
        //         //printf("generate: %lu %lu %lu\n", t.dims_[i], t.base_[i], t.strides_[i]);
        //     }
        //     return ;
        // }

    };



}


#endif /* MPQC_TENSOR_BASE_HPP */
