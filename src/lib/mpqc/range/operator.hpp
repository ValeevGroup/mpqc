#ifndef MPQC_RANGE_OPERATOR_HPP
#define MPQC_RANGE_OPERATOR_HPP

#include "mpqc/range.hpp"
#include <boost/tuple/tuple.hpp>
#include <boost/fusion/include/vector_tie.hpp>
#include <boost/fusion/include/make_vector.hpp>

#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/enum_params.hpp>
#include <boost/preprocessor/tuple.hpp>

namespace mpqc {

    /// boost::tuple tie wrapper
    template<class S>
    struct range::tie : S {
        static const int N = boost::tuples::length<S>::value;
        tie(const S &s) : S(s) {}
        /// cast range tuple as boost::array
        operator boost::array<range,N>() const {
            boost::array<range,N> a;
            cast(a, int_<N>());
        }
        /// cast range tuple as std::vector
        operator std::vector<range>() const {
            boost::array<range,N> a = *this;
            return std::vector<range>(a.begin(), a.end());
        }
    private:
        template<int I>
        struct int_ {};
        template<class V, int I>
        void cast(V &v, int_<I>) const {
            v[I-1] = range_cast(boost::tuples::get<I-1>(*this));
            cast(v, int_<I-1>());
        }
        template<class V>
        void cast(V &v, int_<0>) const {}
    };

}


namespace mpqc {
namespace detail {
namespace _range {

    // wrap boost::tuple T in range::tie
    template<class T>
    range::tie<T> tie(const T &t) {
        return range::tie<T>(t);
    }
    
}
}
}


/// @addtogroup Range
/// @{ 

#define MPQC_RANGE_CONST_OPERATOR(Z, N, DATA)                                   \
        template<BOOST_PP_ENUM_PARAMS(N, class _R)>                             \
        BOOST_PP_TUPLE_ELEM(0, DATA)                                            \
            operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, const _R, &r)) const {    \
            return BOOST_PP_TUPLE_ELEM(1, DATA)                                 \
                (mpqc::detail::_range::tie                                      \
                 (boost::tie(BOOST_PP_ENUM_PARAMS(N, r))));                     \
        }

#define MPQC_RANGE_OPERATOR(Z, N, DATA)                                         \
        template<BOOST_PP_ENUM_PARAMS(N, class _R)>                             \
        BOOST_PP_TUPLE_ELEM(0, DATA)                                            \
            operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, const _R, &r)) {          \
            return BOOST_PP_TUPLE_ELEM(1, DATA)                                 \
                (mpqc::detail::_range::tie                                      \
                 (boost::tie(BOOST_PP_ENUM_PARAMS(N, r))));                     \
        }

/** Generates <c>operators()(const R1 r1, ...)</c> of arity 1 to N.
    The parameters R are packed into mpqc::range::tie and passed to Function
    @param N Maximum operator arity
    @param Type Operator return type
    @param Function Function that accepts mpqc::range::tie and returns Type
 */
#define MPQC_RANGE_OPERATORS(N, Type, Function)                                 \
    BOOST_PP_REPEAT_FROM_TO(1, N, MPQC_RANGE_OPERATOR, (Type, Function))

/** Generates <c>const</c> version of MPQC_RANGE_OPERATORS */
#define MPQC_RANGE_CONST_OPERATORS(N, Type, Function)                           \
    BOOST_PP_REPEAT_FROM_TO(1, N, MPQC_RANGE_CONST_OPERATOR, (Type, Function))

/// @}

#endif // MPQC_RANGE_OPERATOR_HPP
