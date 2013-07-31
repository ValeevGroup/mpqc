#ifndef MPQC_RANGE_OPERATOR_HPP
#define MPQC_RANGE_OPERATOR_HPP

#include "mpqc/range.hpp"
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/tuple.hpp>

#define MPQC_RANGE_OPERATOR_DIMENSIONS

namespace mpqc {

#define MPQC_RANGIFY1(Z, N, TEXT) r[N] = rangify1(BOOST_PP_CAT(TEXT, N))

#define MPQC_RANGIFY(Z, N, TEXT)                                                \
    template<BOOST_PP_ENUM_PARAMS(N, class R)>                                  \
    std::vector<range> rangify(BOOST_PP_ENUM_BINARY_PARAMS(N, const R, &r)) {   \
        std::vector<range> r(N);                                                \
        BOOST_PP_ENUM(N, MPQC_RANGIFY1, r);                                     \
        return r;                                                               \
    }
    
    BOOST_PP_REPEAT_FROM_TO(1, 5, MPQC_RANGIFY, nil)

#undef MPQC_RANGIFY
#undef MPQC_RANGIFY1

}

#define MPQC_RANGE_CONST_OPERATOR(Z, N, DATA)                                   \
        template<BOOST_PP_ENUM_PARAMS(N, class _R)>                             \
        BOOST_PP_TUPLE_ELEM(0, DATA)                                            \
            operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, const _R, &r)) const {    \
            return BOOST_PP_TUPLE_ELEM(1, DATA)                                 \
                (rangify(BOOST_PP_ENUM_PARAMS(N, r)));                          \
        }

#define MPQC_RANGE_OPERATOR(Z, N, DATA)                                         \
        template<BOOST_PP_ENUM_PARAMS(N, class _R)>                             \
        BOOST_PP_TUPLE_ELEM(0, DATA)                                            \
            operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, const _R, &r)) {          \
            return BOOST_PP_TUPLE_ELEM(1, DATA)                                 \
                (rangify(BOOST_PP_ENUM_PARAMS(N, r)));                          \
        }

#define MPQC_RANGE_OPERATORS(N, Type, Function)                                 \
    BOOST_PP_REPEAT_FROM_TO(1, N, MPQC_RANGE_OPERATOR, (Type, Function))

#define MPQC_RANGE_CONST_OPERATORS(N, Type, Function)                           \
    BOOST_PP_REPEAT_FROM_TO(1, N, MPQC_RANGE_CONST_OPERATOR, (Type, Function))


#endif // MPQC_RANGE_OPERATOR_HPP
