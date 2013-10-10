#ifndef MPQC_RANGE_OPERATOR_HPP
#define MPQC_RANGE_OPERATOR_HPP

#include "mpqc/range.hpp"
#include "mpqc/range/tie.hpp"

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/enum_params.hpp>
#include <boost/preprocessor/tuple.hpp>

/// @addtogroup MathRange
/// @{ 

#define MPQC_RANGE_CONST_OPERATOR(Z, N, DATA)                                   \
        template<BOOST_PP_ENUM_PARAMS(N, class _R)>                             \
        BOOST_PP_TUPLE_ELEM(2, 0, DATA)                                         \
            operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, const _R, &r)) const {    \
            return BOOST_PP_TUPLE_ELEM(2, 1, DATA)                              \
                (mpqc::detail::range_tie                                      \
                 (boost::tie(BOOST_PP_ENUM_PARAMS(N, r))));                     \
        }

#define MPQC_RANGE_OPERATOR(Z, N, DATA)                                         \
        template<BOOST_PP_ENUM_PARAMS(N, class _R)>                             \
        BOOST_PP_TUPLE_ELEM(2, 0, DATA)                                         \
            operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, const _R, &r)) {          \
            return BOOST_PP_TUPLE_ELEM(2, 1, DATA)                              \
                (mpqc::detail::range_tie                                      \
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
