#ifndef MPQC_RANGE_OPERATOR_HPP
#define MPQC_RANGE_OPERATOR_HPP

#include "mpqc/range.hpp"
#include <boost/preprocessor/repetition.hpp>

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


    template<class _Object, class _Type, class _Const>
    struct range_operator_base {
    protected:
        _Object object_;
        explicit range_operator_base(_Object object) : object_(object) {}

    public:

#define MPQC_RANGE_CONST_OPERATOR(Z, N, TYPE)                                   \
        template<BOOST_PP_ENUM_PARAMS(N, class _R)>                             \
        TYPE operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, const _R, &r)) const {   \
            return object_->operator()(rangify(BOOST_PP_ENUM_PARAMS(N, r)));    \
        }

#define MPQC_RANGE_OPERATOR(Z, N, TYPE)                                         \
        template<BOOST_PP_ENUM_PARAMS(N, class _R)>                             \
        TYPE operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, const _R, &r)) {         \
            return object_->operator()(rangify(BOOST_PP_ENUM_PARAMS(N, r)));    \
        }

        BOOST_PP_REPEAT_FROM_TO(1, 5, MPQC_RANGE_OPERATOR, _Type)
        BOOST_PP_REPEAT_FROM_TO(1, 5, MPQC_RANGE_CONST_OPERATOR, _Const)

#undef MPQC_RANGE_CONST_OPERATOR
#undef MPQC_RANGE_OPERATOR

    };

} // namespace mpqc

#endif // MPQC_RANGE_OPERATOR_HPP
