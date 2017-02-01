
#ifndef MPQC4_SRC_MPQC_UTIL_META_PREDICATES_H_
#define MPQC4_SRC_MPQC_UTIL_META_PREDICATES_H_

#include "mpqc/util/meta/get_type.h"

namespace mpqc {
namespace utility {
namespace meta {

template <typename... Args>
struct is_homogeneous_parameter_pack;

template <typename FirstArg>
struct is_homogeneous_parameter_pack<FirstArg> : std::true_type {};

template <typename FirstArg, typename... RestOfArgs>
struct is_homogeneous_parameter_pack<FirstArg, RestOfArgs...>
    : std::is_same<FirstArg, first_type_t<RestOfArgs...>>::type {};

#if __cplusplus >= 201402L
// C++14 only
template <typename... Args>
bool is_homogeneous_parameter_pack_v =
    is_homogeneous_parameter_pack<Args...>::value;
#endif
}
}
}

#endif  // MPQC4_SRC_MPQC_UTIL_META_PREDICATES_H_
