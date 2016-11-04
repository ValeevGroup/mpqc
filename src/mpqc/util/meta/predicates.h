
#ifndef UTILITY_META_PREDICATES_H_
#define UTILITY_META_PREDICATES_H_

namespace mpqc {
namespace utility {
namespace meta {

#include "mpqc/util/meta/get_type.h"

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

#endif /* UTILITY_META_PREDICATES_H_ */
