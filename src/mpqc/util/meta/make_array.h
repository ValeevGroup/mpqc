#pragma once
#ifndef MPQC_UTILITY_MAKEARRAY_H
#define MPQC_UTILITY_MAKEARRAY_H

#include <array>
#include <functional>

#include "mpqc/util/meta/get_type.h"
#include "mpqc/util/meta/predicates.h"

namespace mpqc {
namespace utility {

template <typename... Args>
std::array<meta::first_type_t<Args...>, sizeof...(Args)> make_array(Args... args) {
    return std::array<meta::first_type_t<Args...>, sizeof...(Args)>{{args...}};
}

template <typename... Args>
std::array<std::reference_wrapper<meta::first_type_t<Args...>>, sizeof...(Args)> make_array_of_refs(Args&... args) {
    static_assert(meta::is_homogeneous_parameter_pack<Args...>::value, "inhomogeneous parameter packs not allowed");
    return std::array<std::reference_wrapper<meta::first_type_t<Args...>>, sizeof...(Args)>{{std::ref(args)...}};
}

} // namespace utility
} // namespace mpqc

#endif /* end of include guard: MPQC_UTILITY_MAKEARRAY_H */
