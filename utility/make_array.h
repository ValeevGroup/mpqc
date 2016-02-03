#pragma once
#ifndef MPQC_UTILITY_MAKEARRAY_H
#define MPQC_UTILITY_MAKEARRAY_H

#include <array>
#include "meta/get_type.h"

namespace mpqc {
namespace utility {

template <typename... Args>
std::array<meta::first_type_t<Args...>, sizeof...(Args)> make_array(Args... args) {
    return std::array<meta::first_type_t<Args...>, sizeof...(Args)>{{args...}};
}

} // namespace utility
} // namespace mpqc

#endif /* end of include guard: MPQC_UTILITY_MAKEARRAY_H */
