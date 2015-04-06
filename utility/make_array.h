#pragma once
#ifndef TCC_UTILITY_MAKEARRAY_H
#define TCC_UTILITY_MAKEARRAY_H

#include <array>
#include "meta/get_type.h"

namespace tcc {
namespace utility {

template <typename... Args>
std::array<meta::first_type_t<Args...>, sizeof...(Args)> make_array(Args... args) {
    return std::array<meta::first_type_t<Args...>, sizeof...(Args)>{{args...}};
}

} // namespace utility
} // namespace tcc

#endif /* end of include guard: TCC_UTILITY_MAKEARRAY_H */
