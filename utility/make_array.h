#pragma once
#ifndef TCC_UTILITY_MAKEARRAY_H
#define TCC_UTILITY_MAKEARRAY_H

#include <array>
#include "get_type.h"

namespace tcc {
namespace utility {

template <typename... Args>
std::array<first_type_t<Args...>, sizeof...(Args)> make_array(Args... args) {
    return std::array<first_type_t<Args...>, sizeof...(Args)>{{args...}};
}

} // namespace utility
} // namespace tcc

#endif /* end of include guard: TCC_UTILITY_MAKEARRAY_H */
