#pragma once
#ifndef TCC_UTILITY_GETTYPE_H
#define TCC_UTILITY_GETTYPE_H

namespace tcc {
namespace utility {

template <typename first, typename... Rest>
struct first_type {
    using type = first;
};

template <typename... Args>
using first_type_t = typename first_type<Args...>::type;

} // namespace utility
} // namespace tcc


#endif /* end of include guard: TCC_UTILITY_GETTYPE_H */
