#pragma once
#ifndef TCC_UTILITY_META_STOREVARIADICPACK_H
#define TCC_UTILITY_META_STOREVARIADICPACK_H

#include <tuple>

namespace mpqc {
namespace utility {
namespace meta {

namespace detail {

struct any {
    template <typename T>
    any(T &&) {}
};
template <typename T, typename U>
struct fst {
    typedef T type;
};
template <typename... Ts>
struct select_impl {
    template <typename U>
    static U &&select(typename fst<any, Ts>::type..., U &&u) {
        return static_cast<U &&>(u);
    }
};

} // namespace detail

template <typename T, typename... Ts>
auto back(T &&t, Ts &&... ts) -> decltype(detail::select_impl<Ts...>::select(
      static_cast<T &&>(t), static_cast<Ts &&>(ts)...)) {
    return detail::select_impl<Ts...>::select(static_cast<T &&>(t),
                                              static_cast<Ts &&>(ts)...);
}

/// How to get the last type out of a parameter pack
/// adapted taken from http://stackoverflow.com/questions/4691657/
///      is-it-possible-to-store-a-template-parameter-pack-without-expanding-it
template <typename... Args>
struct last_type {
    static constexpr auto size = sizeof...(Args);
    using type
          = decltype(std::get<size - 1>(std::declval<std::tuple<Args...>>()));
};

// I need this specialization for cases when the pack is empty. It shouldn't
// return void, but for now I am using it to avoid extra work in tcc_tile
template <>
struct last_type<> {
    using type = void;
};

} // namespace meta
} // namespace utility
} // namespace mpqc


#endif // TCC_UTILITY_META_STOREVARIADICPACK_H
