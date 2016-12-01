
#ifndef MPQC4_SRC_MPQC_UTIL_META_GET_TYPE_H_
#define MPQC4_SRC_MPQC_UTIL_META_GET_TYPE_H_

#include <type_traits>

namespace mpqc {
namespace utility {
namespace meta {

template <typename first, typename... Rest>
struct first_type {
  using type = first;
};

template <typename... Args>
using first_type_t = typename first_type<Args...>::type;

template <typename... Args>
struct last_type;

template <typename T>
struct last_type<T> {
  using type = T;
};

template <>
struct last_type<> {
  using type = std::false_type;
};

template <typename T, typename... Rest>
struct last_type<T, Rest...> {
  using type = typename last_type<Rest...>::type;
};

// The following code down till back was taken from
// http://stackoverflow.com/questions/13050555/getting-access-to-the-back-of-a-template-parameter-pack-in-o1
// and uses some Eric Niebler magic.
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

}  // namespace detail

// Gets the last element in a variadic parameter pack.
template <typename T, typename... Ts>
auto back(T &&t, Ts &&... ts)
    -> decltype(detail::select_impl<Ts...>::select(static_cast<T &&>(t),
                                                   static_cast<Ts &&>(ts)...)) {
  return detail::select_impl<Ts...>::select(static_cast<T &&>(t),
                                            static_cast<Ts &&>(ts)...);
}

}  // namespace meta
}  // namespace utility
}  // namespace MPQC

#endif  // MPQC4_SRC_MPQC_UTIL_META_GET_TYPE_H_
