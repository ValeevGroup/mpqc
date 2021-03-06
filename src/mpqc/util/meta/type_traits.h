#ifndef SRC_MPQC_UTIL_EXTERNAL_CPP_TYPE_TRAITS
#define SRC_MPQC_UTIL_EXTERNAL_CPP_TYPE_TRAITS

#include "type_traits.h"

namespace mpqc {
namespace meta {

// C++17 features
#if __cplusplus <= 201402L

// GNU stdlibc++ provides void_t if -gnu++11 or -gnu++14 are given
#if __GNUC__ && defined(__GLIBCXX__) && !__STRICT_ANSI__ && __cplusplus >= 201103L
#define HAVE_VOID_T
#endif

#ifndef HAVE_VOID_T
template <typename... Ts>
struct make_void {
  using type = void;
};
template <typename... Ts>
using void_t = typename make_void<Ts...>::type;
#endif

#endif  // C++17 features

}  // namespace meta
}  // namespace mpqc

#endif  // SRC_MPQC_UTIL_EXTERNAL_CPP_TYPE_TRAITS
