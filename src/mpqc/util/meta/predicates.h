
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

template <typename Class>
struct is_shared_ptr : std::false_type {};

template <typename T>
struct is_shared_ptr<std::shared_ptr<T>> : std::true_type {};

/// checks if T has a ctor that accepts U...
/// @tparam T a type whose ctors are checked
/// @tparam U the arguments to the ctor
template <typename T, typename ... U>
class can_construct {
  /// true case
  template <typename _T, typename ... _U>
  static auto __test(_U ... args) -> decltype(_T(args...), std::true_type());
  /// false case
  template <typename _T, typename ... _U>
  static std::false_type __test(...);

 public:
  static constexpr const bool value = decltype(can_construct<T,U...>::__test<T,U...>(std::declval<U>()...))::value;
};


}
}
}

#endif  // MPQC4_SRC_MPQC_UTIL_META_PREDICATES_H_
