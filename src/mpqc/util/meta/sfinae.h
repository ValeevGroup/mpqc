
#ifndef MPQC4_SRC_MPQC_UTIL_META_SFINAE_H_
#define MPQC4_SRC_MPQC_UTIL_META_SFINAE_H_

#include <type_traits>

namespace mpqc {
namespace meta {

/// true if decayed T is Base, or is derived from it
template <typename Base, typename T>
using disable_if_same_or_derived = typename std::enable_if<
    !std::is_base_of<Base, typename std::decay<T>::type>::value>::type;

}
}

#endif  // MPQC4_SRC_MPQC_UTIL_META_SFINAE_H_
