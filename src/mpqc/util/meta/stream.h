//
// Created by Eduard Valeyev on 4/25/18.
//

#ifndef MPQC4_SRC_MPQC_UTIL_META_STREAM_H
#define MPQC4_SRC_MPQC_UTIL_META_STREAM_H

#include <madness/world/type_traits.h>

namespace mpqc {
namespace utility {

/// forwards value to std::ostream if possible
template <typename T>
std::enable_if_t<madness::is_ostreammable<T>::value, T>
    to_ostream(T&& value) {
  return std::forward<T>(value);
}

template <typename T>
std::string
to_ostream(T&& value, std::enable_if_t<!madness::is_ostreammable<T>::value>* = nullptr) {
  return "<unknown>";
}

}
}

#endif  // MPQC4_SRC_MPQC_UTIL_META_STREAM_H
