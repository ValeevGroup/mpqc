/*
 * type_traits.h
 *
 *  Created on: Apr 27, 2016
 *      Author: Drew Lewis
 */

#ifndef MPQC_UTIL_MISC_TYPETRAITS_H_
#define MPQC_UTIL_MISC_TYPETRAITS_H_

#include <type_traits>

namespace mpqc {

/*! \brief enable_if_t typedef
 *
 * This is a standard C++ 14 feature, but we don't want to require C++ 14
 * compilers yet
 */
template <bool B, class T = void>
using enable_if_t = typename std::enable_if<B, T>::type;

/*! \brief void_t typedef
 *
 * This is a standard C++ 17 feature, but we don't want to require C++ 17
 * compilers yet
 */
template <typename... Ts>
struct make_void {
  using type = void;
};

template <typename... Ts>
using void_t = typename make_void<Ts...>::type;

}  // namespace mpqc

#endif  // MPQC_UTIL_MISC_TYPETRAITS_H_
