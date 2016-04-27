/* Collapse to base
 *
 * A set of functions to help Clusterable
 * 2016 Drew Lewis
 */
#pragma once
#ifndef MPQC_CLUSTERING_COLLAPSETOBASE_H
#define MPQC_CLUSTERING_COLLAPSETOBASE_H

#include <vector>
#include <type_traits>

#include "mpqc/util/misc/type_traits.h"

namespace mpqc {
namespace clustering {

template <typename, class = void_t<>>
struct has_function_collapse_to_base : std::false_type {};

template <typename T>
struct has_function_collapse_to_base<
    T, void_t<decltype(std::declval<T>().collapse_to_base())>>
    : std::true_type {};

template <typename Base>
struct collapse_to_base {
  template <typename T, enable_if_t<std::is_same<T, Base>::value> * = nullptr>
  std::vector<Base> operator()(Base const &b) {
    return {b};
  }

  template <
      typename T,
      enable_if_t<not std::is_same<T, Base>::value &&
                  not has_function_collapse_to_base<T>::value> * = nullptr>
  std::vector<Base> operator()(T const &t) {
    std::vector<Base> elems;
    for (auto const &elem : t) {
      using elem_type = typename std::remove_reference<decltype(elem)>::type;
      using nonconst_type = typename std::remove_const<elem_type>::type;
      auto temp_elems = this->operator()<nonconst_type>(elem);
      elems.insert(elems.end(), temp_elems.begin(), temp_elems.end());
    }

    return elems;
  }

  template <typename T,
            enable_if_t<not std::is_same<T, Base>::value &&
                        has_function_collapse_to_base<T>::value> * = nullptr>
  std::vector<Base> operator()(T const &t) {
    return t.collapse_to_base();
  }
};

}  // namespace clustering
}  // namespace mpqc
#endif  //  MPQC_CLUSTERING_COLLAPSETOBASE_H
