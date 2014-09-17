// Modified from stackoverflow question
// http://stackoverflow.com/questions/17032310/how-to-make-a-variadic-is-same

#ifndef TCC_INCLUDE_TCC_ALL_SAME_H
#define TCC_INCLUDE_TCC_ALL_SAME_H

#include <type_traits>

namespace tcc {

template <typename T, typename... Rest>
struct all_same : std::false_type {};

template <typename T, typename First>
struct all_same<T, First> : std::is_same<T, First> {};

template <typename T, typename First, typename... Rest>
struct all_same<T, First, Rest...> : std::integral_constant
                                     <bool, std::is_same<T, First>::value
                                      &&all_same<First, Rest...>::value> {};
}
#endif // TCC_INCLUDE_TCC_ALL_SAME_H
