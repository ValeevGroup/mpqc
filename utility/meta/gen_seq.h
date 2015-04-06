#pragma once
#ifndef TCC_UTLITIY_EXPANDCONTAINER_H
#define TCC_UTLITIY_EXPANDCONTAINER_H

#include <utility>
#include <tuple>
#include <array>

namespace tcc {
namespace utility {

/// seq and gen_seq are taken from stackoverflow's pretty print tuple answer.
template <std::size_t...>
struct seq {};

template <std::size_t N, std::size_t... Is>
struct gen_seq : gen_seq<N - 1, N - 1, Is...> {};

template <std::size_t... Is>
struct gen_seq<0, Is...> : seq<Is...> {};

} // namespace utility
} // namespace tcc

#endif /* end of include guard: TCC_UTLITIY_EXPANDCONTAINER_H */
