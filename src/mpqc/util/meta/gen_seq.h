
#ifndef MPQC4_SRC_MPQC_UTIL_META_GEN_SEQ_H_
#define MPQC4_SRC_MPQC_UTIL_META_GEN_SEQ_H_

#include <array>
#include <tuple>
#include <utility>

namespace mpqc {
namespace utility {

/// seq and gen_seq are taken from stackoverflow's pretty print tuple answer.
template <std::size_t...>
struct seq {};

template <std::size_t N, std::size_t... Is>
struct gen_seq : gen_seq<N - 1, N - 1, Is...> {};

template <std::size_t... Is>
struct gen_seq<0, Is...> : seq<Is...> {};

}  // namespace utility
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_UTIL_META_GEN_SEQ_H_
