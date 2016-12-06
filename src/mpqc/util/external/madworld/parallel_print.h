#ifndef MPQC4_SRC_MPQC_UTIL_EXTERNAL_MADWORLD_PARALLEL_PRINT_H_
#define MPQC4_SRC_MPQC_UTIL_EXTERNAL_MADWORLD_PARALLEL_PRINT_H_

#include <fstream>
#include <iostream>
#include <string>
#include <tuple>

#include <madness/world/world.h>

#include "mpqc/util/misc/exenv.h"

namespace mpqc {
namespace utility {
namespace aux {

// Pretty print tuple from
// http://stackoverflow.com/questions/6245735/pretty-print-stdtuple
template <std::size_t...>
struct seq {};

template <std::size_t N, std::size_t... Is>
struct gen_seq : gen_seq<N - 1, N - 1, Is...> {};

template <std::size_t... Is>
struct gen_seq<0, Is...> : seq<Is...> {};

template <class Ch, class Tr, class Tuple, std::size_t... Is>
void print_tuple(std::basic_ostream<Ch, Tr> &os, Tuple const &t, seq<Is...>) {
  using swallow = int[];
  (void)swallow{0, (void(os << std::get<Is>(t)), 0)...};
}
}  // aux::

template <class Ch, class Tr, class... Args>
auto operator<<(std::basic_ostream<Ch, Tr> &os, std::tuple<Args...> const &t)
    -> std::basic_ostream<Ch, Tr> & {
  aux::print_tuple(os, t, aux::gen_seq<sizeof...(Args)>());
  return os;
}

template <typename... Args>
void print_par(madness::World &world, Args &&... args) {
  auto t = std::make_tuple<Args...>(std::forward<Args>(args)...);
  ExEnv::out0() << t << std::flush;
}

inline void print_file(madness::World &world, const std::string &file) {
  if (world.rank() == 0) {
    std::ifstream file_stream(file);

    std::string line;
    while (getline(file_stream, line)) {
      std::cout << '\t' << line << std::endl;
    }

    file_stream.close();
  }
}

}  // namespace utility
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_UTIL_EXTERNAL_MADWORLD_PARALLEL_PRINT_H_
