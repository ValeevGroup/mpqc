#pragma once
#ifndef MPQC_UTLITY_PARALLELPRINT_H
#define MPQC_UTLITY_PARALLELPRINT_H

#include <tuple>
#include <string>
#include <iostream>
#include <fstream>
#include "../include/tiledarray.h"
#include <mpqc/chemistry/qc/basis/basis.h>
#include "cc_utility.h"

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
    (void)swallow{0,
                  (void(os << std::get<Is>(t)), 0)...};
}
} // aux::

template <class Ch, class Tr, class... Args>
auto operator<<(std::basic_ostream<Ch, Tr> &os, std::tuple<Args...> const &t)
    -> std::basic_ostream<Ch, Tr> & {
    aux::print_tuple(os, t, aux::gen_seq<sizeof...(Args)>());
    return os; 
}


template <typename... Args>
void print_par(madness::World &world, Args&&... args) {
    auto t = std::make_tuple<Args...>(std::forward<Args>(args)...);
    if (world.rank() == 0) {
        std::cout << t << std::flush;
    }
}

template <typename... Args>
void wprint_par(madness::World &world, Args&&... args) {
    auto t = std::make_tuple<Args...>(std::forward<Args>(args)...);
    if (world.rank() == 0) {
        std::wcout << t << std::flush;
    }
}

inline void print_file(madness::World &world, const std::string& file){

    if (world.rank() == 0){
        std::ifstream file_stream(file);

        std::string line;
        while(getline(file_stream,line)){
            std::cout << '\t' << line << std::endl;
        }

        file_stream.close();
    }

}

inline void parallel_print_range_info(madness::World &world, const TA::TiledRange1 &bs_range, const std::string &name){

    if(world.rank() == 0){
        auto minmax_block = cc::minmax_blocksize(bs_range);
        auto average_block = cc::average_blocksize(bs_range);
        std::cout <<  name << " Trange " << std::endl;
        std::cout << bs_range << std::endl;
        std::cout << "Min and Max block size: " << minmax_block.first << " " << minmax_block.second << std::endl;
        std::cout << "Average: " << average_block << std::endl;
        std::cout << std::endl;
    }
}


} // namespace utility
} // namespace mpqc

#endif // MPQC_UTLITY_PARALLELPRINT_H
