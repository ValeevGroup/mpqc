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

namespace detail {
struct expand_args_and_apply {
    template <typename Func, typename Tup, std::size_t... Is>
    typename Func::return_type operator()(Func fn, Tup t, seq<Is...>) {
        return fn(std::get<Is>(t)...);
    }
};
} // namespace detail

template <typename Tuple, typename... Args>
struct TupleExpander{

    TupleExpander(Tuple pack, Args &&... leading_args)
        : pack_{pack}, leading_args_{std::forward_as_tuple(leading_args...)} {}

    template <typename Func>
    typename Func::return_type operator()(Func fn){
        auto func_inputs = std::tuple_cat(leading_args_, pack_);
        constexpr auto t_size = std::tuple_size<decltype(func_inputs)>::value;

        return detail::expand_args_and_apply{}(fn, func_inputs,
                                               gen_seq<t_size>{});
    }

    Tuple pack_;
    std::tuple<Args...> leading_args_;
};

template <typename Tuple, typename... Args>
TupleExpander<Tuple, Args...> make_tuple_expander(Tuple &&tup, Args &&... args) {
    return TupleExpander<Tuple, Args...>(std::forward<Tuple>(tup),
                                       std::forward<Args>(args)...);
}

} // namespace utility
} // namespace tcc

#endif /* end of include guard: TCC_UTLITIY_EXPANDCONTAINER_H */
