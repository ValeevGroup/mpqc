#pragma once
#ifndef TCC_INTEGRAL_TILEENGINEHELPER_H
#define TCC_INTEGRAL_TILEENGINEHELPER_H

#include "../include/libint.h"
#include "../include/btas.h"

#include "../utility/meta/gen_seq.h"

// #include "../basis/cluster_shells.h"
#include "../common/typedefs.h"
#include "../basis/basis.h"
#include "../basis/shell_vec_functions.h"

#include "integral_engine_pool.h"

#include <array>

namespace tcc {
namespace integrals {
namespace detail {

template <std::size_t N>
using RangeNd = btas::RangeNd<CblasRowMajor, std::array<long, N>>;

template <std::size_t N>
using TileType = btas::Tensor<double, RangeNd<N>>;

template <std::size_t N>
std::array<std::size_t, N>
ord_to_idx(std::size_t ord, std::array<std::size_t, N> const &extent) {

    std::array<std::size_t, N> idx;

    auto dim_size = 1ul;
    for (auto i = 1ul; i < N; ++i) {
        dim_size *= extent[i];
    }

    for (auto i = 0ul; i < N; ++i) {
        idx[i] = ord / dim_size;
        ord = ord % dim_size;
        if (N - 1 != i) dim_size /= extent[i + 1];
    }

    return idx;
}

template <std::size_t EngineOrder, std::size_t N>
std::array<libint2::Shell const *, EngineOrder>
unpack_shells(std::array<std::vector<Shell>, N> const &shell_vecs,
              std::array<std::size_t, N> const &idx, Shell const &) {
    std::array<Shell const *, EngineOrder> shell_set;
    for (auto i = 0ul; i < N; ++i) {
        const auto index = idx[i];
        shell_set[i] = &(shell_vecs[i][index]);
    }
    return shell_set;
}

template <>
std::array<libint2::Shell const *, 4> unpack_shells<4, 2>(
    std::array<std::vector<libint2::Shell>, 2> const &shell_vecs,
    std::array<std::size_t, 2> const &idx, libint2::Shell const &unit) {

    std::array<libint2::Shell const *, 4> shell_set;
    shell_set[0] = &shell_vecs[0][idx[0]];
    shell_set[1] = &unit;
    shell_set[2] = &shell_vecs[1][idx[1]];
    shell_set[3] = &unit;
    return shell_set;
}

template <>
std::array<libint2::Shell const *, 4> unpack_shells<4, 3>(
    std::array<std::vector<libint2::Shell>, 3> const &shell_vecs,
    std::array<std::size_t, 3> const &idx, libint2::Shell const &unit) {

    std::array<libint2::Shell const *, 4> shell_set;
    shell_set[0] = &shell_vecs[0][idx[0]];
    shell_set[1] = &unit;
    shell_set[2] = &shell_vecs[1][idx[1]];
    shell_set[3] = &shell_vecs[2][idx[2]];
    return shell_set;
}


template <typename Engine, typename View, typename... Shells>
inline void compute_kernel(Engine &engine, View &view, Shells &&... shells) {
    auto buf = engine.compute(std::forward<Shells>(shells)...);
    const auto size = view.size();
    std::copy(buf, buf + size, view.begin());
}

template <typename Engine, typename View, std::size_t N, std::size_t... Is>
void unpack_to_compute_kernel(Engine &engine, View &view,
                              std::array<Shell const *, N> const &shell_ptrs,
                              utility::seq<Is...>) {
    compute_kernel(engine, view, *std::get<Is>(shell_ptrs)...);
}

template <std::size_t N>
std::array<std::vector<Shell>, N>
flattened_shells(std::array<ShellVec, N> const &clustered_shells) {
    std::array<std::vector<Shell>, N> shells;
    for (auto i = 0ul; i < N; ++i) {
        shells[i] = clustered_shells[i];
    }
    return shells;
}

template <std::size_t N>
std::array<std::size_t, N>
get_shell_extent(std::array<std::vector<Shell>, N> const &shells) {
    std::array<std::size_t, N> extent;
    for (auto i = 0ul; i < N; ++i) {
        extent[i] = shells[i].size();
    }
    return extent;
}

template <std::size_t N>
void update_bounds(std::array<std::size_t, N> const &idx,
                   std::array<std::size_t, N> &old_idx,
                   std::array<std::size_t, N> &lower_bound,
                   std::array<std::size_t, N> &upper_bound,
                   std::array<std::vector<Shell>, N> const &shells) {

    for (auto i = 0ul; i < N; ++i) {

        const auto old_i = old_idx[i];
        const auto idx_i = idx[i];

        // Only update bounds for index which has changed.
        if (old_i != idx_i) {
            const auto shell_size = shells[i][idx_i].size();
            if (idx_i != 0) { // index is still growing
                lower_bound[i] = upper_bound[i];
                upper_bound[i] += shell_size;
            } else { // index has rolled over back to 0
                lower_bound[i] = 0;
                upper_bound[i] = shell_size;
            }

            // update the old value with the current value.
            old_idx[i] = idx_i;
        }
    }
}


template <typename Engine, std::size_t N>
TileType<N>
compute_tile(Engine engine,
             std::array<ShellVec, N> const &clusters) {

    // Compute number of functions in each dimension
    std::array<std::size_t, N> tensor_extent;
    for (auto i = 0ul; i < N; ++i) {
        tensor_extent[i] = mpqc::basis::nfunctions(clusters[i]);
    }

    // Flatten the shells, get the number of shells in each dim, and compute
    // number of ordinals.
    const auto shells = flattened_shells(clusters);
    const auto shell_extent = get_shell_extent(shells);
    const auto nshell_ordinals = std::accumulate(
        shell_extent.begin(), shell_extent.end(), 1ul,
        [](std::size_t result, std::size_t val) { return result * val; });

    TileType<N> tensor{btas::Range(tensor_extent)};

    // initialize arrays to hold the upper and lower bounds for each shell set.
    std::array<std::size_t, N> lower_bound;
    lower_bound.fill(0);
    std::array<std::size_t, N> upper_bound;
    for (auto i = 0ul; i < N; ++i) {
        upper_bound[i] = shells[i][0].size();
    }

    // array to keep the old idx values
    std::array<std::size_t, N> idx;
    idx.fill(0);
    std::array<std::size_t, N> old_idx;
    old_idx.fill(0);

    Shell unit_ = Shell::unit();

    for (auto ord = 0ul; ord < nshell_ordinals; ++ord) {
        idx = ord_to_idx(ord, shell_extent);

        update_bounds(idx, old_idx, lower_bound, upper_bound, shells);

        constexpr std::size_t engine_dim = pool_order<decltype(engine)>();

        auto shell_set
            = unpack_shells<engine_dim>(shells, idx, unit_);

        auto view = btas::make_view(
            tensor.range().slice(lower_bound, upper_bound), tensor.storage());

        unpack_to_compute_kernel(engine, view, shell_set,
                            utility::gen_seq<shell_set.size()>{});
    }

    return tensor;
}


} // namespace detail
} // namespace integrals
} // namespace tcc


#endif /* end of include guard: TCC_INTEGRAL_TILEENGINEHELPER_H */
