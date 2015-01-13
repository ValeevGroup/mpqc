#pragma once
#ifndef TCC_INTEGRALS_TACOMPUTEFUNCTORS_H
#define TCC_INTEGRALS_TACOMPUTEFUNCTORS_H

#include "../include/libint.h"
#include "../include/btas.h"
#include "../include/tiledarray.h"
#include <chrono>

#include "../molecule/cluster_collapse.h"

#include "../basis/cluster_shells.h"
#include "../basis/basis.h"

#include "integral_engine_pool.h"

#include "../tensor/tile_pimpl_devel.h"

#include "../utility/get_type.h"
#include "../utility/expand_container.h"


namespace tcc {
namespace integrals {
namespace compute_functors {
namespace ta_detail {

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
unpack_shells(std::array<std::vector<libint2::Shell>, N> const &shell_vec,
              std::array<std::size_t, N> const &idx, libint2::Shell const &) {
    std::array<libint2::Shell const *, EngineOrder> shells;
    for (auto i = 0ul; i < N; ++i) {
        const auto index = idx[i];
        shells[i] = &(shell_vec[i][index]);
    }
    return shells;
}

template <>
std::array<libint2::Shell const *, 4> unpack_shells<4, 2>(
    std::array<std::vector<libint2::Shell>, 2> const &shell_vec,
    std::array<std::size_t, 2> const &idx, libint2::Shell const &unit) {

    std::array<libint2::Shell const *, 4> shells;
    shells[0] = &shell_vec[0][idx[0]];
    shells[1] = &unit;
    shells[2] = &shell_vec[1][idx[1]];
    shells[3] = &unit;
    return shells;
}

template <>
std::array<libint2::Shell const *, 4> unpack_shells<4, 3>(
    std::array<std::vector<libint2::Shell>, 3> const &shell_vec,
    std::array<std::size_t, 3> const &idx, libint2::Shell const &unit) {

    std::array<libint2::Shell const *, 4> shells;
    shells[0] = &shell_vec[0][idx[0]];
    shells[1] = &unit;
    shells[2] = &shell_vec[1][idx[1]];
    shells[3] = &shell_vec[2][idx[2]];
    return shells;
}

template <typename Engine, typename View, typename... Shells>
void ta_btas_compute_kernel(Engine &engine, View &view, Shells &&... shells) {

    auto buf = engine.compute(std::forward<Shells>(shells)...);
    // const auto size = view.size();
    // std::copy(buf, buf + size, view.begin());

    auto i = 0;
    const auto end = view.end();
    for (auto it = view.begin(); it != end; ++it) {
        *it = buf[i++];
    }
}

template <typename Engine, typename View, std::size_t N, std::size_t... Is>
void wrap_compute_kernel(Engine &engine, View &view,
                         std::array<libint2::Shell const *, N> const &shell_array,
                         utility::seq<Is...>) {
    ta_btas_compute_kernel(engine, view, *std::get<Is>(shell_array)...);
}

template <typename T, unsigned long>
struct ta_integral_wrapper;

template <typename T>
struct ta_integral_wrapper<T, 2ul> {
    template <typename Engine>
    btas::Tensor<T>
    operator()(Engine engine,
               std::vector<basis::ClusterShells> const &clusters) const {

        btas::Tensor<T> tensor{clusters[0].flattened_nfunctions(),
                               clusters[1].flattened_nfunctions()};

        if (engine.integral_type() == libint2::OneBodyEngine::nuclear) {
            std::vector<std::pair<double, std::array<double, 3>>> q;
            for (auto const &shell_cluster : clusters) {
                for (auto const &atom :
                     molecule::collapse_to_atoms(shell_cluster.cluster())) {
                    auto c = atom.center();
                    q.push_back({static_cast<double>(atom.charge()),
                                 {{c[0], c[1], c[2]}}});
                }
            }
            engine.set_q(q);
        }

        auto sh1_start = 0;
        for (auto const &sh1 : clusters[0].flattened_shells()) {
            auto sh1_size = sh1.size();
            auto sh2_start = 0;
            for (auto const &sh2 : clusters[1].flattened_shells()) {
                auto sh2_size = sh2.size();

                auto lower_bound = {sh1_start, sh2_start};
                auto upper_bound = {sh1_start + sh1_size, sh2_start + sh2_size};

                auto view = btas::make_view(
                    tensor.range().slice(lower_bound, upper_bound),
                    tensor.storage());

                ta_btas_compute_kernel(engine, view, sh1, sh2);

                sh2_start += sh2_size;
            }
            sh1_start += sh1_size;
        }

        return tensor;
    }

}; // integralwrapper<T, 2ul> specialization

template <typename T>
struct ta_integral_wrapper<T, 3ul> {

    template <std::size_t N>
    using RangeNd = btas::RangeNd<CblasRowMajor, std::array<long, N>>;

    template <typename Engine>
    btas::Tensor<T, RangeNd<3>>
    operator()(Engine engine,
               std::array<basis::ClusterShells, 3> const &clusters) const {

        btas::Tensor<T, RangeNd<3>> tensor{clusters[0].flattened_nfunctions(),
                                           clusters[1].flattened_nfunctions(),
                                           clusters[2].flattened_nfunctions()};

        libint2::Shell unit{
            {0.0}, // exponent
            {{0, false, {1.0}}},
            {{0.0, 0.0, 0.0}} // placed at origin
        };
        unit.renorm();

        const auto shell_set1 = clusters[0].flattened_shells();
        const auto shell_set2 = clusters[1].flattened_shells();
        const auto shell_set3 = clusters[2].flattened_shells();

        auto sh1_start = 0;
        for (auto const &sh1 : shell_set1) {
            auto sh1_size = sh1.size();
            const auto sh1_upper = sh1_start + sh1_size;

            auto sh2_start = 0;
            for (auto const &sh2 : shell_set2) {
                auto sh2_size = sh2.size();
                const auto sh2_upper = sh2_start + sh2_size;

                auto sh3_start = 0;
                for (auto const &sh3 : shell_set3) {
                    auto sh3_size = sh3.size();
                    const auto sh3_upper = sh3_start + sh3_size;

                    auto lower_bound = {sh1_start, sh2_start, sh3_start};
                    auto upper_bound = {sh1_upper, sh2_upper, sh3_upper};

                    auto view = btas::make_view(
                        tensor.range().slice(lower_bound, upper_bound),
                        tensor.storage());

                    ta_btas_compute_kernel(engine, view, sh1, unit, sh2, sh3);

                    sh3_start = sh3_upper;
                }
                sh2_start = sh2_upper;
            }
            sh1_start = sh1_upper;
        }

        return tensor;
    }

}; // integralwrapper<T, 2ul> specialization

template <typename T>
struct new_ta_integral_wrapper {

    template <std::size_t N>
    using RangeNd = btas::RangeNd<CblasRowMajor, std::array<long, N>>;

    template <typename Engine, std::size_t N>
    btas::Tensor<T, RangeNd<N>>
    operator()(Engine engine,
               std::array<basis::ClusterShells, N> const &clusters) const {

        auto nshell_combos = 1ul;
        std::array<std::vector<libint2::Shell>, N> shells;
        std::array<std::size_t, N> shell_extent;
        for (auto i = 0; i < N; ++i) {
            shells[i] = clusters[i].flattened_shells();
            const auto cluster_size = shells[i].size();
            shell_extent[i] = cluster_size;
            nshell_combos *= cluster_size;
        }

        std::array<std::size_t, N> tensor_extent;
        for (auto i = 0ul; i < N; ++i) {
            tensor_extent[i] = clusters[i].flattened_nfunctions();
        }

        btas::Tensor<T, RangeNd<N>> tensor{btas::Range(tensor_extent)};

        std::array<std::size_t, N> lower_bound;
        lower_bound.fill(0);
        std::array<std::size_t, N> upper_bound;
        for (auto i = 0ul; i < N; ++i) {
            upper_bound[i] = shells[i][0].size();
        }

        libint2::Shell unit{
            {0.0}, // exponent
            {{0, false, {1.0}}},
            {{0.0, 0.0, 0.0}} // placed at origin
        };
        unit.renorm();

        std::array<std::size_t, N> old_idx;
        old_idx.fill(0);
        for (auto i = 0ul; i < nshell_combos; ++i) {

            const auto idx = ord_to_idx(i, shell_extent);

            for (auto i = 0; i < N; ++i) {
                const auto idx_i = idx[i];
                if (old_idx[i] != idx_i) {
                    const auto shell_size = shells[i][idx_i].size();
                    if (idx_i != 0) {
                        lower_bound[i] = upper_bound[i];
                        upper_bound[i] += shell_size;
                    } else {
                        lower_bound[i] = 0;
                        upper_bound[i] = shell_size;
                    }
                    old_idx[i] = idx_i;
                }
            }

            auto shell_set
                = unpack_shells<integrals::pool_order<decltype(engine)>()>(
                    shells, idx, unit);

            auto view = btas::make_view(
                tensor.range().slice(lower_bound, upper_bound),
                tensor.storage());

            wrap_compute_kernel(engine, view, shell_set,
                                utility::gen_seq<shell_set.size()>{});
        }

        return tensor;
    }
};


template <typename T>
struct ta_integral_wrapper<T, 4ul> {

    template <typename Engine>
    btas::Tensor<T>
    operator()(Engine engine,
               std::vector<basis::ClusterShells> const &clusters) const {

        btas::Tensor<T> tensor{clusters[0].flattened_nfunctions(),
                               clusters[1].flattened_nfunctions(),
                               clusters[2].flattened_nfunctions(),
                               clusters[3].flattened_nfunctions()};

        auto sh1_start = 0;
        for (auto const &sh1 : clusters[0].flattened_shells()) {
            auto sh1_size = sh1.size();
            auto sh2_start = 0;

            for (auto const &sh2 : clusters[1].flattened_shells()) {
                auto sh2_size = sh2.size();
                auto sh3_start = 0;
                for (auto const &sh3 : clusters[2].flattened_shells()) {
                    auto sh3_size = sh3.size();
                    auto sh4_start = 0;

                    for (auto const &sh4 : clusters[3].flattened_shells()) {
                        auto sh4_size = sh4.size();

                        auto lower_bound
                            = {sh1_start, sh2_start, sh3_start, sh4_start};
                        auto upper_bound
                            = {sh1_start + sh1_size, sh2_start + sh2_size,
                               sh3_start + sh3_size, sh4_start + sh4_size};

                        auto view = btas::make_view(
                            tensor.range().slice(lower_bound, upper_bound),
                            tensor.storage());

                        ta_btas_compute_kernel(engine, view, sh1, sh2, sh3,
                                               sh4);

                        sh4_start += sh4_size;
                    }
                    sh3_start += sh3_size;
                }
                sh2_start += sh2_size;
            }
            sh1_start += sh1_size;
        }

        return tensor;
    }

}; // integralwrapper<T, 2ul> specialization
} // namespace ta_detail

template <typename Index, std::size_t N>
std::array<basis::ClusterShells, N>
get_shells(Index idx,
           std::array<std::shared_ptr<std::vector<basis::ClusterShells>>,
                      N> const &shell_ptrs) {
    std::array<basis::ClusterShells, N> clusters;

    for (auto i = 0ul; i < N; ++i) {
        const auto cluster_number = idx[i];
        auto shell_ptr = shell_ptrs[i];
        clusters[i] = ((*shell_ptr)[cluster_number]);
    }

    return clusters;
}

template <typename T>
struct TaTileFunctor {
    using TileType = TiledArray::Tensor<T>;

    template <typename Index, typename SharedEnginePool, std::size_t N>
    TileType
    operator()(TiledArray::Range range, Index index, SharedEnginePool engines,
               std::array<std::shared_ptr<std::vector<basis::ClusterShells>>,
                          N> const &shell_ptrs) const {

        // auto btas_tensor = ta_detail::ta_integral_wrapper<T, N>{}(
        //    engines->local(), get_shells(index, shell_ptrs));
        auto btas_tensor = ta_detail::new_ta_integral_wrapper<T>{}(
            engines->local(), get_shells(index, shell_ptrs));

        // copy to TA
        TileType ta_tensor{range};
        auto data = ta_tensor.data();
        auto const size = ta_tensor.size();
        for (auto i = 0ul; i < size; ++i) {
            *(data + i) = btas_tensor.storage()[i];
        }

        return ta_tensor;
    }

}; // TaTileFunctor


} // namespace compute_functors
} // namespace integrals
} // namespace tcc

#endif /* end of include guard: TCC_INTEGRALS_TACOMPUTEFUNCTORS_H */
