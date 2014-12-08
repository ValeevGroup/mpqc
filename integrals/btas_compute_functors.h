#pragma once
#ifndef TCC_INTEGRALS_BTASCOMPUTEFUNCTORS_H
#define TCC_INTEGRALS_BTASCOMPUTEFUNCTORS_H

#include "../include/libint.h"
#include "../include/btas.h"
#include "../basis/cluster_shells.h"
#include "../include/tiledarray.h"
#include "../basis/basis.h"
#include "integral_engine_pool.h"

#include "../tensor/tile_pimpl_devel.h"


namespace tcc {
namespace integrals {
namespace compute_functors {
namespace detail {

template <typename Engine, typename T, typename... Shells>
void btas_compute_kernel(Engine &engine, btas::TensorView<T> &view,
                         Shells... shells) {

    auto buf = engine.compute(shells...);

    auto i = 0;
    auto end = view.end();
    for (auto it = view.begin(); it != end; ++it) {
        *it = buf[i++];
        //        *it = std::max({shells.O[2]...});
    }
}

template <typename T, unsigned long>
struct integral_wrapper;

template <typename T>
struct integral_wrapper<T, 2ul> {
    template <typename Engine>
    btas::Tensor<T>
    operator()(Engine engine,
               std::vector<basis::ClusterShells> const &clusters) const {

        btas::Tensor<T> tensor{clusters[0].flattened_nfunctions(),
                               clusters[1].flattened_nfunctions()};

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

                btas_compute_kernel(engine, view, sh1, sh2);

                sh2_start += sh2_size;
            }
            sh1_start += sh1_size;
        }

        return tensor;
    }

}; // integralwrapper<T, 2ul> specialization

template <typename T>
struct integral_wrapper<T, 4ul> {
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

            for (auto const &sh2 : clusters[0].flattened_shells()) {
                auto sh2_size = sh2.size();
                auto sh3_start = 0;
                for (auto const &sh3 : clusters[0].flattened_shells()) {
                    auto sh3_size = sh3.size();
                    auto sh4_start = 0;

                    for (auto const &sh4 : clusters[1].flattened_shells()) {
                        auto sh4_size = sh4.size();

                        auto lower_bound
                            = {sh1_start, sh2_start, sh3_start, sh4_start};
                        auto upper_bound
                            = {sh1_start + sh1_size, sh2_start + sh2_size,
                               sh3_start + sh3_size, sh4_start + sh4_size};

                        auto view = btas::make_view(
                            tensor.range().slice(lower_bound, upper_bound),
                            tensor.storage());

                        btas_compute_kernel(engine, view, sh1, sh2, sh3, sh4);

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
} // namespace detail


template <typename T>
struct BtasTileFunctor {
    using TileType = tensor::TilePimplDevel<T>;

    template <typename It, typename SharedEnginePool>
    TileType operator()(It it, basis::Basis const *basis,
                        SharedEnginePool engines) const {

        auto idx = it.index();
        std::vector<basis::ClusterShells> clusters;
        clusters.reserve(idx.size());

        for (auto const &i : idx) {
            clusters.push_back(basis->cluster_shells()[i]);
        }

        auto btas_tensor
            = detail::integral_wrapper<T, pool_order<SharedEnginePool>()>{}(
                engines->local(), clusters);

        auto range = it.make_range();
        return TileType{
            std::move(range), tensor::TileVariantDevel<double>{std::move(btas_tensor)}};
    }

}; // BtasTileFunctor


} // namespace compute_functors
} // namespace integrals
} // namespace tcc

#endif /* end of include guard: TCC_INTEGRALS_BTASCOMPUTEFUNCTORS_H */
