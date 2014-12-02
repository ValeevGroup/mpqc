#pragma once
#ifndef TCC_INTEGRALS_LOWRANKTILEFUNCTORS_H
#define TCC_INTEGRALS_LOWRANKTILEFUNCTORS_H

#include "../include/libint.h"
#include "../include/btas.h"
#include "../basis/cluster_shells.h"
#include "../include/tiledarray.h"
#include "../basis/basis.h"
#include "../basis/cluster_shells.h"
#include "../molecule/cluster_collapse.h"
#include "integral_engine_pool.h"

#include "../tensor/tile_pimpl.h"

namespace tcc {
namespace integrals {
namespace compute_functors {

namespace low_rank_tile_detail {

template <typename T, unsigned long>
struct integral_wrapper;

template <typename T>
struct integral_wrapper<T, 2ul> {
    template <typename Engine>
    tensor::TileVariant<T>
    operator()(Engine engine, std::vector<basis::ClusterShells> const &clusters,
               double cut) const {
        Eigen::MatrixXd tile(clusters[0].flattened_nfunctions(),
                             clusters[1].flattened_nfunctions());

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

                Eigen::MatrixXd block(sh1_size, sh2_size);

                auto buf = engine.compute(sh1, sh2);

                auto size = sh1_size * sh2_size;
                for (auto i = 0ul; i < size; ++i) {
                    *(block.data() + i) = buf[i];
                }

                tile.block(sh1_start, sh2_start, sh1_size, sh2_size) = block;

                sh2_start += sh2_size;
            }
            sh1_start += sh1_size;
        }

        Eigen::MatrixXd L, R;
        if (algebra::lapack::Decompose_Matrix(tile, L, R, cut)) {
            return tensor::TileVariant<T>{
                tensor::FullRankTile<T>{std::move(tile)}};
        } else {
            return tensor::TileVariant<T>{
                tensor::LowRankTile<T>{std::move(L), std::move(R)}};
        }
    }
}; // integral_wrapper<T,2ul> specialilization

} // namespace detail

template <typename T>
struct LowRankTileFunctor {
    double cut_ = 1e-7;
    LowRankTileFunctor() = default;
    LowRankTileFunctor(double cut) : cut_(cut) {}

    using TileType = tensor::TilePimpl<T>;

    template <typename It, typename SharedEnginePool>
    TileType operator()(It it, basis::Basis const *basis,
                        SharedEnginePool engines) const {
        auto idx = it.index();
        std::vector<basis::ClusterShells> clusters;
        clusters.reserve(idx.size());

        for (auto const &i : idx) {
            clusters.push_back(basis->cluster_shells()[i]);
        }

        auto tile = low_rank_tile_detail::
            integral_wrapper<T, pool_order<SharedEnginePool>()>{}(
                engines->local(), clusters, cut_);

        return TileType{it.make_range(), std::move(tile), cut_};
    }

}; // LowRankTileFunctor
}
}
}


#endif /* end of include guard: TCC_INTEGRALS_LOWRANKTILEFUNCTORS_H */
