//
// task_integrals.h
//
// Copyright (C) 2014 Drew Lewis
// Maintainer Drew Lewis
//
// Based on file task_integrals.hpp from mpqc
//

#ifndef TCC_INTEGRALS_SPARSETASKINTEGRALS_H
#define TCC_INTEGRALS_SPARSETASKINTEGRALS_H

#include "../include/tiledarray.h"
#include "../include/btas.h"
#include "../include/libint2.h"

#include "../tensor/tile_algebra.h"

#include "../basis/cluster_shells.h"
#include "integral_engine_pool.h"
#include "btas_compute_functors.h"

#include "../basis/basis.h"

#include "../tensor/tile_pimpl_devel.h"

namespace tcc {
namespace integrals {

using TaTRange1 = TiledArray::TiledRange1;
using TaTRange = TiledArray::TiledRange;

template <typename It, typename SharedEnginePool, typename TileFunctor,
          typename... BasesPtrs>
void call_tf(std::shared_ptr<std::vector<std::array<double, 3>>> tile_vec_ptr,
             SharedEnginePool engines, TaTRange trange, It it, TileFunctor tf,
             BasesPtrs &&... bases_ptrs) {

    auto &vec = *tile_vec_ptr;
    auto const ordinal = *it;
    auto const idx = trange.tiles().idx(ordinal);
    auto tile = tf(trange.make_tile_range(ordinal), idx, engines, bases_ptrs...);

    auto const &extent = tile.range().size();
    auto X = extent[0]; 
    auto kl = extent[1] * extent[2];

    double dense_storage = double(X * kl) * 8.0 / 1.0e9;
    double sparse_storage = dense_storage;
    double low_rank_storage = dense_storage;

    auto norm = tile.norm();
    if (norm < TiledArray::SparseShape<float>::threshold()) {
        sparse_storage = 0.0;
        low_rank_storage = 0.0;
    } else {
        Eigen::MatrixXd eig_tile = TiledArray::eigen_map(tile, X, kl);
        Eigen::MatrixXd L, R;
        bool full_rank = algebra::Decompose_Matrix(eig_tile, L, R, 1e-8);

        if (!full_rank) {
            low_rank_storage = double(L.size() + R.size()) * 8.0 / 1.0e9;
        }
    }

    vec[ordinal] = {{low_rank_storage, sparse_storage, dense_storage}};
}

template <typename Pmap, typename SharedEnginePool, typename TileFunctor,
          typename... Bases>
std::vector<std::array<double, 3>>
compute_tiles(madness::World &world, Pmap &p,
              TiledArray::TiledRange const &trange, SharedEnginePool engines,
              TileFunctor tf, Bases &&... bases) {

    std::vector<std::array<double, 3>> tiles(p->size());
    auto shared_tiles = std::make_shared<decltype(tiles)>(std::move(tiles));

    madness::Range<decltype(p->begin())> m_range{p->begin(), p->end()};

    world.taskq.for_each(m_range, [=](decltype(p->begin()) it) {
                             call_tf(shared_tiles, engines, trange, it, tf,
                                     (&bases)...);
                             return true;
                         }).get();

    return *shared_tiles;
}

template <typename... Bases>
TiledArray::TiledRange get_trange(Bases &&... bases) {
    std::vector<TiledArray::TiledRange1> trange1_collector{
        bases.create_flattend_trange1()...};
    return TiledArray::TiledRange(trange1_collector.begin(),
                                  trange1_collector.end());
}

template <unsigned long order>
std::shared_ptr<TiledArray::Pmap>
create_pmap(madness::World &world, TiledArray::TiledRange const &trange) {
    return TiledArray::SparsePolicy::default_pmap(world,
                                                  trange.tiles().volume());
}

template <typename SharedEnginePool, typename TileFunctor, typename... Bases>
void Compute_Storage(madness::World &world, SharedEnginePool engines,
                     TileFunctor tf, Bases &&... bases) {

    constexpr auto tensor_order = sizeof...(bases);
    using TileType = typename TileFunctor::TileType;

    auto trange = get_trange(bases...);

    auto pmap = create_pmap<tensor_order>(world, trange);

    auto storage_info = compute_tiles(world, pmap, trange, std::move(engines),
                                      tf, std::forward<Bases>(bases)...);

    std::array<double, 3> out = {{0.0, 0.0, 0.0}};
    for (auto const &group : storage_info) {
        out[0] += group[0];
        out[1] += group[1];
        out[2] += group[2];
    }
    world.gop.sum(&out[0], 3);

    if (world.rank() == 0) {
        std::cout << "Low rank storage was " << out[0] << "GB\n"
                  << "Sparse storage was " << out[1] << "GB\n"
                  << "Dense storage was " << out[2] << "GB" << std::endl;
    }
}

} // namespace integrals
} // namespace tcc


#endif // TCC_INTEGRALS_SPARSETASKINTEGRALS_H
