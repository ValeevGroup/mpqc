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

#include "../basis/cluster_shells.h"
#include "integral_engine_pool.h"
#include "btas_compute_functors.h"

#include "../basis/basis.h"

#include "../tensor/tile_pimpl_devel.h"

namespace tcc {
namespace integrals {

template <typename Pmap, typename SharedEnginePool, std::size_t N, typename TF>
std::vector<std::pair<std::size_t, typename TF::TileType>>
compute_tiles(Pmap &p, TiledArray::TiledRange const &trange,
              SharedEnginePool engines,
              std::array<basis::Basis, N> const &bases, TF tf) {


    std::array<std::shared_ptr<std::vector<basis::ClusterShells>>, N>
        shell_ptrs;

    for (auto i = 0ul; i < N; ++i) {
        shell_ptrs[i] = std::make_shared<std::vector<basis::ClusterShells>>(
            bases[i].cluster_shells());
    }

    auto call_tf = [=](TiledArray::Range range, TiledArray::Range::index idx,
                       decltype(shell_ptrs) shells) {
        return std::make_pair(trange.tiles().ordinal(idx),
                              tf(std::move(range), idx, engines, shells));
    };

    std::vector<std::pair<std::size_t, typename TF::TileType>> tiles;
    const auto end = p->end();
    for (auto it = p->begin(); it != end; ++it) {
        tiles.push_back(call_tf(trange.make_tile_range(*it),
                                trange.tiles().idx(*it), shell_ptrs));
    }

    return tiles;
}

template <std::size_t N>
TiledArray::TiledRange create_trange(std::array<basis::Basis, N> basis_array) {
    std::vector<TiledArray::TiledRange1> trange1_collector;
    trange1_collector.reserve(N);

    for (auto i = 0ul; i < N; ++i) {
        trange1_collector.push_back(basis_array[i].create_flattend_trange1());
    }

    return TiledArray::TiledRange(trange1_collector.begin(),
                                  trange1_collector.end());
}


template <unsigned long order>
std::shared_ptr<TiledArray::Pmap>
create_pmap(madness::World &world, TiledArray::TiledRange const &trange) {
    return TiledArray::SparsePolicy::default_pmap(world,
                                                  trange.tiles().volume());
}

template <typename SharedEnginePool, std::size_t N,
          typename TF = compute_functors::BtasTileFunctor<double>>
TiledArray::Array<double, N, typename TF::TileType, TiledArray::SparsePolicy>
SparseIntegrals(madness::World &world, SharedEnginePool engines,
                std::array<basis::Basis, N> const &bases, TF tf = TF{}) {

    using TileType = typename TF::TileType;

    auto basis = bases[0];
    auto trange = create_trange(bases);
    auto pmap = create_pmap<N>(world, trange);

    auto computed_tiles
        = compute_tiles(pmap, trange, std::move(engines), bases, tf);

    TiledArray::Tensor<float> tile_norms(trange.tiles(), 0.0);
    for (auto &ord_tile_pair : computed_tiles) {
        tile_norms[ord_tile_pair.first] = ord_tile_pair.second.norm();
    }
    TiledArray::SparseShape<float> sparse_S(world, tile_norms, trange);

    TiledArray::Array<double, N, typename TF::TileType,
                      TiledArray::SparsePolicy> ta_array(world, trange,
                                                         sparse_S);

    for (auto &pair : computed_tiles) {
        if (!ta_array.is_zero(pair.first)) {
            ta_array.set(pair.first, pair.second);
        }
    }

    return ta_array;
}

} // namespace integrals
} // namespace tcc


#endif // TCC_INTEGRALS_SPARSETASKINTEGRALS_H
