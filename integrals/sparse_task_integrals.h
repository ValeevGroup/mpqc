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

template <typename Pmap, typename SharedEnginePool, typename TileFunctor>
std::vector<std::pair<std::size_t, typename TileFunctor::TileType>>
compute_tiles(Pmap &p, TiledArray::TiledRange const &trange,
              SharedEnginePool engines, basis::Basis const &basis,
              TileFunctor tf) {
    std::vector<std::pair<std::size_t, typename TileFunctor::TileType>> tiles;

    auto call_tf = [=](TiledArray::Range range, TiledArray::Range::index idx,
                       basis::Basis const *basis_) {
        return std::make_pair(trange.tiles().ordinal(idx),
                              tf(std::move(range), idx, basis_, engines));
    };

    const auto end = p->end();
    for (auto it = p->begin(); it != end; ++it) {
        tiles.push_back(call_tf(trange.make_tile_range(*it), trange.tiles().idx(*it), &basis));
    }
    return tiles;
}

template <std::size_t order>
TiledArray::TiledRange get_trange(basis::Basis const &basis) {
    std::vector<TiledArray::TiledRange1> trange1_collector;
    trange1_collector.reserve(order);

    const auto trange1 = basis.create_flattend_trange1();
    for (auto i = 0ul; i < order; ++i) {
        trange1_collector.push_back(trange1);
    }

    return TiledArray::TiledRange(trange1_collector.begin(),
                                  trange1_collector.end());
}


template <unsigned long order>
std::shared_ptr<TiledArray::Pmap>
create_pmap(madness::World &world, TiledArray::TiledRange const &trange) {
    return TiledArray::SparsePolicy::default_pmap(world, trange.tiles().volume());
}

template <typename SharedEnginePool,
          typename TileFunctor = compute_functors::BtasTileFunctor<double>>
TiledArray::Array<double, integrals::pool_order<SharedEnginePool>(),
                  typename TileFunctor::TileType, TiledArray::SparsePolicy>
SparseIntegrals(madness::World &world, SharedEnginePool engines,
                basis::Basis const &basis, TileFunctor tf = TileFunctor{}) {

    constexpr auto tensor_order = integrals::pool_order<decltype(engines)>();
    using TileType = typename TileFunctor::TileType;

    auto trange = get_trange<tensor_order>(basis);

    auto pmap = create_pmap<tensor_order>(world, trange);

    auto computed_tiles
        = compute_tiles(pmap, trange, std::move(engines), basis, tf);

    TiledArray::Tensor<float> tile_norms(trange.tiles(), 0.0);
    for (auto &ord_tile_pair : computed_tiles) {

        tile_norms[ord_tile_pair.first] = ord_tile_pair.second.norm();
    }
    TiledArray::SparseShape<float> sparse_S(world, tile_norms, trange);

    TiledArray::Array<double, tensor_order, typename TileFunctor::TileType,
                      TiledArray::SparsePolicy> array(world, trange, sparse_S);

    for (auto &pair : computed_tiles) {
        if(!array.is_zero(pair.first)){
            array.set(pair.first, pair.second);
        }
    }

    return array;
}

} // namespace integrals
} // namespace tcc


#endif // TCC_INTEGRALS_SPARSETASKINTEGRALS_H
