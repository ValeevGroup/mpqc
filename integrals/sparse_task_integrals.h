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
compute_tiles(Pmap &p, SharedEnginePool engines, basis::Basis const &basis,
                 TileFunctor tf) {
    std::vector<std::pair<std::size_t, typename TileFunctor::TileType>> tiles;

    auto call_tf = [=](decltype(p.begin()) it, basis::Basis const *basis_) {
        return std::make_pair(it.ordinal(), tf(it, basis_, engines));
    };

    const auto end = p.end();
    for (auto it = p.begin(); it != end; ++it) {
        tiles.push_back(call_tf(it, &basis));
        // a.get_world().taskq.add(call_tf, &basis);
    }
    return tiles;
}

template <unsigned long order>
std::shared_ptr<TiledArray::Pmap>
create_pmap(madness::World &world, basis::Basis const &basis) {

    std::vector<TiledArray::TiledRange1> trange1_collector;
    trange1_collector.reserve(order);

    const auto trange1 = basis.create_flattend_trange1();
    for (auto i = 0ul; i < order; ++i) {
        trange1_collector.push_back(trange1);
    }

    TiledArray::TiledRange trange(trange1_collector.begin(),
                                  trange1_collector.end());

    return TiledArray::SparsePolicy::default_pmap(world, trange.data().size());
}

template <typename SharedEnginePool,
          typename TileFunctor = compute_functors::BtasTileFunctor<double>>
TiledArray::Array<double, integrals::pool_order<SharedEnginePool>(),
                  typename TileFunctor::TileType, TiledArray::SparsePolicy>
SparseIntegrals(madness::World &world, SharedEnginePool engines,
                basis::Basis const &basis, TileFunctor tf = TileFunctor{}) {

    constexpr auto tensor_order = integrals::pool_order<decltype(engines)>();
    using TileType = typename TileFunctor::TileType;

    auto pmap = create_pmap<tensor_order>(world, basis);

    auto computed_tiles = compute_tiles(pmap, std::move(engines), basis, tf);

    TiledArray::Tensor<float> array_shape(pmap.trange().tiles(), 0.0);
    for (auto &ord_tile_pair : computed_tiles) {
        array_shape[ord_tile_pair.first] = ord_tile_pair.second.norm();
    }
    TiledArray::SparseShape<float> sparse_S(world, array_shape, pmap.trange());

    TiledArray::Array<double, tensor_order, typename TileFunctor::TileType,
                      TiledArray::SparsePolicy> array(world, pmap.trange(),
                                                      sparse_S);

    for (auto &pair : computed_tiles) {
        array.set(pair.first, pair.second);
    }

    return array;
}

} // namespace integrals
} // namespace tcc


#endif // TCC_INTEGRALS_SPARSETASKINTEGRALS_H
