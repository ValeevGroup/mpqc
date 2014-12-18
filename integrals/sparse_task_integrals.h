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
std::vector<std::tuple<std::size_t, double, typename TileFunctor::TileType>>
compute_tiles(madness::World &world, Pmap &p,
              TiledArray::TiledRange const &trange, SharedEnginePool engines,
              basis::Basis const &basis, TileFunctor tf) {
    std::vector<std::tuple<std::size_t, double, typename TileFunctor::TileType>>
        tiles(p->size());
    auto shared_tiles = std::make_shared<decltype(tiles)>(std::move(tiles));

    auto call_tf = [=](std::shared_ptr<decltype(tiles)> tile_vec_ptr,
                       TiledArray::Range range, TiledArray::Range::index idx,
                       basis::Basis const *basis_) {
        auto &vec = *tile_vec_ptr;
        auto const ordinal = trange.tiles().ordinal(idx);
        auto tile = tf(std::move(range), idx, basis_, engines);
        auto norm = tile.norm();
        if (norm >= TiledArray::SparseShape<float>::threshold()) {
            vec[ordinal] = std::make_tuple(ordinal, norm, std::move(tile));
        } else {
            vec[ordinal] = std::make_tuple(
                ordinal, norm, typename TileFunctor::TileType{tile.range()});
        }
    };

    madness::Range<decltype(p->begin())> m_range{p->begin(), p->end()};

    world.taskq.for_each(m_range, [=](decltype(p->begin()) it) {
                             call_tf(shared_tiles, trange.make_tile_range(*it),
                                     trange.tiles().idx(*it), &basis);
                             return true;
                         }).get();

    return *shared_tiles;
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
    return TiledArray::SparsePolicy::default_pmap(world,
                                                  trange.tiles().volume());
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
        = compute_tiles(world, pmap, trange, std::move(engines), basis, tf);

    TiledArray::Tensor<float> tile_norms(trange.tiles(), 0.0);
    for (auto &group : computed_tiles) {
        tile_norms[std::get<0>(group)] = std::get<1>(group);
    }

    TiledArray::SparseShape<float> sparse_S(world, tile_norms, trange);

    TiledArray::Array<double, tensor_order, typename TileFunctor::TileType,
                      TiledArray::SparsePolicy> array(world, trange, sparse_S);

    for (auto &group : computed_tiles) {
        auto ordinal = std::get<0>(group);
        if (!array.is_zero(ordinal)) {
            array.set(ordinal, std::move(std::get<2>(group)));
        }
    }

    return array;
}

} // namespace integrals
} // namespace tcc


#endif // TCC_INTEGRALS_SPARSETASKINTEGRALS_H
