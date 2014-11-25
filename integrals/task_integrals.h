//
// task_integrals.h
//
// Copyright (C) 2014 Drew Lewis
// Maintainer Drew Lewis
//
// Based on file task_integrals.hpp from mpqc
//

#ifndef TCC_INTEGRALS_TASKINTEGRALS_H
#define TCC_INTEGRALS_TASKINTEGRALS_H

#include "../include/tiledarray.h"
#include "../include/btas.h"
#include "../include/libint2.h"

#include "../basis/cluster_shells.h"
#include "integral_engine_pool.h"
#include "compute_functors.h"

#include "../basis/basis.h"

namespace tcc {
namespace integrals {

template <typename Array, typename SharedEnginePool, typename TileFunctor>
void initialize_tiles(Array &a, SharedEnginePool engines,
                      basis::Basis const &basis, TileFunctor tf) {
    const auto end = a.end();
    for (auto it = a.begin(); it != end; ++it) {
        auto call_tf = [=](basis::Basis const *basis_) {
            *it = (tf(it, basis_, engines));
        };
        a.get_world().taskq.add(call_tf, &basis);
    }
}

template <unsigned long order, typename TileType>
TiledArray::Array<double, order, TileType>
create_array(madness::World &world, basis::Basis const &basis) {

    std::vector<TiledArray::TiledRange1> trange1_collector;
    trange1_collector.reserve(order);

    const auto trange1 = basis.create_flattend_trange1();
    for (auto i = 0ul; i < order; ++i) {
        trange1_collector.push_back(trange1);
    }

    TiledArray::TiledRange trange(trange1_collector.begin(),
                                  trange1_collector.end());

    return TiledArray::Array<double, order, TileType>(world, trange);
}

template <typename SharedEnginePool,
          typename TileFunctor = compute_functors::BtasTileFunctor<double>>
TiledArray::Array<double, integrals::pool_order<SharedEnginePool>(),
                  typename TileFunctor::TileType>
Integrals(madness::World &world, SharedEnginePool engines,
          basis::Basis const &basis, TileFunctor tf = TileFunctor{}) {

    constexpr auto tensor_order = integrals::pool_order<decltype(engines)>();
    using TileType = typename TileFunctor::TileType;

    auto array = create_array<tensor_order, TileType>(world, basis);

    initialize_tiles(array, std::move(engines), basis, tf);
    return array;
}

} // namespace integrals
} // namespace tcc


#endif // TCC_INTEGRALS_TASKINTEGRALS_H
