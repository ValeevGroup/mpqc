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
#include "../include/libint2.h"

#include "../basis/basis.h"

namespace tcc {
namespace integrals {

namespace detail {
class BtasTileFunctor {
  public:
    // using TileType = btas::Tensor<double>;
    using TileType = TiledArray::Tensor<double>;

    template <typename Index, typename Engine>
    TileType operator()(Index index, basis::Basis const *basis, Engine engine) {
        TileType tile{};
        return tile;
    }
}; // class BtasTileFunctor
} // namespace detail

template <unsigned long order, typename TileFunctor>
TiledArray::Array<double, order, typename TileFunctor::TileType>
create_array(madness::World &world, basis::Basis const &basis,
             TileFunctor const &tf) {

    std::vector<TiledArray::TiledRange1> trange1_collector;
    trange1_collector.reserve(order);

    const auto trange1 = basis.create_trange1();
    for (auto i = 0ul; i < order; ++i) {
        trange1_collector.push_back(trange1);
    }

    TiledArray::TiledRange trange(trange1_collector.begin(),
                                  trange1_collector.end());

    using TileType = typename TileFunctor::TileType;
    return TiledArray::Array<double, order, TileType>(world, trange);
}


template <typename It, typename Engine, typename TileFunctor>
void create_tile(It it, EnginePool<Engine> *engines, basis::Basis const *basis,
                 TileFunctor tf) {
    auto local_engine = engines->local();
    *it = tf(it.index(), basis, local_engine);
}


template <typename Array, typename Engine, typename TileFunctor>
void initialize_tiles(Array &a, EnginePool<Engine> engines,
                      basis::Basis const &basis, TileFunctor tf) {

    const auto end = a.end();
    for (auto it = a.begin(); it != end; ++it) {
        create_tile(it, &engines, &basis, tf);
    }
}


template <typename Engine,
          typename TileFunctor = typename detail::BtasTileFunctor>
TiledArray::Array<double, integrals::pool_order<EnginePool<Engine>>(),
                  typename TileFunctor::TileType>
Integrals(madness::World &world, EnginePool<Engine> engines,
          basis::Basis const &basis, TileFunctor tf = TileFunctor{}) {
    constexpr auto tensor_order = integrals::pool_order<EnginePool<Engine>>();
    auto array = create_array<tensor_order>(world, basis, tf);
    initialize_tiles(array, std::move(engines), basis, tf);
    return array;
}

} // namespace integrals
} // namespace tcc


#endif // TCC_INTEGRALS_TASKINTEGRALS_H
