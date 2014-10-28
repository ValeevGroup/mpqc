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

namespace tcc {
namespace integrals {

template <typename It, typename TileFunctor>
void tile_task(It first, It last, engine, TileFunctor tf);

    template <typename It, typename TileFunctor>
    void generate_workpile(It begin, It end, engine, TileFunctor tf);


template <typename A, typename TileFunctor>
void create_tiles(A array, engine, TileFunctor tf);

template <typename TileType>
TiledArray::Array<double, order, TileType> create_array(world, trange_functor) {
    // Generate trange
    // create array
    // return array
}

template <typename TileFunctor>
TiledArray::Array<double, order, TileFunctor::TileType>
integrals(madness::World &world, engine, trange_functor, TileFunctor tf) {
    // Get integral engine
    // determine rank

    // create a TA::Array
    // fill array with integrals

    // return array
}

} // namespace integrals
} // namespace tcc


#endif // TCC_INTEGRALS_TASKINTEGRALS_H
