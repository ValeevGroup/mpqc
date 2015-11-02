//
// task_integrals.h
//
// Copyright (C) 2015 Drew Lewis
// Maintainer Drew Lewis
//
// Based on file task_integrals.hpp from mpqc
//

#pragma once

#ifndef MPQC_INTEGRALS_TASKINTEGRALS_H
#define MPQC_INTEGRALS_TASKINTEGRALS_H

#include "task_integrals_common.h"
// #include "task_integrals_op_invoker.h"
#include "screening/screen_base.h"
#include "integral_builder.h"

namespace mpqc {
namespace integrals {

/*! \brief Construct sparse integral tensors in parallel.
 *
 * \param op needs to be a function or functor that takes a TA::TensorD && and 
 * returns any valid tile type. 
 *
 * \param screen needs to be a type derived from Screener
 */
template <typename ScreenerType, typename E, unsigned long N, typename Op>
DArray<N, detail::Ttype<Op>, SpPolicy>
sparse_integrals(mad::World &world, E const &engine, Barray<N> const &bases,
                 Op op, ScreenerType screen) {

    using Tile = detail::Ttype<Op>;

    auto trange = detail::create_trange(bases);
    const auto tvolume = trange.tiles().volume();
    std::vector<std::pair<unsigned long, Tile>> tiles(tvolume);
    TA::TensorF tile_norms(trange.tiles(), 0.0);

    auto builder_ptr = std::make_shared<IntegralBuilder<N, E, Op>>(
          make_integral_builder(world, engine, bases, std::move(screen), op));

    auto task_f = [=](int64_t ord, detail::IdxVec idx, TA::Range rng,
                      TA::TensorF *tile_norms_ptr, Tile *out_tile) {

        auto &builder = *builder_ptr;
        auto ta_tile = builder.integrals(idx, std::move(rng));

        const auto tile_volume = ta_tile.range().volume();
        const auto tile_norm = ta_tile.norm();

        // Keep tile if it was significant.
        bool save_norm = tile_norm >= tile_volume * SpShapeF::threshold();
        if (save_norm) {
            *out_tile = builder.op(std::move(ta_tile));

            auto &norms = *tile_norms_ptr;
            norms[ord] = tile_norm;
        }
    };

    auto pmap = SpPolicy::default_pmap(world, tvolume);
    for (auto const ord : *pmap) {
        tiles[ord].first = ord;
        detail::IdxVec idx = trange.tiles().idx(ord);
        world.taskq.add(task_f, ord, idx, trange.make_tile_range(ord),
                        &tile_norms, &tiles[ord].second);
    }
    world.gop.fence();

    TA::SparseShape<float> shape(world, tile_norms, trange);
    DArray<N, Tile, SpPolicy> out(world, trange, shape, pmap);

    detail::set_array(tiles, out);
    out.truncate();

    return out;
}

/*! \brief Construct a dense integral tensor in parallel.
 *
 */
template <typename E, unsigned long N, typename Op>
DArray<N, detail::Ttype<Op>, DnPolicy>
dense_integrals(mad::World &world, E const &engine, Barray<N> const &bases,
                Op op) {

    using Tile = detail::Ttype<Op>;
    DArray<N, Tile, DnPolicy> out(world, detail::create_trange(bases));

    auto builder = std::make_shared<IntegralBuilder<N, E, Op>>(
          make_integral_builder(world, engine, bases, Screener(), op));

    // builder is shared_ptr so just capture it by copy.
    auto task_func = [=](detail::IdxVec const &idx, TA::Range rng) {
        return builder->operator()(idx, std::move(rng));
    };

    auto const &trange = out.trange();
    auto const &pmap = *(out.get_pmap());
    for (auto const ord : pmap) {
        detail::IdxVec idx = trange.tiles().idx(ord);

        auto range = trange.make_tile_range(ord);
        mad::Future<Tile> tile
              = world.taskq.add(task_func, idx, range);

        out.set(ord, tile);
    }

    return out;
}

} // namespace integrals
} // namespace mpqc


#endif // MPQC_INTEGRALS_TASKINTEGRALS_H
