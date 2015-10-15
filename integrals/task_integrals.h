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
#include "task_integrals_op_invoker.h"
#include "screening/screen_base.h"

namespace mpqc {
namespace integrals {

namespace detail {

template <typename E, typename Op, unsigned long N, typename ScreenOp>
void sparse_integral_task_function(int64_t ord, IdxVec idx, ShrPool<E> engs,
                                   TA::Range rng, ShrBases<N> shr_bases, Op op,
                                   TA::TensorF *tile_norms_ptr,
                                   detail::Ttype<Op> *out_tile) {

    auto shell_vecs = get_shells(idx, shr_bases);


    auto screen = ScreenOp()(idx, shr_bases, engs);

    auto op_invoker = make_op_invoke(std::move(idx), std::move(engs),
                                     std::move(shell_vecs), std::move(op),
                                     std::move(screen));

    auto ta_tile = op_invoker.integrals(std::move(rng));

    const auto tile_volume = ta_tile.range().volume();
    const auto tile_norm = ta_tile.norm();
    bool save_norm = tile_norm >= tile_volume * SpShapeF::threshold();

    auto &norms = *tile_norms_ptr;
    if (save_norm) {
        auto tile = op_invoker.apply(std::move(ta_tile));
        norms[ord] = tile_norm;
        *out_tile = std::move(tile);
    }
}

} // namespace detail

/*! \brief Construct sparse integral tensors in parallel.
 *
 */
template <typename ScreenOp = init_base_screen, typename E,
          unsigned long N, typename Op>
DArray<N, detail::Ttype<Op>, SpPolicy>
sparse_integrals(mad::World &world, ShrPool<E> const &engines,
                 Barray<N> const &bases, Op op) {

    using Tile = detail::Ttype<Op>;

    auto shr_bases = std::make_shared<Barray<N>>(bases);

    auto trange = detail::create_trange(bases);
    const auto tvolume = trange.tiles().volume();
    std::vector<std::pair<unsigned long, Tile>> tiles(tvolume);
    TA::TensorF tile_norms(trange.tiles(), 0.0);

    auto pmap = SpPolicy::default_pmap(world, tvolume);
    for (auto const ord : *pmap) {

        tiles[ord].first = ord;
        detail::IdxVec idx = trange.tiles().idx(ord);

        world.taskq.add(
              detail::sparse_integral_task_function<E, Op, N, ScreenOp>,
              ord, idx, engines, trange.make_tile_range(ord), shr_bases, op,
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
dense_integrals(mad::World &world, ShrPool<E> const &engines,
                Barray<N> const &bases, Op op) {

    using Tile = detail::Ttype<Op>;
    DArray<N, Tile, DnPolicy> out(world, detail::create_trange(bases));

    // Shared ptr to bases
    auto shared_bases = std::make_shared<Barray<N>>(bases);

    auto const &trange = out.trange();
    auto const &pmap = *(out.get_pmap());
    for (auto const ord : pmap) {

        detail::IdxVec idx = trange.tiles().idx(ord);
        auto shell_vecs = detail::get_shells(idx, shared_bases);
        auto op_wrapper
              = detail::make_op_invoke(idx, engines, std::move(shell_vecs), op);

        auto range = trange.make_tile_range(ord);
        mad::Future<Tile> tile = world.taskq.add(op_wrapper, std::move(range));

        out.set(ord, tile);
    }

    return out;
}

} // namespace integrals
} // namespace mpqc


#endif // MPQC_INTEGRALS_TASKINTEGRALS_H
