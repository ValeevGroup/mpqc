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

namespace mpqc {
namespace integrals {

namespace detail {

// Specialize construction based on policy.
template <typename E, unsigned long N, typename Op, typename Policy>
struct compute_integrals;

// Dense array specialization
template <typename E, unsigned long N, typename Op>
struct compute_integrals<E, N, Op, DnPolicy> {
    DArray<N, Ttype<Op>, DnPolicy>
    operator()(mad::World &world, ShrPool<E> const &engines,
               Barray<N> const &bases, Op op) {

        using Tile = Ttype<Op>;
        DArray<N, Tile, DnPolicy> out(world, create_trange(bases));

        // Get Trange ptr for tasks will be used within tasks.  This assumes
        // that the trange outlives all of the tasks, otherwise this should be
        // a shared ptr.
        auto t_ptr = &(out.trange());

        // Shared ptr to bases
        auto shared_bases = std::make_shared<Barray<N>>(bases);


        auto op_wrapper =
               make_unscreened_op_invoke(t_ptr, engines, shared_bases, op);

        auto const &pmap = *(out.get_pmap());
        for (auto const ord : pmap) {
            mad::Future<Tile> tile = world.taskq.add(op_wrapper, ord);
            out.set(ord, tile);
        }

        return out;
    }
};

// Sparse Policy Specialization
template <typename E, unsigned long N, typename Op>
struct compute_integrals<E, N, Op, SpPolicy> {

    DArray<N, Ttype<Op>, SpPolicy>
    operator()(mad::World &world, ShrPool<E> const &engines,
               Barray<N> const &bases, Op op) {

        using Tile = Ttype<Op>;

        auto trange = create_trange(bases);
        const auto tvolume = trange.tiles().volume();

        auto shr_bases = std::make_shared<Barray<N>>(bases);

        // TODO make vector only hold as many tiles are on this node.
        std::vector<std::pair<unsigned long, Tile>> tiles(tvolume);
        TA::TensorF tile_norms(trange.tiles(), 0.0);

        auto const trange_ptr = &trange;
        auto norm_ptr = &tile_norms;
        auto tiles_ptr = &tiles;
        auto op_wrapper = [=](int64_t ord) { // copy ptrs by value

            auto op_invoker = make_unscreened_op_invoke(trange_ptr, engines,
                                                        shr_bases, op);
            auto ta_tile = op_invoker.integrals(ord);

            const auto tile_volume = ta_tile.range().volume();
            const auto tile_norm = ta_tile.norm();
            bool save_norm = tile_norm >= tile_volume * SpShapeF::threshold();

            auto &norms = *norm_ptr;
            auto &tiles = *tiles_ptr;
            if (save_norm) {
                auto tile = op_invoker.apply(std::move(ta_tile));
                norms[ord] = tile_norm;
                tiles[ord] = std::make_pair(ord, std::move(tile));
            } else {
                tiles[ord] = std::make_pair(ord, Tile());
            }
        };

        auto pmap = SpPolicy::default_pmap(world, tvolume);

        for (auto const ord : *pmap) {
            world.taskq.add(op_wrapper, ord);
        }
        world.gop.fence();

        TA::SparseShape<float> shape(world, tile_norms, trange);
        DArray<N, Tile, SpPolicy> out(world, trange, shape, pmap);

        set_array(tiles, out);
        out.truncate();

        return out;
    }
};

} // namespace detail

/*! \brief Construct integral tensors in parallel.
 *
 */
template <typename Policy = SpPolicy, typename E, unsigned long N, typename Op>
DArray<N, detail::Ttype<Op>, Policy>
TaskInts(mad::World &world, ShrPool<E> const &engines, Barray<N> const &bases,
         Op op) {
    return detail::compute_integrals<E, N, Op, Policy>()(world, engines, bases,
                                                         op);
}


} // namespace integrals
} // namespace mpqc


#endif // MPQC_INTEGRALS_TASKINTEGRALS_H
