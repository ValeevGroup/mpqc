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

namespace mpqc {
namespace integrals {

namespace detail {

// Task wrapper which computes a TA::Tensor<double> containing integrals and
// passes it to Op. Takes the tile ordinal as a parameter.
template <typename E, unsigned long N, typename Op>
struct op_invoke {
    const TRange *const trange_ptr_;
    ShrPool<E> engines_;
    ShrBases<N> bases_;
    Op op_;

    op_invoke(const TRange *const tp, ShrPool<E> const &es,
              ShrBases<N> const &bs, Op op)
            : trange_ptr_(tp), engines_(es), bases_(bs), op_(op) {}

    TA::TensorD ta_integrals(int64_t ord) {

        std::array<basis::ShellVec const *, N> cluster_ptrs;
        auto const &idx = trange_ptr_->tiles().idx(ord);
        for (auto i = 0ul; i < N; ++i) {
            cluster_ptrs[i] = &(bases_->operator[](i).cluster_shells()[idx[i]]);
        }

        return integral_kernel(engines_->local(),
                               trange_ptr_->make_tile_range(ord), cluster_ptrs);
    }

    Ttype<Op> operator()(int64_t ord) { return op_(ta_integrals(ord)); }

    Ttype<Op> apply(TA::TensorD &&t) { return op_(std::move(t)); }
};

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

        // Get reference to pmap
        auto const &pmap = *(out.get_pmap());

        auto op_wrapper = op_invoke<E, N, Op>(t_ptr, engines, shared_bases, op);

        // For each ordinal create a task
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

        // Shared ptr to bases
        auto shared_bases = std::make_shared<Barray<N>>(bases);

        std::vector<std::pair<unsigned long, Tile>> tiles(tvolume);
        TA::TensorF tile_norms(trange.tiles(), 0.0);

        // Need to pass ptrs to the task function
        auto t_ptr = &trange;
        auto tn_ptr = &tile_norms;
        auto tile_vec_ptr = &tiles;
        auto op_wrapper = [=](int64_t ord) {

            // Compute Tile in TA::Tensor Form
            auto op_invoker
                  = op_invoke<E, N, Op>(t_ptr, engines, shared_bases, op);

            // Doing it this way potentially saves us from computing Op, which
            // could save a lot when Op is expensive. It also lets us get
            // better norm estimates.
            auto ta_tile = op_invoker.ta_integrals(ord);

            // Get volume and norm
            const auto tile_volume = ta_tile.range().volume();
            const auto tile_norm = ta_tile.norm();
            bool save_norm = tile_norm >= tile_volume * SpShapeF::threshold();

            if (save_norm) {
                tn_ptr->operator[](ord) = tile_norm;
                auto tile = op_invoker.apply(std::move(ta_tile));
                tile_vec_ptr->operator[](ord)
                      = std::make_pair(ord, std::move(tile));
            } else {
                tile_vec_ptr->operator[](ord) = std::make_pair(ord, Tile());
            }
        };

        auto pmap = SpPolicy::default_pmap(world, tvolume);

        // For each ordinal create a task
        for (auto const ord : *pmap) {
            world.taskq.add(op_wrapper, ord);
        }
        // Make sure all tiles are evaluated.
        world.gop.fence();

        TA::SparseShape<float> shape(world, tile_norms, trange);
        DArray<N, Tile, SpPolicy> out(world, trange, shape, pmap);

        for (auto &&tile : tiles) {
            const auto ord = tile.first;
            if(!out.is_zero(ord)){
                assert(!tile.second.empty());
                out.set(ord, std::move(tile.second));
            }
        }

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
