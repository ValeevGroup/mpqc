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

#include "task_integrals_helper.h"

#include "../common/typedefs.h"

#include "../include/tiledarray.h"
#include "../include/tbb.h"

#include "../basis/basis.h"

#include <memory>
#include <array>

namespace mpqc {
namespace integrals {

template <typename E>
using ShrPool = std::shared_ptr<Epool<E>>;

template <unsigned long N>
using Barray = std::array<tcc::basis::Basis, N>;

namespace detail {

template <unsigned long N>
using ShrBases = std::shared_ptr<Barray<N>>;

template <typename Op>
using Ttype = decltype(std::declval<Op>()(std::declval<TA::TensorD>()));

// Capture by value since tasks must capture by value
template <typename E, unsigned long N, typename Op>
struct op_pass_through {
    const TRange *const trange_ptr_;
    ShrPool<E> engines_;
    ShrBases<N> bases_;
    Op op_;

    op_pass_through() = delete;
    op_pass_through(const TRange *const tp, ShrPool<E> const &es,
                    ShrBases<N> const &bs, Op op)
            : trange_ptr_(tp), engines_(es), bases_(bs), op_(op) {}

    op_pass_through(op_pass_through const &) = default;
    op_pass_through(op_pass_through &&) = default;

    Ttype<Op> operator()(int64_t ord) {

        // Get pointers to the shells that are needed for this tile.
        std::array<tcc::basis::ClusterShells const *, N> cluster_ptrs;
        auto const &idx = trange_ptr_->tiles().idx(ord);
        for (auto i = 0ul; i < N; ++i) {
            cluster_ptrs[i] = &(bases_->operator[](i).cluster_shells()[idx[i]]);
        }

        return op_(integral_kernel(engines_->local(),
                                   trange_ptr_->make_tile_range(ord),
                                   cluster_ptrs));
    }
};

// Create TRange from bases
template <std::size_t N>
TRange create_trange(Barray<N> const &basis_array) {

    std::vector<TRange1> trange1s;
    trange1s.reserve(N);

    for (auto i = 0ul; i < N; ++i) {
        trange1s.emplace_back(basis_array[i].create_flattend_trange1());
    }

    return TRange(trange1s.begin(), trange1s.end());
}

// Specialize construction based on policy.
template <typename E, unsigned long N, typename Op, typename Policy>
struct compute_integrals;

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

        auto op_wrapper
              = op_pass_through<E, N, Op>(t_ptr, engines, shared_bases, op);

        // For each ordinal create a task
        for (auto const ord : pmap) {
            mad::Future<Tile> tile = world.taskq.add(op_wrapper, ord);
            out.set(ord, tile);
        }


        return out;
    }
};

static tbb::spin_mutex task_int_mutex;

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
        auto t_ptr = &trange;

        std::vector<std::pair<unsigned long, Tile>> tiles;
        TA::TensorF tile_norms(trange.tiles(), 0.0);

        auto tn_ptr = &tile_norms;
        auto tile_vec_ptr = &tiles;
        auto op_wrapper = [=](int64_t ord) {
            auto tile = op_pass_through<E, N, Op>(t_ptr, engines, shared_bases,
                                                  op)(ord);

            const auto tile_volume = tile.range().volume();
            const auto tile_norm = tile.norm();

            if (tile_norm >= tile_volume
                             * TA::SparseShape<float>::threshold()) {
                tn_ptr->operator[](ord) = tile.norm();
                // Lock the vector to push back, I don't expect much contention
                tbb::spin_mutex::scoped_lock lock(task_int_mutex);
                tile_vec_ptr->emplace_back(std::make_pair(ord, tile));
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
        
        for(auto &&tile : tiles){
            const auto ord = tile.first;
            out.set(ord, std::move(tile.second));
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
         Op op, bool screening = false) {
        return detail::compute_integrals<E, N, Op, Policy>()(
              world, engines, bases, std::move(op));
}


} // namespace integrals
} // namespace mpqc


#endif // MPQC_INTEGRALS_TASKINTEGRALS_H
