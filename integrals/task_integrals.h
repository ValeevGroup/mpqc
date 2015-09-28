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

#include "../common/typedefs.h"
#include "../include/tiledarray.h"
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
            : trange_ptr_{tp}, engines_{es}, bases_{bs}, op_{op} {}

    Ttype<Op> operator()(int64_t ord) {
        return op_(TA::TensorD(trange_ptr_->make_tile_range(ord), 0.0));
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

        // Get reference to pmap
        auto const &pmap = *(out.get_pmap());
        // For each ordinal create a task
        for (auto const ord : pmap) {
            mad::Future<Tile> tile = world.taskq.add(op_pass_through<E, N, Op>(
                  t_ptr, engines, std::make_shared<Barray<N>>(bases), op), ord);
            out.set(ord, tile);
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
