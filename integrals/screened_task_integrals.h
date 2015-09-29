// Copyright (C) 2015 Drew Lewis
// Maintainer Drew Lewis
//
// Based on file task_integrals.hpp from mpqc
//

#pragma once

#ifndef MPQC_INTEGRALS_SCREENEDTASKINTEGRALS_H
#define MPQC_INTEGRALS_SCREENEDTASKINTEGRALS_H

#include "task_integrals_common.h"
#include "integral_screening_matrices.h"

#include "../integrals/integral_engine_pool.h"
#include "../include/eigen.h"

#include <memory>
#include <algorithm>

namespace mpqc {
namespace integrals {

namespace detail {

// Task wrapper which computes a TA::Tensor<double> containing integrals and
// passes it to Op. Takes the tile ordinal as a parameter.
template <typename E, unsigned long N, typename Op>
struct op_invoke_sc;

template <typename E, typename Op>
struct op_invoke_sc<E, 3, Op> {
    const TRange *const trange_ptr_;
    ShrPool<E> engines_;
    ShrBases<3> bases_;
    const MatrixD *X_shells_;
    const MatrixD *ab_shells_;
    Op op_;

    op_invoke_sc(const TRange *const tp, ShrPool<E> const &es,
                 ShrBases<3> const &bs, const MatrixD *xshs,
                 const MatrixD *abshs, Op op)
            : trange_ptr_(tp),
              engines_(es),
              bases_(bs),
              X_shells_(xshs),
              ab_shells_(abshs),
              op_(op) {}

    TA::TensorD ta_integrals(int64_t ord) {

        std::array<tcc::basis::ClusterShells const *, 3> cluster_ptrs;
        auto const &idx = trange_ptr_->tiles().idx(ord);
        for (auto i = 0ul; i < 3ul; ++i) {
            cluster_ptrs[i] = &(bases_->operator[](i).cluster_shells()[idx[i]]);
        }

        return integral_kernel(engines_->local(),
                               trange_ptr_->make_tile_range(ord), cluster_ptrs,
                               *X_shells_, *ab_shells_);
    }

    Ttype<Op> operator()(int64_t ord) { return op_(ta_integrals(ord)); }

    Ttype<Op> apply(TA::TensorD &&t) { return op_(std::move(t)); }
};

template <typename E, typename Op>
DArray<3, Ttype<Op>, SpPolicy>
compute_screened_integrals(mad::World &world, ShrPool<E> &engines,
                           Barray<3> const &bases, Op op) {
    // Depends on auxilary basis being in position 0.
    const auto Q_X = screening_matrix_X(engines, bases[0].cluster_shells());
    auto const &cls_X = Q_X.cluster_screening;

    const auto Q_ab = screening_matrix_ab(engines, bases[1].cluster_shells());
    auto const &cls_ab = Q_ab.cluster_screening;

    using Tile = Ttype<Op>;

    auto trange = create_trange(bases);
    const auto tvolume = trange.tiles().volume();

    // Shared ptr to bases
    auto shared_bases = std::make_shared<Barray<3>>(bases);

    std::vector<std::pair<unsigned long, Tile>> tiles(tvolume);
    TA::TensorF tile_norms(trange.tiles(), 0.0);

    // Need to pass ptrs to the task function
    auto t_ptr = &trange;
    auto tn_ptr = &tile_norms;
    auto tile_vec_ptr = &tiles;
    auto const *sh_X = &(Q_X.shell_screening);
    auto const *sh_ab = &(Q_ab.shell_screening);

    auto op_wrapper = [=](int64_t ord) {

        // Compute Tile in TA::Tensor Form
        auto op_invoker = op_invoke_sc<E, 3, Op>(t_ptr, engines, shared_bases,
                                                 sh_X, sh_ab, op);

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

            // Only compute Op if norm was large. Also don't move this
            // inside of the lock since it may be expensive and take a long
            // time.
            auto tile = op_invoker.apply(std::move(ta_tile));

            // Lock to save tile, I don't expect much contention
            tile_vec_ptr->operator[](ord)
                  = std::make_pair(ord, std::move(tile));
        } else {
            tile_vec_ptr->operator[](ord) = std::make_pair(ord, Tile());
        }
    };

    auto pmap = SpPolicy::default_pmap(world, tvolume);

    // For each ordinal create a task
    for (auto const ord : *pmap) {
        auto const &idx = trange.tiles().idx(ord);
        const auto guess_norm = cls_X(idx[0]) * cls_ab(idx[1], idx[2]);
        if (guess_norm > 1e-10) {
            world.taskq.add(op_wrapper, ord);
        } else {
            tiles[ord] = std::make_pair(ord, Tile());
        }
    }
    // Make sure all tiles are evaluated.
    world.gop.fence();

    TA::SparseShape<float> shape(world, tile_norms, trange);
    DArray<3, Tile, SpPolicy> out(world, trange, shape, pmap);

    for (auto &&tile : tiles) {
        const auto ord = tile.first;
        if (!out.is_zero(ord)) {
            assert(!tile.second.empty());
            out.set(ord, std::move(tile.second));
        }
    }

    return out;
}

template <typename E, typename Op>
DArray<4, Ttype<Op>, SpPolicy>
compute_screened_integrals(mad::World &world, ShrPool<E> &engines,
                           Barray<4> const &bases, Op op) {
    // Assumes basis 1 and basis 2 are the same
    auto Q_ab = screening_matrix_ab(engines, bases[0].cluster_shells());
}

} // namespace detail

/*! \brief Construct integral tensors in parallel with screening.
 *
 */
template <typename E, unsigned long N, typename Op>
DArray<N, detail::Ttype<Op>, SpPolicy>
ScreenedTaskInts(mad::World &world, ShrPool<E> &engines, Barray<N> const &bases,
                 Op op) {
    static_assert(N == 3 || N == 4,
                  "Screening only avalible for 3 and 4 center integrals.");
    return detail::compute_screened_integrals(world, engines, bases, op);
}


} // namespace integrals
} // namespace mpqc


#endif // MPQC_INTEGRALS_SCREENEDTASKINTEGRALS_H
