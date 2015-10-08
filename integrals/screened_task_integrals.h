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

template <unsigned long N>
struct compute_cluster;

template <>
struct compute_cluster<3> {
    template <typename Index>
    bool
    operator()(Index const &idx, int64_t ord, MatrixD const &Cx,
               MatrixD const &Cab, int64_t volume, TA::TensorF &norms) const {

        const auto estimate = Cx(idx[0]) * Cab(idx[1], idx[2]);

        norms[ord] = estimate;
        if (estimate >= volume * SpShapeF::threshold()) {
            return true;
        } else {
            return false;
        }
    }
};

template <typename E, typename Op>
DArray<3, Ttype<Op>, SpPolicy>
compute_screened_integrals(mad::World &world, ShrPool<E> &engines,
                           Barray<3> const &bases, Op op) {
    // Depends on auxiliary basis being in position 0.
    const auto Q_X
          = screening_matrix_X(world, engines, bases[0].cluster_shells());
    auto const &cls_X = Q_X.cluster_screening;

    const auto Q_ab
          = screening_matrix_ab(world, engines, bases[1].cluster_shells());
    auto const &cls_ab = Q_ab.cluster_screening;

    using Tile = Ttype<Op>;

    auto trange = create_trange(bases);
    const auto tvolume = trange.tiles().volume();

    // Shared ptr to bases
    auto shared_bases = std::make_shared<Barray<3>>(bases);

    auto pmap = SpPolicy::default_pmap(world, tvolume);

    OrdTileVec<Tile> tiles(tvolume);
    TA::TensorF tile_norms(trange.tiles(), 0.0);


    const auto estimate = compute_cluster<3>();
    const auto t_ptr = &trange;
    auto sh_vecs_X = std::make_shared<decltype(Q_X.shell_screenings)>(
          Q_X.shell_screenings);
    auto sh_vecs_ab = std::make_shared<decltype(Q_ab.shell_screenings)>(
          Q_ab.shell_screenings);
    auto norms_ptr = &tile_norms;
    auto tiles_ptr = &tiles;

    for (auto const ord : *pmap) {
        auto const &idx = trange.tiles().idx(ord);
        const auto tile_volume = trange.make_tile_range(ord).volume();

        if (estimate(idx, ord, cls_X, cls_ab, tile_volume, tile_norms)) {

            auto sh_vec_X_ptr = &(sh_vecs_X->operator[](idx[0]).back());
            auto sh_X_ptr = std::shared_ptr<MatrixD>(sh_vecs_X, sh_vec_X_ptr);

            auto sh_vec_ab_ptr = &(sh_vecs_ab->operator[](idx[1])[idx[2]]);
            auto sh_ab_ptr
                  = std::shared_ptr<MatrixD>(sh_vecs_ab, sh_vec_ab_ptr);

            auto op_wrapper = [=](int64_t tile_ord) {

                IdxVec idx = t_ptr->tiles().idx(ord);

                auto shell_vecs = get_shells(idx, shared_bases);

                const auto invoker = make_schwartz_op_invoke(
                      std::move(idx), engines, std::move(shell_vecs), op,
                      std::move(sh_X_ptr), std::move(sh_ab_ptr));

                auto range = t_ptr->make_tile_range(ord);

                auto tile = invoker.integrals(std::move(range));

                const auto vol = tile.range().volume();
                const auto norm = tile.norm();
                bool save = norm >= vol * SpShapeF::threshold();

                auto &norms = *norms_ptr;
                auto &tiles = *tiles_ptr;
                if (save) {
                    auto op_tile = invoker.apply(std::move(tile));
                    norms[tile_ord] = norm;
                    tiles[tile_ord] = std::make_pair(ord, std::move(op_tile));
                } else {
                    norms[tile_ord] = 0.0;
                    tiles[tile_ord] = std::make_pair(ord, Tile());
                }
            };

            world.taskq.add(op_wrapper, ord);
        } else {
            tiles[ord] = std::make_pair(ord, Tile());
        }
    }
    world.gop.fence();

    TA::SparseShape<float> shape(world, tile_norms, trange);
    DArray<3, Tile, SpPolicy> out(world, trange, shape, pmap);

    set_array(tiles, out);
    out.truncate();
    world.gop.fence();

    return out;
}

// template <typename E, typename Op>
// DArray<4, Ttype<Op>, SpPolicy>
// compute_screened_integrals(mad::World &world, ShrPool<E> &engines,
//                            Barray<4> const &bases, Op op) {
//     // Assumes basis 1 and basis 2 are the same
//     auto Q_ab = screening_matrix_ab(engines, bases[0].cluster_shells());
// }

} // namespace detail

/*! \brief Construct integral tensors in parallel with screening.
 *
 */
template <typename E, unsigned long N, typename Op>
DArray<N, detail::Ttype<Op>, SpPolicy>
ScreenedTaskInts(mad::World &world, ShrPool<E> &engines, Barray<N> const &bases,
                 Op op) {
    static_assert(N == 3, "Screening only avalible for 3 center integrals.  At "
                          "the moment");
    return detail::compute_screened_integrals(world, engines, bases, op);
}


} // namespace integrals
} // namespace mpqc


#endif // MPQC_INTEGRALS_SCREENEDTASKINTEGRALS_H
