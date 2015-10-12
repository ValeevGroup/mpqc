#pragma once
#ifndef MPQC_INTEGRALS_DIRECTTASKINTEGRALS_H
#define MPQC_INTEGRALS_DIRECTTASKINTEGRALS_H

#include "direct_tile.h"
#include "task_integrals_op_invoker.h"
#include "task_integrals_common.h"
#include <limits>

namespace mpqc {
namespace integrals {

namespace detail {

template <typename E, unsigned long N, typename Op, typename ScreenOp>
void direct_tile_task(
      int64_t ord, TA::Range range, IdxVec idx, ShrPool<E> engines,
      ShrBases<N> bases, Op op, TA::TensorF *norm_ptr,
      std::vector<std::pair<int64_t, DirectTile<E, N, Op>>> *tiles) {

    const auto volume = range.volume();
    auto tile = make_direct_tile(std::move(range), std::move(idx),
                                 std::move(engines), std::move(bases),
                                 std::move(op));

    const auto temp_tile = typename DirectTile<E, N, Op>::eval_type(tile);
    const auto norm = temp_tile.norm();
    (*norm_ptr)[ord] = norm;

    if (norm >= volume * SpShapeF::threshold()) {
        (*tiles)[ord] = std::make_pair(ord, std::move(tile));
    } else {
        (*tiles)[ord].first = ord;
    }
}

} // namespace detail

/*! \brief Construct direct integral tensors in parallel with screening.
 *
 */
template <typename ScreenOp = init_base_screen, typename E, unsigned long N,
          typename Op>
DArray<N, DirectTile<E, N, Op>, SpPolicy>
direct_sparse_integrals(mad::World &world, ShrPool<E> &engines,
                        Barray<N> const &bases, Op op) {

    const auto trange = detail::create_trange(bases);
    const auto tvolume = trange.tiles().volume();
    TA::TensorF tile_norms(trange.tiles(), std::numeric_limits<double>::max());

    std::vector<std::pair<int64_t, DirectTile<E, N, Op>>> tiles(tvolume);

    auto shared_bases = std::make_shared<Barray<N>>(bases);

    auto pmap = SpPolicy::default_pmap(world, tvolume);
    for (auto const &ord : *pmap) {
        detail::IdxVec idx = trange.tiles().idx(ord);
        auto range = trange.make_tile_range(ord);
        world.taskq.add(detail::direct_tile_task<E, N, Op, ScreenOp>, ord,
                        range, idx, engines, shared_bases, op, &tile_norms,
                        &tiles);
    }
    world.gop.fence();

    SpShapeF shape(world, tile_norms, trange);
    DArray<N, DirectTile<E, N, Op>, SpPolicy> out(world, trange, shape, pmap);

    for (auto &&tile : tiles) {
        if (!out.is_zero(tile.first)) {
            out.set(tile.first, std::move(tile.second));
        }
    }


    return out;
}

} // namespace integrals
} // namespace mpqc

#endif // MPQC_INTEGRALS_DIRECTTASKINTEGRALS_H
