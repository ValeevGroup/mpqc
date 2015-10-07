#pragma once
#ifndef MPQC_INTEGRALS_DIRECTTASKINTEGRALS_H
#define MPQC_INTEGRALS_DIRECTTASKINTEGRALS_H

#include "direct_tile.h"
#include "task_integrals_op_invoker.h"
#include "task_integrals_common.h"

namespace mpqc {
namespace integrals {

namespace detail {} // namespace detail

/*! \brief Construct direct integral tensors in parallel without screening.
 *
 */
template <typename E, unsigned long N, typename Op>
DArray<N, DirectTile<E, N, Op, detail::no_screen_t<E, N>>, DnPolicy>
direct_task_ints_unscreeened(mad::World &world, ShrPool<E> &engines,
                             Barray<N> const &bases, Op op) {

    auto tile = DirectTile<E, N, Op, detail::no_screen_t<E,N>>();

    DArray<N, DirectTile<E, N, Op, detail::no_screen_t<E, N>>, DnPolicy> out(
          world, detail::create_trange(bases));

    auto shared_bases = std::make_shared<Barray<N>>(bases);

    auto const &pmap = *(out.get_pmap());
    auto const &trange = out.trange();
    for (auto const &ord : pmap) {
        detail::IdxVec idx = trange.tiles().idx(ord);
        auto range = trange.make_tile_range(ord);
        out.set(ord,
                make_direct_tile(std::move(range), std::move(idx), engines,
                                 shared_bases, op, detail::no_screening<E, N>));
    }

    return out;
}

// TODO FINISH THIS TOMORROW 
/*! \brief Construct direct integral tensors in parallel with screening.
 *
 */
template <typename E, unsigned long N, typename Op>
DArray<N, DirectTile<E, N, Op, detail::no_screen_t<E, N>>, DnPolicy>
direct_task_ints_screeened(mad::World &world, ShrPool<E> &engines,
                             Barray<N> const &bases, Op op) {

    auto tile = DirectTile<E, N, Op, detail::no_screen_t<E,N>>();
    const auto Q_X
          = screening_matrix_X(world, engines, bases[0].cluster_shells());
    auto const &cls_X = Q_X.cluster_screening;

    const auto Q_ab
          = screening_matrix_ab(world, engines, bases[1].cluster_shells());
    auto const &cls_ab = Q_ab.cluster_screening;

    DArray<N, DirectTile<E, N, Op, detail::no_screen_t<E, N>>, DnPolicy> out(
          world, detail::create_trange(bases));

    auto shared_bases = std::make_shared<Barray<N>>(bases);

    auto const &pmap = *(out.get_pmap());
    auto const &trange = out.trange();
    for (auto const &ord : pmap) {
        detail::IdxVec idx = trange.tiles().idx(ord);
        auto range = trange.make_tile_range(ord);
        out.set(ord,
                make_direct_tile(std::move(range), std::move(idx), engines,
                                 shared_bases, op, detail::no_screening<E, N>));
    }

    return out;
}

} // namespace integrals
} // namespace mpqc

#endif // MPQC_INTEGRALS_DIRECTTASKINTEGRALS_H
