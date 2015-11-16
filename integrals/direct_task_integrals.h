#pragma once
#ifndef MPQC_INTEGRALS_DIRECTTASKINTEGRALS_H
#define MPQC_INTEGRALS_DIRECTTASKINTEGRALS_H

#include "direct_tile.h"
#include "task_integrals_common.h"
#include "integral_builder.h"

#include <limits>

namespace mpqc {
namespace integrals {

/*! \brief Construct direct integral tensors in parallel with screening.
 *
 * Same requirements on Op as those in Integral Builder
 */
template <typename E, unsigned long N, typename Op = TensorPassThrough>
DArray<N, DirectTile<IntegralBuilder<N, E, Op>>, SpPolicy>
soad_direct_integrals(mad::World &world, ShrPool<E> shr_pool,
                      Barray<N> const &bases, Op op = Op{}) {

    const auto trange = detail::create_trange(bases);
    const auto tvolume = trange.tiles().volume();
    TA::TensorF tile_norms(trange.tiles(), 0.0);

    // Copy the Bases for the Integral Builder
    auto shr_bases = std::make_shared<Barray<N>>(bases);

    // No screening
    auto screen = std::make_shared<Screener>();

    // Make a pointer to an Integral builder.  Doing this because we want to use
    // it in Tasks.
    auto builder_ptr = std::make_shared<IntegralBuilder<N, E, Op>>(
          make_integral_builder(world, std::move(shr_pool),
                                std::move(shr_bases), std::move(screen),
                                std::move(op)));

    using TileType = DirectTile<IntegralBuilder<N, E, Op>>;
    std::vector<std::pair<int64_t, TileType>> tiles(tvolume);

    auto task_f = [=](int64_t ord, detail::IdxVec const &idx, TA::Range rng,
                      TA::TensorF *tile_norms_ptr, TileType *out_tile) {
        // No truncating for soad due to cost.
        const auto tile_norm = 10 * rng.volume();

        // Since no truncating we must be prepaired for every tile
        *out_tile = TileType(idx, std::move(rng), std::move(builder_ptr));
        auto &norms = *tile_norms_ptr;
        norms[ord] = tile_norm;
    };

    auto pmap = SpPolicy::default_pmap(world, tvolume);
    for (auto const &ord : *pmap) {
        detail::IdxVec idx = trange.tiles().idx(ord);
        tiles[ord].first = ord;
        auto range = trange.make_tile_range(ord);
        world.taskq.add(task_f, ord, idx, range, &tile_norms,
                        &tiles[ord].second);
    }
    world.gop.fence();

    SpShapeF shape(world, tile_norms, trange);
    DArray<N, TileType, SpPolicy> out(world, trange, shape, pmap);

    for (auto &&tile : tiles) {
        if (out.is_local(tile.first) && !out.is_zero(tile.first)) {
            out.set(tile.first, std::move(tile.second));
        }
    }

    return out;
}

/*! \brief Construct direct integral tensors in parallel with screening.
 *
 * Same requirements on Op as those in Integral Builder
 */
template <typename E, unsigned long N, typename Op = TensorPassThrough>
DArray<N, DirectTile<IntegralBuilder<N, E, Op>>, SpPolicy>
direct_sparse_integrals(mad::World &world, ShrPool<E> shr_pool,
                        Barray<N> const &bases,
                        std::shared_ptr<Screener> screen
                        = std::make_shared<Screener>(Screener{}),
                        Op op = Op{}) {

    const auto trange = detail::create_trange(bases);
    const auto tvolume = trange.tiles().volume();
    TA::TensorF tile_norms(trange.tiles(), 0.0);

    // Copy the Bases for the Integral Builder
    auto shr_bases = std::make_shared<Barray<N>>(bases);

    // Make a pointer to an Integral builder.  Doing this because we want to use
    // it in Tasks.
    auto builder_ptr = std::make_shared<IntegralBuilder<N, E, Op>>(
          make_integral_builder(world, std::move(shr_pool),
                                std::move(shr_bases), std::move(screen),
                                std::move(op)));

    using TileType = DirectTile<IntegralBuilder<N, E, Op>>;
    std::vector<std::pair<int64_t, TileType>> tiles(tvolume);

    auto task_f = [=](int64_t ord, detail::IdxVec const &idx, TA::Range rng,
                      TA::TensorF *tile_norms_ptr, TileType *out_tile) {

        auto &builder = *builder_ptr;
        auto ta_tile = builder.integrals(idx, rng);

        const auto tile_volume = ta_tile.range().volume();
        const auto tile_norm = ta_tile.norm();

        // Keep tile if it was significant.
        bool save_norm = tile_norm >= tile_volume * SpShapeF::threshold();
        if (save_norm) {
            *out_tile = TileType(idx, std::move(rng), std::move(builder_ptr));

            auto &norms = *tile_norms_ptr;
            norms[ord] = tile_norm;
        }
    };

    auto pmap = SpPolicy::default_pmap(world, tvolume);
    for (auto const &ord : *pmap) {
        detail::IdxVec idx = trange.tiles().idx(ord);
        tiles[ord].first = ord;
        auto range = trange.make_tile_range(ord);
        world.taskq.add(task_f, ord, idx, range, &tile_norms,
                        &tiles[ord].second);
    }
    world.gop.fence();

    SpShapeF shape(world, tile_norms, trange);
    DArray<N, TileType, SpPolicy> out(world, trange, shape, pmap);

    for (auto &&tile : tiles) {
        if (out.is_local(tile.first) && !out.is_zero(tile.first)) {
            out.set(tile.first, std::move(tile.second));
        }
    }

    return out;
}

} // namespace integrals
} // namespace mpqc

#endif // MPQC_INTEGRALS_DIRECTTASKINTEGRALS_H
