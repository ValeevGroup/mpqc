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

#include <mpqc/chemistry/qc/integrals/task_integrals_common.h>
#include <mpqc/chemistry/qc/integrals/screening/screen_base.h>
#include <mpqc/chemistry/qc/integrals/integral_builder.h>

namespace mpqc {
namespace integrals {

/*! \brief Construct sparse integral tensors from sets in parallel.
 *
 * This is needed for integrals such as the dipole integrals that come as a
 * set.
 *
 * \param shr_pool should be a std::shared_ptr to an IntegralEnginePool
 * \param bases should be a std::array of Basis, which will be copied.
 * \param op needs to be a function or functor that takes a TA::TensorD && and
 * returns any valid tile type. Op is copied so it can be moved.
 * ```
 * auto t = [](TA::TensorD &&ten){return std::move(ten);};
 * ```
 */
template <typename E, typename Op = TensorPassThrough>
std::vector<DArray<2, detail::Ttype<Op>, SpPolicy>>
sparse_xyz_integrals(mad::World &world, ShrPool<E> shr_pool,
                     Barray<2> const &bases, Op op = Op{}) {

    using Tile = detail::Ttype<Op>;

    // Build the Trange and Shape Tensor
    auto trange = detail::create_trange(bases);
    const auto tvolume = trange.tiles().volume();

    using TileVec = std::vector<std::pair<unsigned long, Tile>>;
    std::vector<TileVec> tiles(3, TileVec(tvolume));

    using NormVec = std::vector<TA::TensorF>;
    NormVec tile_norms(3, TA::TensorF(trange.tiles(), 0.0));

    // Capture by ref since we are going to fence after loops.
    auto task_f = [&](int64_t ord, detail::IdxVec idx, TA::Range rng) {
        auto &eng = shr_pool->local();

        const auto idx0 = idx[0];
        const auto idx1 = idx[1];

        const auto t_volume = rng.volume();

        auto const &lobound = rng.lobound();
        std::array<std::size_t, 2> lb = {{lobound[0], lobound[1]}};
        std::array<std::size_t, 2> ub = {{lobound[0], lobound[1]}};

        std::vector<TA::TensorD> t_xyz
              = {TA::TensorD(rng, 0.0), TA::TensorD(rng, 0.0),
                 TA::TensorD(std::move(rng), 0.0)};

        const double dummy = 0.0;
        auto map = TA::make_map(&dummy, {0, 0}, {1, 1});

        auto const &shells0 = bases[0].cluster_shells()[idx0];
        auto const &shells1 = bases[1].cluster_shells()[idx1];

        for (auto const &s0 : shells0) {
            const auto n0 = s0.size();
            ub[0] += n0;

            lb[1] = ub[1] = lobound[1];
            for (auto const &s1 : shells1) {
                const auto n1 = s1.size();
                ub[1] += n1;

                const auto n1n2 = n0 * n1;
                const auto& bufs = eng.compute(s0, s1);
                TA_USER_ASSERT(bufs.size() >= 4,
                               "unexpected result from Engine::compute()");

                for (auto i = 1; i < 4; ++i) {
                    TA::remap(map, bufs[i], lb, ub);
                    t_xyz[i-1].block(lb, ub) = map;
                }

                lb[1] = ub[1];
            }
            lb[0] = ub[0];
        }

        std::array<double, 3> norm
              = {{t_xyz[0].norm(), t_xyz[1].norm(), t_xyz[2].norm()}};

        for (auto i = 0; i < 3; ++i) {
            tile_norms[i][ord] = norm[i];
            if (TA::SparseShape<float>::threshold() <= t_volume * norm[i]) {
                tiles[i][ord].second = op(std::move(t_xyz[i]));
            }
        }
    };

    auto pmap = SpPolicy::default_pmap(world, tvolume);
    for (auto const ord : *pmap) {
        tiles[0][ord].first = ord;
        tiles[1][ord].first = ord;
        tiles[2][ord].first = ord;
        detail::IdxVec idx = trange.tiles().idx(ord);
        world.taskq.add(task_f, ord, idx, trange.make_tile_range(ord));
    }
    world.gop.fence();

    std::vector<DArray<2, Tile, SpPolicy>> arrays(3);

    for (auto i = 0; i < 3; ++i) {
        TA::SparseShape<float> shape(world, tile_norms[i], trange);
        arrays[i] = DArray<2, Tile, SpPolicy>(world, trange, shape, pmap);
        detail::set_array(tiles[i], arrays[i]);
    }
    world.gop.fence();

    return arrays;
}

/*! \brief Construct sparse integral tensors in parallel.
 *
 * \param shr_pool should be a std::shared_ptr to an IntegralEnginePool
 * \param bases should be a std::array of Basis, which will be copied.
 * \param op needs to be a function or functor that takes a TA::TensorD && and
 * returns any valid tile type. Op is copied so it can be moved.
 * ```
 * auto t = [](TA::TensorD &&ten){return std::move(ten);};
 * ```
 *
 * \param screen should be a std::shared_ptr to a Screener.
 */
template <typename E, unsigned long N, typename Op = TensorPassThrough>
DArray<N, detail::Ttype<Op>, SpPolicy>
sparse_integrals(mad::World &world, ShrPool<E> shr_pool, Barray<N> const &bases,
                 std::shared_ptr<Screener> screen
                 = std::make_shared<Screener>(Screener{}),
                 Op op = Op{}) {

    using Tile = detail::Ttype<Op>;

    // Build the Trange and Shape Tensor
    auto trange = detail::create_trange(bases);
    const auto tvolume = trange.tiles().volume();
    std::vector<std::pair<unsigned long, Tile>> tiles(tvolume);
    TA::TensorF tile_norms(trange.tiles(), 0.0);
    
    // Copy the Bases for the Integral Builder
    auto shr_bases = std::make_shared<Barray<N>>(bases);

    // Make a pointer to an Integral builder.  Doing this because we want to use
    // it in Tasks.
    auto builder_ptr = 
          make_integral_builder(world, std::move(shr_pool),
                                std::move(shr_bases), std::move(screen),
                                std::move(op));

    auto task_f = [=](int64_t ord, detail::IdxVec idx, TA::Range rng,
                      TA::TensorF *tile_norms_ptr, Tile *out_tile) {

        // This is why builder was made into a shared_ptr.
        auto &builder = *builder_ptr;
        auto ta_tile = builder.integrals(idx, std::move(rng));

        const auto tile_volume = ta_tile.range().volume();
        const auto tile_norm = ta_tile.norm();

        //BUG cc-pv5z-ri doesn't work trying to figure out why.
        // if(tile_norm != tile_norm){
        //     std::cout << "Tile norm = " << tile_norm << std::endl;
        //     std::cout << "Tile norm broken idx = " << idx[0] << " " << idx[1] << " " << idx[2] << std::endl;
        //     assert(false);
        // }

        // Keep tile if it was significant.
        bool save_norm = tile_norm >= tile_volume * SpShapeF::threshold();
        if (save_norm) {
            *out_tile = builder.op(std::move(ta_tile));

            auto &norms = *tile_norms_ptr;
            norms[ord] = tile_norm;
        }
    };

    auto pmap = SpPolicy::default_pmap(world, tvolume);
    for (auto const ord : *pmap) {
        tiles[ord].first = ord;
        detail::IdxVec idx = trange.tiles().idx(ord);
        world.taskq.add(task_f, ord, idx, trange.make_tile_range(ord),
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
template <typename E, unsigned long N, typename Op = TensorPassThrough>
DArray<N, detail::Ttype<Op>, DnPolicy>
dense_integrals(mad::World &world, ShrPool<E> shr_pool, Barray<N> const &bases,
                std::shared_ptr<Screener> screen
                = std::make_shared<Screener>(Screener{}),
                Op op = Op{}) {

    using Tile = detail::Ttype<Op>;
    DArray<N, Tile, DnPolicy> out(world, detail::create_trange(bases));

    // Copy the Bases for the Integral Builder
    auto shr_bases = std::make_shared<Barray<N>>(bases);

    // Make a pointer to a builder which can be shared by different tasks.
    auto builder_ptr = make_integral_builder(world, std::move(shr_pool),
                                std::move(shr_bases), std::move(screen),
                                std::move(op));

    // builder is shared_ptr so just capture it by copy.
    auto task_func = [=](detail::IdxVec const &idx, TA::Range rng) {
        return builder_ptr->operator()(idx, std::move(rng));
    };

    auto const &trange = out.trange();
    auto const &pmap = *(out.get_pmap());
    for (auto const ord : pmap) {
        detail::IdxVec idx = trange.tiles().idx(ord);

        auto range = trange.make_tile_range(ord);
        mad::Future<Tile> tile = world.taskq.add(task_func, idx, range);

        out.set(ord, tile);
    }
    world.gop.fence();

    return out;
}

} // namespace integrals
} // namespace mpqc


#endif // MPQC_INTEGRALS_TASKINTEGRALS_H
