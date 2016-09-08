#pragma once
#ifndef MPQC_INTEGRALS_DIRECTTASKINTEGRALS_H
#define MPQC_INTEGRALS_DIRECTTASKINTEGRALS_H

#include <mpqc/chemistry/qc/integrals/direct_tile.h>
#include <mpqc/chemistry/qc/integrals/integral_builder.h>
#include <mpqc/chemistry/qc/integrals/task_integrals_common.h>

#include "../../../../../ta_routines/tile_convert.h"
#include <limits>

namespace mpqc {
namespace integrals {

/*! \brief Construct direct integral tensors in parallel with screening.
 *
 * Same requirements on Op as those in Integral Builder
 */
template <typename Tile=TA::TensorD, typename Engine, unsigned long N>
DirectArray<Tile, TA::SparsePolicy, Engine> soad_direct_integrals(
    mad::World &world, ShrPool<Engine> shr_pool,
    Barray<N> const &bases,
    std::function<Tile(TA::TensorD &&)> op = mpqc::ta_routines::TensorDPassThrough()) {
  const auto trange = detail::create_trange(bases);
  const auto tvolume = trange.tiles().volume();
  TA::TensorF tile_norms(trange.tiles(), 0.0);

  // Copy the Bases for the Integral Builder
  auto shr_bases = std::make_shared<Barray<N>>(bases);

  // No screening
  auto screen = std::make_shared<Screener>();

  auto builder =
      make_integral_builder(world, std::move(shr_pool), std::move(shr_bases),
                            std::move(screen), std::move(op));

  auto dir_array = DirectArray<Tile, TA::SparsePolicy, Engine>(std::move(builder));
  auto builder_ptr = dir_array.builder();

  using DirectTileType = DirectTile<Tile,Engine>;
  std::vector<std::pair<int64_t, DirectTileType>> tiles(tvolume);
  

  auto task_f = [=](int64_t ord, detail::IdxVec const &idx, TA::Range rng,
                    TA::TensorF *tile_norms_ptr, DirectTileType *out_tile) {
    // No truncating for soad due to cost.
    const auto tile_norm = 10 * rng.volume();

    // Since no truncating we must be prepaired for every tile
    *out_tile = DirectTileType(idx, std::move(rng), std::move(builder_ptr));
    auto &norms = *tile_norms_ptr;
    norms[ord] = tile_norm;
  };

  auto pmap = SpPolicy::default_pmap(world, tvolume);
  for (auto const &ord : *pmap) {
    detail::IdxVec idx = trange.tiles().idx(ord);
    tiles[ord].first = ord;
    auto range = trange.make_tile_range(ord);
    world.taskq.add(task_f, ord, idx, range, &tile_norms, &tiles[ord].second);
  }
  world.gop.fence();

  SpShapeF shape(world, tile_norms, trange);
  TA::DistArray<DirectTileType, SpPolicy> out(world, trange, shape, pmap);

  for (auto it : *out.get_pmap()) {
    if (!out.is_zero(it)) {
      out.set(it, std::move(tiles[it].second));
    }
  }

  dir_array.set_array(std::move(out));

  return dir_array;
}

/*! \brief Construct direct integral tensors in parallel with screening.
 *
 * Same requirements on Op as those in Integral Builder
 */
template <typename Tile=TA::TensorD, typename Engine, unsigned long N>
DirectArray<Tile, TA::SparsePolicy, Engine> direct_sparse_integrals(
    mad::World &world, ShrPool<Engine> shr_pool, Barray<N> const &bases,
    std::shared_ptr<Screener> screen = std::make_shared<Screener>(Screener{}),
    std::function<Tile(TA::TensorD &&)> op = mpqc::ta_routines::TensorDPassThrough())

{

  const auto trange = detail::create_trange(bases);
  const auto tvolume = trange.tiles().volume();
  TA::TensorF tile_norms(trange.tiles(), 0.0);

  // Copy the Bases for the Integral Builder
  auto shr_bases = std::make_shared<Barray<N>>(bases);

  auto builder =
      make_integral_builder(world, std::move(shr_pool), std::move(shr_bases),
                            std::move(screen), std::move(op));

  auto dir_array = DirectArray<Tile, TA::SparsePolicy,Engine>(std::move(builder));
  auto builder_ptr = dir_array.builder();
  using DirectTileType = DirectTile<Tile,Engine>;

  std::vector<std::pair<int64_t, DirectTileType>> tiles(tvolume);

  auto task_f = [=](int64_t ord, detail::IdxVec const &idx, TA::Range rng,
                    TA::TensorF *tile_norms_ptr, DirectTileType *out_tile) {

    auto &builder = *builder_ptr;
    auto ta_tile = builder.integrals(idx, rng);

    const auto tile_volume = ta_tile.range().volume();
    const auto tile_norm = ta_tile.norm();

    // Keep tile if it was significant.
    bool save_norm = tile_norm >= tile_volume * SpShapeF::threshold();
    if (save_norm) {
      *out_tile = DirectTileType(idx, std::move(rng), std::move(builder_ptr));

      auto &norms = *tile_norms_ptr;
      norms[ord] = tile_norm;
    }
  };

  auto pmap = SpPolicy::default_pmap(world, tvolume);
  for (auto const &ord : *pmap) {
    detail::IdxVec idx = trange.tiles().idx(ord);
    tiles[ord].first = ord;
    auto range = trange.make_tile_range(ord);
    world.taskq.add(task_f, ord, idx, range, &tile_norms, &tiles[ord].second);
  }
  world.gop.fence();

  SpShapeF shape(world, tile_norms, trange);
  TA::DistArray<DirectTileType, SpPolicy> out(world, trange, shape, pmap);

  for (auto it : *out.get_pmap()) {
    if (!out.is_zero(it)) {
      out.set(it, std::move(tiles[it].second));
    }
  }

  dir_array.set_array(std::move(out));
  return dir_array;
}

/*! \brief Construct direct integral tensors in parallel with screening.
 *
 * Same requirements on Op as those in Integral Builder.
 *
 * I only plan to use this for CADF, no point in truncating tiles.
 */
template <typename Tile=TA::TensorD, typename Engine, unsigned long N>
DirectArray<Tile, TA::SparsePolicy, Engine> untruncated_direct_sparse_integrals(
    mad::World &world, ShrPool<Engine> shr_pool, Barray<N> const &bases,
    std::shared_ptr<Screener> screen = std::make_shared<Screener>(Screener{}),
    std::function<Tile(TA::TensorD&&)> op = mpqc::ta_routines::TensorDPassThrough()) {
  const auto trange = detail::create_trange(bases);
  const auto tvolume = trange.tiles().volume();
  TA::TensorF tile_norms(trange.tiles(), 0.0);

  // Copy the Bases for the Integral Builder
  auto shr_bases = std::make_shared<Barray<N>>(bases);

  auto builder =
      make_integral_builder(world, std::move(shr_pool), std::move(shr_bases),
                            std::move(screen), std::move(op));

  auto dir_array = DirectArray<Tile, TA::SparsePolicy,Engine>(std::move(builder));
  auto builder_ptr = dir_array.builder();

  using DirectTileType = DirectTile<Tile,Engine>;

  std::vector<std::pair<int64_t, DirectTileType>> tiles(tvolume);

  auto task_f = [=](int64_t ord, detail::IdxVec const &idx, TA::Range rng,
                    TA::TensorF *tile_norms_ptr, DirectTileType *out_tile) {
    *out_tile = DirectTileType(idx, std::move(rng), std::move(builder_ptr));

    auto &norms = *tile_norms_ptr;
    norms[ord] = std::numeric_limits<float>::max();

  };

  auto pmap = SpPolicy::default_pmap(world, tvolume);
  for (auto const &ord : *pmap) {
    detail::IdxVec idx = trange.tiles().idx(ord);
    tiles[ord].first = ord;
    auto range = trange.make_tile_range(ord);
    world.taskq.add(task_f, ord, idx, range, &tile_norms, &tiles[ord].second);
  }
  world.gop.fence();

  SpShapeF shape(world, tile_norms, trange);
  TA::DistArray<DirectTileType, SpPolicy> out(world, trange, shape, pmap);

  for (auto it : *out.get_pmap()) {
    if (!out.is_zero(it)) {
      out.set(it, std::move(tiles[it].second));
    }
  }

  dir_array.set_array(std::move(out));
  return dir_array;
}

/*! \brief Construct direct dense integral tensors in parallel with screening.
 *
 * Same requirements on Op as those in Integral Builder
 */
template <typename Tile=TA::TensorD, typename Engine, unsigned long N>
DirectArray<Tile, TA::DensePolicy, Engine> direct_dense_integrals(
    mad::World &world, ShrPool<Engine> shr_pool, Barray<N> const &bases,
    std::shared_ptr<Screener> screen = std::make_shared<Screener>(Screener{}),
    std::function<Tile(TA::TensorD&&)> op = mpqc::ta_routines::TensorDPassThrough()) {
  const auto trange = detail::create_trange(bases);

  // Copy the Bases for the Integral Builder
  auto shr_bases = std::make_shared<Barray<N>>(bases);

  // Make a pointer to an Integral builder.
  auto builder =
      make_integral_builder(world, std::move(shr_pool), std::move(shr_bases),
                            std::move(screen), std::move(op));

  auto dir_array = DirectArray<Tile, TA::DensePolicy,Engine>(std::move(builder));
  auto builder_ptr = dir_array.builder();

  using DirectTileType = DirectTile<Tile,Engine>;
  
  TA::DistArray<DirectTileType, TA::DensePolicy> out(world, trange);

  auto pmap = out.get_pmap();
  for (auto const &ord : *pmap) {
    detail::IdxVec idx = trange.tiles().idx(ord);
    auto range = trange.make_tile_range(ord);
    out.set(ord, DirectTileType(idx, std::move(range), builder_ptr));
  }
  world.gop.fence();

  dir_array.set_array(std::move(out));

  return dir_array;
}

}  // namespace integrals
}  // namespace mpqc

#endif  // MPQC_INTEGRALS_DIRECTTASKINTEGRALS_H
