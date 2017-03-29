
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_DIRECT_TASK_INTEGRALS_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_DIRECT_TASK_INTEGRALS_H_

#include <limits>

#include "mpqc/chemistry/qc/lcao/integrals/direct_tile.h"
#include "mpqc/chemistry/qc/lcao/integrals/integral_builder.h"
#include "mpqc/chemistry/qc/lcao/integrals/task_integrals_common.h"
#include "mpqc/util/misc/time.h"

namespace mpqc {
namespace lcao {
namespace gaussian {

/*! \brief Construct direct integral tensors in parallel with screening.
 *
 * Same requirements on Op as those in Integral Builder
 */
template <typename Tile = TA::TensorD, typename Engine>
DirectArray<Tile, TA::SparsePolicy, Engine> direct_sparse_integrals(
    madness::World &world, ShrPool<Engine> shr_pool, BasisVector const &bases,
    std::shared_ptr<Screener> screen = std::make_shared<Screener>(Screener{}),
    std::function<Tile(TA::TensorD &&)> op = TA::Noop<TA::TensorD, true>(),
    std::shared_ptr<const math::PetiteList> plist =
        math::PetiteList::make_trivial()) {
  const auto trange = detail::create_trange(bases);

  auto pmap =
      TA::SparsePolicy::default_pmap(world, trange.tiles_range().volume());

  TA::TensorF tile_norms = screen->norm_estimate(world, bases, *pmap);

  TA::SparseShape<float> shape(world, tile_norms, trange);

  // Copy the Bases for the Integral Builder
  auto shr_bases = std::make_shared<BasisVector>(bases);

  auto builder = make_direct_integral_builder(
      world, std::move(shr_pool), std::move(shr_bases), std::move(screen),
      std::move(op), std::move(plist));

  auto dir_array =
      DirectArray<Tile, TA::SparsePolicy, Engine>(std::move(builder));

  auto builder_ptr = dir_array.builder();

  using DirectTileType = DirectTile<Tile, Engine>;

  auto task_f = [=](int64_t ord, detail::IdxVec const &idx, TA::Range rng) {
    return DirectTileType(idx, std::move(rng), std::move(builder_ptr));
  };

  TA::DistArray<DirectTileType, TA::SparsePolicy> out(world, trange, shape,
                                                      pmap);

  for (auto const &ord : *pmap) {
    if (!out.is_zero(ord)) {
      detail::IdxVec idx = trange.tiles_range().idx(ord);
      auto range = trange.make_tile_range(ord);
      auto tile_fut = world.taskq.add(task_f, ord, idx, range);
      out.set(ord, tile_fut);
    }
  }
  world.gop.fence();

  dir_array.set_array(std::move(out));
  return dir_array;
}

/*! \brief Construct direct integral tensors in parallel with screening.
 *
 * \param user_provided_norms is a user supplied replicated tensor with the
 * norm estimates for the array. It is of type vector<std::pair<Idx, float>>
 * where Idx is a type that provides random access via operator[] and the
 * values are the indices of the tile. The float is then the norm estimate for
 * that tile.
 */
template <typename Tile = TA::TensorD, typename Engine, typename Idx>
DirectArray<Tile, TA::SparsePolicy, Engine> direct_sparse_integrals(
    madness::World &world, ShrPool<Engine> shr_pool, BasisVector const &bases,
    std::vector<std::pair<Idx, float>> const &user_provided_norms,
    std::shared_ptr<Screener> screen = std::make_shared<Screener>(Screener{}),
    std::function<Tile(TA::TensorD &&)> op = TA::Noop<TA::TensorD, true>(),
    std::shared_ptr<const math::PetiteList> plist =
        math::PetiteList::make_trivial()) {
  const auto trange = detail::create_trange(bases);
  auto pmap =
      TA::SparsePolicy::default_pmap(world, trange.tiles_range().volume());

  // Don't use world constructor since norms must be replicated
  TA::SparseShape<float> shape(user_provided_norms, trange);

  // Copy the Bases for the Integral Builder
  auto shr_bases = std::make_shared<BasisVector>(bases);

  auto builder = make_direct_integral_builder(
      world, std::move(shr_pool), std::move(shr_bases), std::move(screen),
      std::move(op), std::move(plist));

  auto dir_array =
      DirectArray<Tile, TA::SparsePolicy, Engine>(std::move(builder));

  auto builder_ptr = dir_array.builder();

  using DirectTileType = DirectTile<Tile, Engine>;

  auto task_f = [=](detail::IdxVec idx, TA::Range rng) {
    return DirectTileType(std::move(idx), std::move(rng),
                          std::move(builder_ptr));
  };

  TA::DistArray<DirectTileType, TA::SparsePolicy> out(world, trange, shape,
                                                      pmap);
  for (auto const &ord : *pmap) {
    if (!out.is_zero(ord)) {
      auto range = trange.make_tile_range(ord);
      detail::IdxVec idx = trange.tiles_range().idx(ord);
      auto tile_fut = world.taskq.add(task_f, idx, range);
      out.set(ord, tile_fut);
    }
  }
  world.gop.fence();

  dir_array.set_array(std::move(out));
  return dir_array;
}

/*! \brief Construct direct integral tensors in parallel with screening.
 *
 * Same requirements on Op as those in Integral Builder.
 *
 * I only plan to use this for CADF, no point in truncating tiles.
 */
template <typename Tile = TA::TensorD, typename Engine>
DirectArray<Tile, TA::SparsePolicy, Engine> untruncated_direct_sparse_integrals(
    madness::World &world, ShrPool<Engine> shr_pool, BasisVector const &bases,
    std::shared_ptr<Screener> screen = std::make_shared<Screener>(Screener{}),
    std::function<Tile(TA::TensorD &&)> op = TA::Noop<TA::TensorD, true>(),
    std::shared_ptr<const math::PetiteList> plist =
        math::PetiteList::make_trivial()) {
  const auto trange = detail::create_trange(bases);
  const auto tvolume = trange.tiles_range().volume();
  TA::TensorF tile_norms(trange.tiles_range(), 0.0);

  // Copy the Bases for the Integral Builder
  auto shr_bases = std::make_shared<BasisVector>(bases);

  auto builder = make_direct_integral_builder(
      world, std::move(shr_pool), std::move(shr_bases), std::move(screen),
      std::move(op), std::move(plist));

  auto dir_array =
      DirectArray<Tile, TA::SparsePolicy, Engine>(std::move(builder));
  auto builder_ptr = dir_array.builder();

  using DirectTileType = DirectTile<Tile, Engine>;

  std::vector<std::pair<int64_t, DirectTileType>> tiles(tvolume);

  auto task_f = [=](int64_t ord, detail::IdxVec const &idx, TA::Range rng,
                    TA::TensorF *tile_norms_ptr, DirectTileType *out_tile) {
    *out_tile = DirectTileType(idx, std::move(rng), std::move(builder_ptr));

    auto &norms = *tile_norms_ptr;
    norms[ord] = std::numeric_limits<float>::max();

  };

  auto pmap = TA::SparsePolicy::default_pmap(world, tvolume);
  for (auto const &ord : *pmap) {
    detail::IdxVec idx = trange.tiles_range().idx(ord);
    tiles[ord].first = ord;
    auto range = trange.make_tile_range(ord);
    world.taskq.add(task_f, ord, idx, range, &tile_norms, &tiles[ord].second);
  }
  world.gop.fence();

  TA::SparseShape<float> shape(world, tile_norms, trange);
  TA::DistArray<DirectTileType, TA::SparsePolicy> out(world, trange, shape,
                                                      pmap);

  for (auto it : *out.pmap()) {
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
template <typename Tile = TA::TensorD, typename Engine>
DirectArray<Tile, TA::DensePolicy, Engine> direct_dense_integrals(
    madness::World &world, ShrPool<Engine> shr_pool, BasisVector const &bases,
    std::shared_ptr<Screener> screen = std::make_shared<Screener>(Screener{}),
    std::function<Tile(TA::TensorD &&)> op = TA::Noop<TA::TensorD, true>(),
    std::shared_ptr<const math::PetiteList> plist =
        math::PetiteList::make_trivial()) {
  const auto trange = detail::create_trange(bases);

  // Copy the Bases for the Integral Builder
  auto shr_bases = std::make_shared<BasisVector>(bases);

  // Make a pointer to an Integral builder.
  auto builder = make_direct_integral_builder(
      world, std::move(shr_pool), std::move(shr_bases), std::move(screen),
      std::move(op), std::move(plist));

  auto dir_array =
      DirectArray<Tile, TA::DensePolicy, Engine>(std::move(builder));
  auto builder_ptr = dir_array.builder();

  using DirectTileType = DirectTile<Tile, Engine>;

  TA::DistArray<DirectTileType, TA::DensePolicy> out(world, trange);

  auto pmap = out.pmap();
  for (auto const &ord : *pmap) {
    detail::IdxVec idx = trange.tiles_range().idx(ord);
    auto range = trange.make_tile_range(ord);
    out.set(ord, DirectTileType(idx, std::move(range), builder_ptr));
  }
  world.gop.fence();

  dir_array.set_array(std::move(out));

  return dir_array;
}

/**
 * construct direct dense LCAO integral from density fitting
 */

template <typename Tile>
TA::DistArray<DirectDFTile<Tile, TA::DensePolicy>, TA::DensePolicy>
df_direct_integrals(TA::DistArray<Tile, TA::DensePolicy> &bra,
                          TA::DistArray<Tile, TA::DensePolicy> &ket) {
  std::vector<TA::TiledRange1> vec_tr1(4);
  vec_tr1[0] = bra.trange().data()[1];
  vec_tr1[1] = bra.trange().data()[2];
  vec_tr1[2] = ket.trange().data()[1];
  vec_tr1[3] = ket.trange().data()[2];

  TA::TiledRange trange(vec_tr1.begin(), vec_tr1.end());

  auto builder_ptr =
      std::make_shared<DirectDFIntegralBuilder<Tile, TA::DensePolicy>>(bra,
                                                                       ket);

  auto world = bra.world();

  TA::DistArray<DirectDFTile<Tile, TA::DensePolicy>, TA::DensePolicy> result(
      world, trange);

  auto task_f = [builder_ptr](detail::IdxVec idx, TA::Range &range) {
    return DirectDFTile<Tile, TA::DensePolicy>(std::move(idx), std::move(range),
                                               std::move(builder_ptr));
  };

  auto pmap = result.pmap();
  for (const auto &ord : *pmap) {
    auto idx = trange.tiles_range().idx(ord);
    auto range = trange.make_tile_range(ord);
    auto tile_future = world.taskq.add(task_f, idx, range);
    result.set(ord, tile_future);
  }

  world.gop.fence();

  return result;
};

/**
 * construct direct sparse LCAO integral from density fitting
 */
template <typename Tile>
TA::DistArray<DirectDFTile<Tile, TA::SparsePolicy>, TA::SparsePolicy>
df_direct_integrals(TA::DistArray<Tile, TA::SparsePolicy> &bra,
                          TA::DistArray<Tile, TA::SparsePolicy> &ket) {
  std::vector<TA::TiledRange1> vec_tr1(4);
  vec_tr1[0] = bra.trange().data()[1];
  vec_tr1[1] = bra.trange().data()[2];
  vec_tr1[2] = ket.trange().data()[1];
  vec_tr1[3] = ket.trange().data()[2];

  TA::TiledRange trange(vec_tr1.begin(), vec_tr1.end());

  const auto tvolume = trange.tiles_range().volume();
  TA::TensorF tile_norms(trange.tiles_range(), 0.0);

  auto builder_ptr =
      std::make_shared<DirectDFIntegralBuilder<Tile, TA::SparsePolicy>>(bra,
                                                                        ket);

  auto& world = bra.world();

  auto task_norm = [builder_ptr](int64_t ord, const detail::IdxVec &idx,
                                 const TA::Range &range,
                                 TA::TensorF &tile_norm) {
    auto tile = builder_ptr->operator()(idx, range);
    const auto tile_volume = tile.range().volume();
    const auto norm = tile.norm();

    bool save_norm =
        norm >= tile_volume * TA::SparseShape<float>::threshold();
    if (save_norm) {
      tile_norm[ord] = norm;
    }

  };

  // initialize norm by compute all tiles once
  auto pmap = TA::SparsePolicy::default_pmap(world, tvolume);
  for (const auto &ord : *pmap) {
    detail::IdxVec idx = trange.tiles_range().idx(ord);
    auto range = trange.make_tile_range(ord);
    world.taskq.add(task_norm, ord, idx, range, tile_norms);
  }
  world.gop.fence();

  TA::SparseShape<float> shape(world, tile_norms, trange);

  TA::DistArray<DirectDFTile<Tile, TA::SparsePolicy>, TA::SparsePolicy> result(
      world, trange, shape, pmap);

  auto task_f = [builder_ptr](detail::IdxVec idx, TA::Range &range) {
    return DirectDFTile<Tile, TA::SparsePolicy>(
        std::move(idx), std::move(range), std::move(builder_ptr));
  };

  // set all tiles
  for (const auto &ord : *pmap) {
    if (!result.is_zero(ord)) {
      auto idx = trange.tiles_range().idx(ord);
      auto range = trange.make_tile_range(ord);
      auto tile_future = world.taskq.add(task_f, idx, range);
      result.set(ord, tile_future);
    }
  }
  world.gop.fence();

  return result;
};

}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_DIRECT_TASK_INTEGRALS_H_
