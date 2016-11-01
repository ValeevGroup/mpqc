#pragma once
#ifndef MPQC_TAROUTINES_DIAGONALARRAY_H
#define MPQC_TAROUTINES_DIAGONALARRAY_H

#include <tiledarray.h>

#include "../src/mpqc/math/tensor/clr/decomposed_tensor.h"
#include "../src/mpqc/math/tensor/clr/mpqc_tile.h"
#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/math/tensor/clr/array_to_eigen.h"

namespace mpqc {
namespace array_ops {

template <typename Tile>
TA::DistArray<Tile, TA::SparsePolicy> create_diagonal_array_from_eigen(madness::World& world,
    const TA::TiledRange1 &trange1, const TA::TiledRange1 &trange2, typename Tile::numeric_type val){

    using numeric_type = typename Tile::numeric_type;

    std::size_t x = trange1.elements().second;
    std::size_t y = trange2.elements().second;

    TA_ASSERT(x == y);

    auto diag = Eigen::DiagonalMatrix<numeric_type, Eigen::Dynamic>(x);
    diag.setIdentity();
    diag = val*diag;

//    RowMatrix<numeric_type> matrix = diag;
    return array_ops::eigen_to_array<Tile>(world, diag, trange1, trange2);

};

template <typename T>
void make_diagonal_tile(TiledArray::Tensor<T> &tile, T val) {
  auto const extent = tile.range().extent();
  auto map = TiledArray::eigen_map(tile, extent[0], extent[1]);
  for (auto i = 0ul; i < extent[0]; ++i) {
    map(i, i) = val;
  }
}

template <typename T>
void make_diagonal_tile(tensor::Tile<tensor::DecomposedTensor<T>> &tile,
                        T val) {
  assert(tile.tile().ndecomp() == 1);
  auto extent = tile.range().extent();
  auto local_range = TA::Range(extent[0], extent[1]);
  auto tensor = tensor::DecomposedTensor<T>(tile.tile().cut(),
                                            TA::Tensor<T>(local_range, 0.0));
  auto map = TiledArray::eigen_map(tensor.tensor(0), extent[0], extent[1]);
  for (auto i = 0ul; i < extent[0]; ++i) {
    map(i, i) = val;
  }
  tile.tile() = std::move(tensor);
}

template <typename Tile>
TiledArray::DistArray<Tile, TiledArray::SparsePolicy> create_diagonal_matrix(
    TiledArray::DistArray<Tile, TiledArray::SparsePolicy> const &model,
    double val) {
  using Array = TiledArray::DistArray<Tile, TiledArray::SparsePolicy>;

  TiledArray::Tensor<float> tile_shape(model.trange().tiles_range(), 0.0);

  auto pmap = model.pmap();

  auto pmap_end = pmap->end();
  for (auto it = pmap->begin(); it != pmap_end; ++it) {
    auto idx = model.trange().tiles_range().idx(*it);
    if (idx[0] == idx[1]) {
      tile_shape[*it] = val;
    }
  }

  TiledArray::SparseShape<float> shape(model.world(), tile_shape,
                                       model.trange());

  Array diag(model.world(), model.trange(), shape);

  auto const &trange = diag.trange();
  pmap = diag.pmap();
  auto end = pmap->end();
  for (auto it = pmap->begin(); it != end; ++it) {
    const auto ord = *it;

        auto idx = trange.tiles_range().idx(ord);
        auto diagonal_tile
              = std::all_of(idx.begin(), idx.end(),
        [&](typename Array::size_type const &x) { return x == idx.front(); });

    using TileType = typename Array::value_type;
    if (diagonal_tile && !diag.is_zero(ord)) {
      TileType tile = TileType{trange.make_tile_range(ord), 0.0};
      make_diagonal_tile(tile, val);
      diag.set(ord, std::move(tile));
    }
  }

  return diag;
}
/*!
 * \brief takes a TiledArray::TiledRange and a value and returns a diagonal
 *matrix.
 *
 * This function creates a diagonal matrix given a TiledArray::TiledRange and
 * a value.  The template parameters must be a numeric type followed by a
 * tile type. Finally there must be a create_diagonal_tile overload for the tile
 * type that gets passed in.
 *
 * By default this function will create the array using the default
 * madness::World, but you can pass in a world as an optional third parameter.
 *
 * \todo Finish diagonal matrix this is a little trickier than previous identity
 *functions because it needs to gracefully handle non symmetric tiling.
 */
template <typename Tile>
TiledArray::DistArray<Tile, TiledArray::SparsePolicy> diagonal_matrix(
    TiledArray::TiledRange const &trange, double val,
    madness::World &world = madness::World::get_default()) {
  TA_ASSERT(trange.rank() == 2);

  using Array = TiledArray::DistArray<Tile, TiledArray::SparsePolicy>;

  TiledArray::Tensor<float> tile_norms(trange.tiles_range(), 0.0);

  TiledArray::SparseShape<float> shape(world, tile_norms, trange);

  Array diag(world, trange, shape);

  return diag;
}

template <typename T, unsigned int N, typename Tile>
TiledArray::Array<T, N, Tile, TiledArray::DensePolicy> create_diagonal_matrix(
    TiledArray::Array<T, N, Tile, TiledArray::DensePolicy> const &model,
    double val) {
  using Array = TiledArray::Array<T, N, Tile, TiledArray::DensePolicy>;

  Array diag(model.world(), model.trange());

  auto pmap_ptr = diag.pmap();
  const auto end = pmap_ptr->end();
  for (auto it = pmap_ptr->begin(); it != end; ++it) {
    const auto ord = *it;
    auto const &idx = diag.trange().tiles_range().idx(ord);
    auto tile = Tile{diag.trange().make_tile_range(ord)};
    auto const extent = tile.range().extent();
    auto map = TiledArray::eigen_map(tile, extent[0], extent[1]);

    if (idx[0] == idx[1]) {
      for (auto i = 0ul; i < extent[0]; ++i) {
        for (auto j = 0ul; j < extent[1]; ++j) {
          if (i != j) {
            map(i, j) = 0;
          } else {
            map(i, i) = val;
          }
        }
      }
    } else {
      for (auto i = 0ul; i < extent[0]; ++i) {
        for (auto j = 0ul; j < extent[1]; ++j) {
          map(i, j) = 0;
        }
      }
    }

    diag.set(ord, std::move(tile));
  }

  return diag;
}

}  // namespace array_ops
}  // namespace mpqc

#endif  // MPQC_TAROUTINES_DIAGONALARRAY_H
