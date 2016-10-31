#pragma once
#ifndef MPQC_ARRAYOPS_ARRAYTOEIGEN_H
#define MPQC_ARRAYOPS_ARRAYTOEIGEN_H

#include "../common/namespaces.h"
#include "../common/typedefs.h"
#include <tiledarray.h>

#include "../src/mpqc/math/tensor/clr/decomposed_tensor.h"
#include "../src/mpqc/math/tensor/clr/decomposed_tensor_nonintrusive_interface.h"
#include "../src/mpqc/math/tensor/clr/mpqc_tile.h"
#include "mpqc/math/external/eigen/eigen.h"

namespace mpqc {
namespace array_ops {

template <typename T>
using Matrix = Eig::Matrix<T, Eig::Dynamic, Eig::Dynamic, Eig::RowMajor>;

template <typename T>
Matrix<T> tile_to_eigen(TA::Tensor<T> const &t) {
  auto const extent = t.range().extent();
  return TA::eigen_map(t, extent[0], extent[1]);
}

template <typename T>
Matrix<T> tile_to_eigen(tensor::Tile<tensor::DecomposedTensor<T>> const &t) {
  return tile_to_eigen(tensor::algebra::combine(t.tile()));
}

/*! \brief copies the tiles fromn TA to eigen using tasks
 *
 * Takes tiles by copy since they are shallow copy and will not be written to.
 */
template <typename Tile>
void write_to_eigen_task(Tile t, Matrix<typename Tile::numeric_type> *mat) {
  auto const &start = t.range().lobound();
  const auto extent = t.range().extent();
  mat->block(start[0], start[1], extent[0], extent[1]) = tile_to_eigen(t);
}

/*! \brief converts a TiledArray::Array to an Eigen Matrix
 *
 * The function needs to know the tile type and policy type to work.
 *
 */

template <typename T, typename Policy>
Matrix<T> array_to_eigen(TA::DistArray<TA::Tensor<T>, Policy> const &A) {
  TA_ASSERT(A.range().rank() == 2);

  auto const &mat_extent = A.trange().elements_range().extent();
  Matrix<T> out_mat = Matrix<T>::Zero(mat_extent[0], mat_extent[1]);

  // Copy A and make it replicated.  Making A replicated is a mutating op.
  auto repl_A = A;
  A.world().gop.fence();
  repl_A.make_replicated();

  // Loop over the array and assign the tiles to blocks of the Eigen Mat.
  auto pmap = repl_A.pmap();
  const auto end = pmap->end();
  for (auto it = pmap->begin(); it != end; ++it) {
    if (!repl_A.is_zero(*it)) {
      auto tile = repl_A.find(*it).get();
      A.world().taskq.add(write_to_eigen_task<TA::Tensor<T>>, tile,
                              &out_mat);
    }
  }
  A.world().gop.fence();  // Can't let M go out of scope

  return out_mat;
}

/*! \brief takes an Eigen matrix and converts it to the type of the template.
 *
 * The idea is that users will provide a specialization which converts to the
 * tile type that they want.
 */
template <typename TileType>
TileType mat_to_tile(TA::Range range, Matrix<typename TileType::numeric_type> const *M, double cut);

template <>
inline tensor::Tile<tensor::DecomposedTensor<double>>
mat_to_tile<tensor::Tile<tensor::DecomposedTensor<double>>>(
    TA::Range range, Matrix<double> const *M, double cut) {
  auto const extent = range.extent();
  auto local_range = TA::Range{extent[0], extent[1]};
  auto tensor =
      tensor::DecomposedTensor<double>(cut, TA::Tensor<double>(local_range));
  auto t_map = TA::eigen_map(tensor.tensor(0), extent[0], extent[1]);

  auto const start = range.lobound();
  t_map = M->block(start[0], start[1], extent[0], extent[1]);
  return tensor::Tile<tensor::DecomposedTensor<double>>(range,
                                                        std::move(tensor));
}

template <typename T>
inline TA::Tensor<T> mat_to_ta_tensor(TA::Range & range,
                                      Matrix<T> const *M)
{
  const auto extent = range.extent();
  auto tensor = TA::Tensor<T>(range);
  auto t_map = TA::eigen_map(tensor, extent[0], extent[1]);

  auto const start = range.lobound();
  t_map = M->block(start[0], start[1], extent[0], extent[1]);
  return tensor;
}

template <>
inline TA::TensorD mat_to_tile<TA::TensorD>(TA::Range range,
                                            Matrix<double> const *M, double) {
  return mat_to_ta_tensor(range,M);
}

template <>
inline TA::TensorZ mat_to_tile<TA::TensorZ>(TA::Range range,
                                            Matrix<std::complex<double>> const *M, double) {
  return mat_to_ta_tensor(range, M);
}

// M must be replicated on all nodes.
template <typename Tile>
TA::DistArray<Tile, TA::SparsePolicy> eigen_to_array(
    madness::World &world, Matrix<typename Tile::numeric_type> const &M, TA::TiledRange1 tr0,
    TA::TiledRange1 tr1, double cut = 1e-7) {
  TA::TiledRange trange{tr0, tr1};
  TA::Tensor<float> norms(trange.tiles_range(),
                          trange.elements_range().volume() / trange.tiles_range().volume());

  TA::SparseShape<float> shape(world, norms, trange);
  TA::DistArray<Tile, TA::SparsePolicy> array(world, trange, shape);

  auto const &pmap = array.pmap();
  const auto end = pmap->end();
  for (auto it = pmap->begin(); it != end; ++it) {
    if (!array.is_zero(*it)) {
      auto range = trange.make_tile_range(*it);
      madness::Future<Tile> tile =
          world.taskq.add(mat_to_tile<Tile>, range, &M, cut);
      array.set(*it, tile);
    }
  }

  world.gop.fence();  // Be careful because M is ref
  array.truncate();   // Be lazy and fix shape like this.
  return array;
}

}  // namespace array_ops
}  // namespace mpqc

#endif  // mpqc_ARRAYOPS_ARRAYTOEIGEN_H
