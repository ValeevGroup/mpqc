
#ifndef MPQC4_SRC_MPQC_MATH_TENSOR_CLR_ARRAY_TO_EIGEN_H_
#define MPQC4_SRC_MPQC_MATH_TENSOR_CLR_ARRAY_TO_EIGEN_H_

#include <tiledarray.h>

#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/math/tensor/clr/decomposed_tensor.h"
#include "mpqc/math/tensor/clr/decomposed_tensor_nonintrusive_interface.h"
#include "tile.h"

namespace mpqc {
namespace math {

template <typename T>
using Matrix =
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template <typename T>
Matrix<T> tile_to_eigen(TA::Tensor<T> const &t) {
  auto const extent = t.range().extent();
  MPQC_ASSERT(extent.size() == 2);
  return TA::eigen_map(t, extent[0], extent[1]);
}

template <typename T>
Matrix<T> tile_to_eigen(Tile<DecomposedTensor<T>> const &t) {
  return tile_to_eigen(combine(t.tile()));
}

/*! \brief converts a TiledArray::Array to a row-major Eigen Matrix
 *  \param[in] A an order-2 tensor
 *  \warning this is a collective operation over the World object returned by \c A.world()
 */
template <typename Tile, typename Policy>
Matrix<typename Tile::value_type>
    array_to_eigen(TA::DistArray<Tile, Policy> const &A) {
  // Copy A and make it replicated.  Making A replicated is a mutating op.
  auto A_repl = A;
  A_repl.make_replicated();
  return TA::array_to_eigen<Tile, Policy, Eigen::RowMajor>(A_repl);
}

/*! \brief takes an Eigen matrix and converts it to the type of the template.
 *
 * The idea is that users will provide a specialization which converts to the
 * tile type that they want.
 */
template <typename TileType>
TileType mat_to_tile(TA::Range range,
                     Matrix<typename TileType::numeric_type> const *M,
                     double cut);

template <>
inline Tile<DecomposedTensor<double>>
mat_to_tile<Tile<DecomposedTensor<double>>>(
    TA::Range range, Matrix<double> const *M, double cut) {
  auto const extent = range.extent();
  auto local_range = TA::Range{extent[0], extent[1]};
  auto tensor =
      DecomposedTensor<double>(cut, TA::Tensor<double>(local_range));
  auto t_map = TA::eigen_map(tensor.tensor(0), extent[0], extent[1]);

  auto const start = range.lobound();
  t_map = M->block(start[0], start[1], extent[0], extent[1]);
  return Tile<DecomposedTensor<double>>(range,
                                                        std::move(tensor));
}

template <typename T>
inline TA::Tensor<T> mat_to_ta_tensor(TA::Range &range, Matrix<T> const *M) {
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
  return mat_to_ta_tensor(range, M);
}

template <>
inline TA::TensorZ mat_to_tile<TA::TensorZ>(
    TA::Range range, Matrix<std::complex<double>> const *M, double) {
  return mat_to_ta_tensor(range, M);
}


/**
 *  convert eigen(replicated) to sparse TA::DistArray
 * @tparam Tile
 * @tparam Policy
 * @param world
 * @param M  Eigen matrix, must be replicated on all nodes
 * @param tr0 the TiledArray::TiledRange1 object for the row dimension of the result
 * @param tr1 the TiledArray::TiledRange1 object for the column dimension of the result
 * @param cut
 * @return
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, typename std::enable_if<std::is_same<Policy, TA::SparsePolicy>::value, TA::SparsePolicy>::type >
eigen_to_array(
    madness::World &world, Matrix<typename Tile::numeric_type> const &M,
    TA::TiledRange1 tr0, TA::TiledRange1 tr1, double cut = 1e-7) {

  // make sure dimensions of M match the dimensions of tiled ranges
  assert(std::size_t(M.rows()) == tr0.extent() && "eigen_to_array(): row dimensions do not match");
  assert(std::size_t(M.cols()) == tr1.extent() && "eigen_to_array(): col dimensions do not match");

  TA::TiledRange trange{tr0, tr1};
  TA::Tensor<float> norms(
      trange.tiles_range(),
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


/**
 *  convert eigen(replicated) to dense TA::DistArray
 * @tparam Tile
 * @tparam Policy TA::DensePolicy
 * @param world
 * @param M  Eigen matrix, must be replicated on all nodes
 * @param tr0 the TiledArray::TiledRange1 object for the row dimension of the result
 * @param tr1 the TiledArray::TiledRange1 object for the column dimension of the result
 * @param cut
 * @return
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, typename std::enable_if<std::is_same<Policy, TA::DensePolicy>::value, TA::DensePolicy>::type >
eigen_to_array(
    madness::World &world, Matrix<typename Tile::numeric_type> const &M,
    TA::TiledRange1 tr0, TA::TiledRange1 tr1, double cut = 1e-7) {

  // make sure dimensions of M match the dimensions of tiled ranges
  assert(M.rows() == tr0.extent() && "eigen_to_array(): row dimensions do not match");
  assert(M.cols() == tr1.extent() && "eigen_to_array(): col dimensions do not match");

  TA::TiledRange trange{tr0, tr1};
  TA::DistArray<Tile, TA::DensePolicy> array(world, trange);

  auto const &pmap = array.pmap();
  const auto end = pmap->end();
  for (auto it = pmap->begin(); it != end; ++it) {
    auto range = trange.make_tile_range(*it);
    madness::Future<Tile> tile =
        world.taskq.add(mat_to_tile<Tile>, range, &M, cut);
    array.set(*it, tile);
  }

  world.gop.fence();  // Be careful because M is ref
  return array;
}


}  // namespace math
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_MATH_TENSOR_CLR_ARRAY_TO_EIGEN_H_
