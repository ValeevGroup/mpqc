#pragma once
#ifndef TCC_ARRAYOPS_ARRAYTOEIGEN_H
#define TCC_ARRAYOPS_ARRAYTOEIGEN_H

#include "../include/tiledarray.h"
#include "../tensor/tcc_tile.h"
#include "../tensor/decomposed_tensor.h"
#include "../tensor/decomposed_tensor_nonintrusive_interface.h"
#include "../include/eigen.h"
#include "../common/typedefs.h"
#include "../common/namespaces.h"

namespace tcc {
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

/*! \brief converts a TiledArray::Array to an Eigen Matrix
 *
 * The function needs to know the tile type and policy type to work.
 *
 * \bug Sometimes this causes a an error when run with multiple mpi process
 */
template <typename T, typename Tile, typename Policy>
Matrix<T> array_to_eigen(TA::Array<T, 2, Tile, Policy> const &A) {

    auto const &mat_extent = A.trange().elements().extent();
    Matrix<T> out_mat = Matrix<T>::Zero(mat_extent[0], mat_extent[1]);

    // Copy A and make it replicated.  Making A replicated is a mutating op.
    auto repl_A = A;
    repl_A.make_replicated();

    // Loop over the array and assign the tiles to blocks of the Eigen Mat.
    auto pmap = repl_A.get_pmap();
    const auto end = pmap->end();
    for (auto it = pmap->begin(); it != end; ++it) {
        if (!A.is_zero(*it)) {
            auto tile = A.find(*it).get();

            auto const &start = tile.range().lobound();
            const auto extent = tile.range().extent();
            out_mat.block(start[0], start[1], extent[0], extent[1])
                  = tile_to_eigen(tile);
        }
    }

    // overflows for large arrays.
    // A.get_world().gop.sum(out_mat.data(), out_mat.size());
    return out_mat;
}

/*! \brief takes an Eigen matrix and converts it to the type of the template.
 *
 * The idea is that users will provide a specialization which converts to the
 * tile type that they want.
 */
template <typename TileType>
TileType mat_to_tile(TA::Range range, Matrix<double> const &M);

template <>
tensor::Tile<tensor::DecomposedTensor<double>>
mat_to_tile<tensor::Tile<tensor::DecomposedTensor<double>>>(
      TA::Range range, Matrix<double> const &M) {
    auto const extent = range.extent();
    auto local_range = TA::Range{extent[0], extent[1]};
    // TODO fix 1e-7 to use cut
    auto tensor = tensor::DecomposedTensor<double>(1e-7, TA::Tensor<double>(
                                                               local_range));
    auto t_map = TA::eigen_map(tensor.tensor(0), extent[0], extent[1]);

    auto const start = range.lobound();
    t_map = M.block(start[0], start[1], extent[0], extent[1]);
    return tensor::Tile<tensor::DecomposedTensor<double>>(range,
                                                          std::move(tensor));
}

template <>
TA::Tensor<double>
mat_to_tile<TA::Tensor<double>>(TA::Range range, Matrix<double> const &M) {
    const auto extent = range.extent();
    auto tensor = TA::Tensor<double>(range);
    auto t_map = TA::eigen_map(tensor, extent[0], extent[1]);

    auto const start = range.lobound();
    t_map = M.block(start[0], start[1], extent[0], extent[1]);
    return tensor;
}

// M must be replicated on all nodes.
template <typename Tile>
TA::Array<double, 2, Tile, TA::SparsePolicy>
eigen_to_array(madness::World &world, Matrix<double> const &M,
               TA::TiledRange1 tr0, TA::TiledRange1 tr1) {
    TA::TiledRange trange{tr0, tr1};
    TA::Tensor<float> norms(trange.tiles(), trange.elements().volume()
                                            / trange.tiles().volume());

    TA::SparseShape<float> shape(world, norms, trange);
    TA::Array<double, 2, Tile, TA::SparsePolicy> array(world, trange, shape);

    auto const &pmap = array.get_pmap();
    const auto end = pmap->end();
    for (auto it = pmap->begin(); it != end; ++it) {
        if (!array.is_zero(*it)) {
            auto range = trange.make_tile_range(*it);
            madness::Future<Tile> tile
                  = world.taskq.add(mat_to_tile<Tile>, range, M);
            array.set(*it, tile);
        }
    }

    world.gop.fence(); // Be careful because M is ref
    array.truncate();  // Be lazy and fix shape like this.
    return array;
}

} // namespace array_ops
} // namespace tcc


#endif // TCC_ARRAYOPS_ARRAYTOEIGEN_H
