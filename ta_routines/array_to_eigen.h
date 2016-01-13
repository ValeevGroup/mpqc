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

/*! \brief copies the tiles fromn TA to eigen using tasks
 *
 * Takes tiles by copy since they are shallow copy and will not be written to.
 */
template <typename Tile>
void write_to_eigen_task(Tile t, Matrix<double> *mat) {
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

        auto const &mat_extent = A.trange().elements().extent();
        Matrix<T> out_mat = Matrix<T>::Zero(mat_extent[0], mat_extent[1]);

        // Copy A and make it replicated.  Making A replicated is a mutating op.
        auto repl_A = A;
        repl_A.make_replicated();

        // Loop over the array and assign the tiles to blocks of the Eigen Mat.
        auto pmap = repl_A.get_pmap();
        const auto end = pmap->end();
        for (auto it = pmap->begin(); it != end; ++it) {
            if (!repl_A.is_zero(*it)) {
                auto tile = repl_A.find(*it).get();
                A.get_world().taskq.add(write_to_eigen_task<TA::Tensor<T>>, tile, &out_mat);
            }
        }
        A.get_world().gop.fence(); // Can't let M go out of scope

        return out_mat;
    }

/*! \brief converts a New TiledArray::Array to an Eigen Matrix
 *
 * The interface for TiledArray::Arrays changed so we need a new function.
 *
 */
template <typename Tile, typename Policy>
MatrixD array_to_eigen(TA::DistArray<Tile, Policy> const &A) {

    auto trange = A.trange();
    assert(trange.tiles().rank() == 2);

    auto const &mat_extent = trange.elements().extent();
    MatrixD out_mat(mat_extent[0], mat_extent[1]);
    out_mat.setZero();

    // Copy A and make it replicated.  Making A replicated is a mutating op.
    auto repl_A = A;
    repl_A.make_replicated();

    // Loop over the array and assign the tiles to blocks of the Eigen Mat.
    auto pmap = repl_A.get_pmap();
    const auto end = pmap->end();
    for (auto it = pmap->begin(); it != end; ++it) {
        if (!repl_A.is_zero(*it)) {
            auto tile = repl_A.find(*it).get();
            A.get_world().taskq.add(write_to_eigen_task<Tile>, tile, &out_mat);
        }
    }
    A.get_world().gop.fence(); // Can't let M go out of scope

    return out_mat;
}

/*! \brief takes an Eigen matrix and converts it to the type of the template.
 *
 * The idea is that users will provide a specialization which converts to the
 * tile type that they want.
 */
template <typename TileType>
TileType mat_to_tile(TA::Range range, Matrix<double> const *M, double cut);

template <>
tensor::Tile<tensor::DecomposedTensor<double>>
mat_to_tile<tensor::Tile<tensor::DecomposedTensor<double>>>(
      TA::Range range, Matrix<double> const *M, double cut) {
    auto const extent = range.extent();
    auto local_range = TA::Range{extent[0], extent[1]};
    auto tensor = tensor::DecomposedTensor<double>(cut, TA::Tensor<double>(
                                                              local_range));
    auto t_map = TA::eigen_map(tensor.tensor(0), extent[0], extent[1]);

    auto const start = range.lobound();
    t_map = M->block(start[0], start[1], extent[0], extent[1]);
    return tensor::Tile<tensor::DecomposedTensor<double>>(range,
                                                          std::move(tensor));
}

template <>
TA::Tensor<double>
mat_to_tile<TA::Tensor<double>>(TA::Range range, Matrix<double> const *M,
                                double) {
    const auto extent = range.extent();
    auto tensor = TA::Tensor<double>(range);
    auto t_map = TA::eigen_map(tensor, extent[0], extent[1]);

    auto const start = range.lobound();
    t_map = M->block(start[0], start[1], extent[0], extent[1]);
    return tensor;
}

// M must be replicated on all nodes.
template <typename Tile>
TA::Array<double, 2, Tile, TA::SparsePolicy>
eigen_to_array(madness::World &world, Matrix<double> const &M,
               TA::TiledRange1 tr0, TA::TiledRange1 tr1, double cut=1e-7) {
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
                  = world.taskq.add(mat_to_tile<Tile>, range, &M, cut);
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
