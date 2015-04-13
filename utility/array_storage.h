#pragma once
#ifndef TCC_UTILITY_ARRAYSTORAGE_H
#define TCC_UTILITY_ARRAYSTORAGE_H

#include "../common/typedefs.h"
#include "../include/tiledarray.h"
#include "../tensor/tile_pimpl.h"
#include "../tensor/tile_algebra.h"
#include "../tensor/tcc_tile.h"
#include "../tensor/decomposed_tensor.h"
#include "../tensor/decomposed_tensor_unary.h"
#include "../tensor/decomposed_tensor_subtraction.h"
#include "../include/eigen.h"

#include "parallel_print.h"

#include <numeric>
#include <string>

namespace tcc {
namespace utility {

inline double tile_size(tensor::TilePimpl<double> const &tile) {
    if (tile.isFull()) {
        return tile.tile().ftile().size();
    } else {
        return tile.tile().lrtile().size();
    }
}

inline double
tile_size(tensor::Tile<tensor::DecomposedTensor<double>> const &tile) {
    auto size = 0.0;
    for (auto const &t : tile.tile().tensors()) {
        size += t.range().volume();
    }
    return size;
}

inline double tile_size(TiledArray::Tensor<double> const &tile) { return 0.0; }

template <typename T, unsigned int DIM, typename TileType, typename Policy>
std::array<double, 3>
array_storage(TA::Array<T, DIM, TileType, Policy> const &A) {

    std::array<double, 3> out = {{0.0, 0.0, 0.0}};
    double &full_size = out[0];
    double &sparse_size = out[1];
    double &low_size = out[2];

    auto const &pmap = A.get_pmap();
    TA::TiledRange const &trange = A.trange();
    const auto end = pmap->end();
    for (auto it = pmap->begin(); it != end; ++it) {
        const TA::Range range = trange.make_tile_range(*it);
        auto const &size_array = range.size();
        auto const size = std::accumulate(size_array.begin(), size_array.end(),
                                          1, std::multiplies<unsigned long>{});
        full_size += size;

        if (!A.is_zero(*it)) {
            sparse_size += size;
            low_size += tile_size(A.find(*it).get());
        }
    }

    A.get_world().gop.sum(&out[0], 3);

    out[0] *= 8 * 1e-9;
    out[1] *= 8 * 1e-9;
    out[2] *= 8 * 1e-9;

    return out;
}

template <typename T, unsigned int DIM, typename TileType, typename Policy>
void print_size_info(TA::Array<T, DIM, TileType, Policy> const &a,
                     std::string name) {
    print_par(a.get_world(), "Printing size information for ", name, "\n");

    auto data = array_storage(a);

    print_par(a.get_world(), "\tFull   = ", data[0], " GB\n", "\tSparse = ",
              data[1], " GB\n", "\tLow Rank = ", data[2], " GB\n");
}

template <typename T, unsigned int DIM, typename TileTypeL, typename TileTypeR,
          typename PolicyL, typename PolicyR>
void print_array_difference(TA::Array<T, DIM, TileTypeL, PolicyL> const &left,
                            TA::Array<T, DIM, TileTypeR, PolicyR> const &right,
                            std::string const &left_name,
                            std::string const &right_name) {
    // This really needs to check that both pmaps are the same size.
    auto const &pmap_left = left.get_pmap();
    auto const &pmap_right = right.get_pmap();

    auto l_end = pmap_left->end();
    auto r_end = pmap_right->end();

    // Compute the norm difference for each tile.
    auto diff_norm_sum = 0.0;
    auto it_r = pmap_right->begin();
    for (auto it_l = pmap_left->begin(); it_l != l_end && it_r != r_end;
         ++it_l, ++it_r) {
        // Ensure we are looking at the same ordinal.
        assert(*it_l == *it_r);

        if (left.is_zero(*it_l)) {
            // do nothing if right is zero
            if (!right.is_zero(*it_r)) {
                diff_norm_sum += right.find(*it_r).get().norm();
            }
        } else {
            if (right.is_zero(*it_r)) {
                diff_norm_sum += left.find(*it_l).get().norm();
            } else {
                const auto l_t = left.find(*it_l).get();
                const auto r_t = right.find(*it_r).get();
                diff_norm_sum += l_t.subt(r_t).norm();
            }
        }
    }

    left.get_world().gop.sum(&diff_norm_sum, 1);

    print_par(left.get_world(), "The difference between array ", left_name,
              " and array ", right_name, " is ", std::sqrt(diff_norm_sum),
              "\n");
}

} // namespace utility
} // namespace tcc

#endif // TCC_UTILITY_ARRAYSTORAGE_H
