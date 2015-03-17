#pragma once
#ifndef TCC_UTILITY_ARRAYSTORAGE_H
#define TCC_UTILITY_ARRAYSTORAGE_H

#include "../common/typedefs.h"
#include "../include/tiledarray.h"
#include "../tensor/tile_pimpl.h"
#include "../tensor/tile_algebra.h"
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
            low_size += tile_size(A.find(*it));
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
    print_par(a.get_world(), "Printing size information for ", name,
                       "\n");

    auto data = array_storage(a);

    print_par(a.get_world(), "\tFull   = ", data[0], " GB\n",
                       "\tSparse = ", data[1], " GB\n", "\tLow Rank = ",
                       data[2], " GB\n");
}

} // namespace utility
} // namespace tcc

#endif // TCC_UTILITY_ARRAYSTORAGE_H
