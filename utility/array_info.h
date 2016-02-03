#pragma once
#ifndef MPQC_UTILITY_ARRAYINFO_H
#define MPQC_UTILITY_ARRAYINFO_H

#include "../common/typedefs.h"
#include "../include/tiledarray.h"

#include <array>

namespace mpqc {
namespace utility {

double tile_clr_storage(Tile<DecompTensorD> const &tile);
double tile_clr_storage(TA::TensorD const &tile);

template <typename TileType, typename Policy>
std::array<double, 3> array_storage(TA::DistArray<TileType, Policy> const &A) {
    std::array<double, 3> out = {{0.0, 0.0, 0.0}};
    double &full_size = out[0];
    double &sparse_size = out[1];
    double &low_size = out[2];

    auto const &pmap = A.get_pmap();
    TA::TiledRange const &trange = A.trange();
    const auto end = pmap->end();

    for (auto it = pmap->begin(); it != end; ++it) {
        const TA::Range range = trange.make_tile_range(*it);
        auto const size = range.volume();

        full_size += size;

        if (!A.is_zero(*it)) {
            sparse_size += size;
            low_size += tile_clr_storage(A.find(*it).get());
        }
    }

    A.get_world().gop.sum(&out[0], 3);

    out[0] *= 8 * 1e-9;
    out[1] *= 8 * 1e-9;
    out[2] *= 8 * 1e-9;

    return out;
}

} // namespace utility
} // namespace mpqc

#endif // MPQC_UTILITY_ARRAYINFO_H
