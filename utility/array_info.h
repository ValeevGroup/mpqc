#pragma once
#ifndef MPQC_UTILITY_ARRAYINFO_H
#define MPQC_UTILITY_ARRAYINFO_H

#include "../common/typedefs.h"
#include "../include/tiledarray.h"

#include <array>
#include <atomic>

namespace mpqc {
namespace utility {

unsigned long tile_clr_storage(Tile<DecompTensorD> const &tile);
unsigned long tile_clr_storage(TA::TensorD const &);

template <typename TileType, typename Policy>
std::array<double, 3> array_storage(TA::DistArray<TileType, Policy> const &A) {
    std::atomic_ulong full_size(0);
    std::atomic_ulong sparse_size(0);
    std::atomic_ulong low_size(0);

    TA::TiledRange const &trange = A.trange();

    auto task_is_zero = [&](unsigned long ord){
        full_size += trange.make_tile_range(ord).volume(); 
    };

    auto task_is_not_zero = [&](unsigned long ord, TileType const &t){
        const TA::Range range = trange.make_tile_range(ord);
        const auto size = range.volume();

        full_size += size;
        sparse_size += size;
        low_size += tile_clr_storage(t);
    };

    auto const &pmap = A.get_pmap();
    const auto end = pmap->end();

    for (auto it = pmap->begin(); it != end; ++it) {
        const auto ord = *it;
        if(A.is_local(ord)){
            if(!A.is_zero(ord)){
                A.get_world().taskq.add(task_is_not_zero, ord, A.find(ord));
            } else {
                A.get_world().taskq.add(task_is_zero, ord);
            }
        }
    }

    A.get_world().gop.fence();

    std::array<double, 3> out;
    out[0] = double(full_size);
    out[1] = double(sparse_size);
    out[2] = double(low_size);

    A.get_world().gop.sum(&out[0], 3);

    out[0] *= 8 * 1e-9;
    out[1] *= 8 * 1e-9;
    out[2] *= 8 * 1e-9;

    return out;
}

template <typename Tile>
double array_size(const TA::DistArray<Tile,TA::DensePolicy>& A){
    std::atomic_ulong full_size(0);
    auto const &pmap = A.get_pmap();
    TA::TiledRange const &trange = A.trange();

    auto task = [&trange, &full_size](unsigned long ord){
        const TA::Range range = trange.make_tile_range(ord);
        const auto size = range.volume();
        full_size += size;
    };


    const auto end = pmap->end();

    for (auto it = pmap->begin(); it != end; ++it) {
        A.get_world().taskq.add(task,*it);
    }
    A.get_world().gop.fence();

    double size = double(full_size)*8* 1e-9;
    A.get_world().gop.sum(size);
    return size;
}

template <typename Tile>
double array_size(const TA::DistArray<Tile,TA::SparsePolicy>& A){
    std::atomic_ulong full_size(0);
    auto const &pmap = A.get_pmap();
    TA::TiledRange const &trange = A.trange();

    auto task = [&trange, &full_size](unsigned long ord){
      const TA::Range range = trange.make_tile_range(ord);
      const auto size = range.volume();
      full_size += size;
    };

    const auto end = pmap->end();

    for (auto it = pmap->begin(); it != end; ++it) {
        if(!A.is_zero(*it)){
            A.get_world().taskq.add(task,*it);
        }
    }
    A.get_world().gop.fence();

    double size = double(full_size)*8* 1e-9;
    A.get_world().gop.sum(size);
    return size;
}

} // namespace utility
} // namespace mpqc

#endif // MPQC_UTILITY_ARRAYINFO_H
