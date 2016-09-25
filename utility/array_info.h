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


/// @brief reports array size (in bytes) on each process
/// @note this is a collective operation
/// @return vector of integers, each element is the aggregate size of array's
///         tiles on the corresponding process
template <typename Tile, typename Policy,
          typename = typename std::enable_if<
              std::is_fundamental<typename Tile::value_type>::value>::type>
std::vector<std::size_t> array_sizes(const TA::DistArray<Tile, Policy> &A) {
  auto const &pmap = A.get_pmap();
  TA::TiledRange const &trange = A.trange();

  const auto end = pmap->end();

  std::size_t full_size(0); // number of elements in nonzero local tiles
  for (auto it = pmap->begin(); it != end; ++it) {
    const auto tile_ord = *it;
    if (!A.is_zero(tile_ord)) {
      full_size += trange.make_tile_range(tile_ord).volume();
    }
  }

  const auto world_size = A.get_world().size();
  std::vector<std::size_t> result(world_size, 0);
  result[A.get_world().rank()] = full_size * sizeof(typename Tile::value_type);
  A.get_world().gop.sum(&result[0], world_size);

  return result;
}

/// @brief reports the total size of an array's tiles (in gigabytes)
/// @note this is a collective operation
/// @return the aggregate size of array's tiles in gigabytes
template <typename Tile, typename Policy>
double array_size(const TA::DistArray<Tile,Policy>& A){
  const auto rank_sizes = array_sizes(A);
  const auto total_size =
      std::accumulate(rank_sizes.begin(), rank_sizes.end(), std::size_t(0));
  return double(total_size) / (1<<30);
}

} // namespace utility
} // namespace mpqc

#endif // MPQC_UTILITY_ARRAYINFO_H
