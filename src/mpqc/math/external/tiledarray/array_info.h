#ifndef MPQC4_SRC_MPQC_MATH_EXTERNAL_TILEDARRAY_ARRAY_INFO_H_
#define MPQC4_SRC_MPQC_MATH_EXTERNAL_TILEDARRAY_ARRAY_INFO_H_

#include <array>
#include <atomic>
#include <string>

#include <tiledarray.h>

#include "mpqc/util/core/exenv.h"
#include "mpqc/util/misc/string.h"

namespace mpqc {
namespace detail {

template <typename T>
unsigned long tile_clr_storage(TA::Tensor<T> const &) {
  return 0ul;
}

template <typename TileType, typename Policy>
std::array<double, 3> array_storage(TA::DistArray<TileType, Policy> const &A) {
  std::atomic_ulong full_size(0);
  std::atomic_ulong sparse_size(0);
  std::atomic_ulong low_size(0);

  TA::TiledRange const &trange = A.trange();

  auto task_is_zero = [&](unsigned long ord) {
    full_size += trange.make_tile_range(ord).volume();
  };

  auto task_is_not_zero = [&](unsigned long ord, TileType const &t) {
    const TA::Range range = trange.make_tile_range(ord);
    const auto size = range.volume();

    full_size += size;
    sparse_size += size;
    low_size += tile_clr_storage(t);
  };

  auto const &pmap = A.pmap();
  const auto end = pmap->end();

  for (auto it = pmap->begin(); it != end; ++it) {
    const auto ord = *it;
    if (A.is_local(ord)) {
      if (!A.is_zero(ord)) {
        A.world().taskq.add(task_is_not_zero, ord, A.find(ord));
      } else {
        A.world().taskq.add(task_is_zero, ord);
      }
    }
  }

  A.world().gop.fence();

  std::array<double, 3> out;
  out[0] = double(full_size);
  out[1] = double(sparse_size);
  out[2] = double(low_size);

  A.world().gop.sum(&out[0], 3);

  // reals -> GBs
  using value_type = typename TileType::value_type;
  out[0] *= sizeof(value_type) * 1e-9;
  out[1] *= sizeof(value_type) * 1e-9;
  out[2] *= sizeof(value_type) * 1e-9;

  return out;
}

/// @brief reports array size (in bytes) on each process
/// @note this is a collective operation
/// @return vector of integers, each element is the aggregate size of array's
///         tiles on the corresponding process
template <typename Tile, typename Policy,
          typename = typename std::enable_if<
              TA::detail::is_numeric<typename Tile::value_type>::value>::type>
std::vector<std::size_t> array_sizes(const TA::DistArray<Tile, Policy> &A) {
  auto const &pmap = A.pmap();
  TA::TiledRange const &trange = A.trange();

  const auto end = pmap->end();

  std::size_t full_size(0);  // number of elements in nonzero local tiles
  for (auto it = pmap->begin(); it != end; ++it) {
    const auto tile_ord = *it;
    if (!A.is_zero(tile_ord)) {
      full_size += trange.make_tile_range(tile_ord).volume();
    }
  }

  const auto world_size = A.world().size();
  std::vector<std::size_t> result(world_size, 0);
  result[A.world().rank()] = full_size * sizeof(typename Tile::value_type);
  A.world().gop.sum(&result[0], world_size);

  return result;
}

/// @brief reports the total size of an array's tiles (in gigabytes)
/// @note this is a collective operation
/// @return the aggregate size of array's tiles in gigabytes
template <typename Tile, typename Policy>
double array_size(const TA::DistArray<Tile, Policy> &A) {
  const auto rank_sizes = array_sizes(A);
  const auto total_size =
      std::accumulate(rank_sizes.begin(), rank_sizes.end(), std::size_t(0));
  return double(total_size) / (1 << 30);
}

template <typename Array>
void print_size_info(Array const &A, std::string const &name) {
  auto sizes = mpqc::detail::array_storage(A);

  ExEnv::out0(A.world()) << indent << "Printing size information for "
                         << utility::to_string(name) << std::endl
                         << incindent;

  ExEnv::out0(A.world()) << indent << "Full     = " << sizes[0] << " GB"
                         << std::endl
                         << indent << "Sparse   = " << sizes[1] << " GB"
                         << std::endl
                         << indent << "Low Rank = " << sizes[2] << " GB"
                         << std::endl
                         << decindent;
}

// average block size
std::size_t average_blocksize(TA::TiledRange1 tr1);

// compute the min and max block size in TiledRange1
std::pair<std::size_t, std::size_t> minmax_blocksize(
    TiledArray::TiledRange1 tr1);

inline void parallel_print_range_info(madness::World &world,
                                      const TA::TiledRange1 &bs_range,
                                      const std::string &name) {
  if (world.rank() == 0) {
    auto minmax_block = minmax_blocksize(bs_range);
    auto average_block = average_blocksize(bs_range);
    std::cout << name << std::endl;
    std::cout << bs_range << std::endl;
    std::cout << "Min and Max block size: " << minmax_block.first << " "
              << minmax_block.second << std::endl;
    std::cout << "Average: " << average_block << std::endl;
    std::cout << std::endl;
  }
}

}  // namespace detail
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_MATH_EXTERNAL_TILEDARRAY_ARRAY_INFO_H_
