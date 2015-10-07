#pragma once
#ifndef TCC_INTEGRALS_TILEENGINE_H
#define TCC_INTEGRALS_TILEENGINE_H

#include "../include/libint.h"
#include "../include/btas.h"

#include "../molecule/cluster_collapse.h"

#include "../basis/cluster_shells.h"
#include "../basis/basis.h"

#include "integral_engine_pool.h"
#include "tile_engine_helper.h"

#include "../common/typedefs.h"

namespace tcc {
namespace integrals {

template <typename T>
class TileEngine {
  public:
    template <std::size_t N>
    using RangeNd = detail::RangeNd<N>;

    template <std::size_t N>
    using TileType = detail::TileType<N>;


    template <std::size_t N>
    using ShellPtrs = std::array<std::shared_ptr<std::vector<ShellVec>>, N>;


    template <typename Index, typename SharedEnginePool, std::size_t N>
    TileType<N> operator()(Index cluster_idx, SharedEnginePool &engines,
                           ShellPtrs<N> const &shell_ptrs) const {
        return detail::compute_tile(engines->local(),
                                    get_shells(cluster_idx, shell_ptrs));
    }


  private:
    template <typename Index, std::size_t N>
    std::array<ShellVec, N>
    get_shells(Index cluster_idx, ShellPtrs<N> const &shell_ptrs) const {
        std::array<ShellVec, N> clusters;

        for (auto i = 0ul; i < N; ++i) {
            const auto cluster_number = cluster_idx[i];
            const auto &shell_vec = *shell_ptrs[i];
            clusters[i] = shell_vec[cluster_number];
        }

        return clusters;
    }

}; // TileEngine

} // namespace integrals
} // namespace tcc

#endif /* end of include guard: TCC_INTEGRALS_TILEENGINE_H */
