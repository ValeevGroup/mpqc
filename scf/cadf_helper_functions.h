#pragma once
#ifndef MPQC_SCF_CADFHELPERFUNCTIONS_H
#define MPQC_SCF_CADFHELPERFUNCTIONS_H

#include "../include/tiledarray.h"
#include "../common/typedefs.h"

#include <unordered_map>

namespace mpqc {
namespace scf {

inline TA::DistArray<TA::TensorD, SpPolicy>
array_from_tile_map(madness::World &world, TA::TiledRange const &trange,
                    std::unordered_map<std::size_t, TA::TensorD> const &tiles) {
    auto shape_data = TA::TensorF(trange.tiles(), 0.0);
    for (auto &pair : tiles) {
        *(shape_data.data() + pair.first) = pair.second.norm();
    }

    TA::SparseShape<float> shape(shape_data, trange);
    TA::DistArray<TA::TensorD, SpPolicy> array(world, trange, shape);

    for (auto &pair : tiles) {
        array.set(pair.first, pair.second);
    }

    array.truncate();
    return array;
}

inline TA::DistArray<TA::TensorD, SpPolicy> reblock_from_atoms(
      TA::DistArray<TA::TensorD, SpPolicy> const &A,
      std::unordered_map<std::size_t, std::size_t> const &output_cluster_obs,
      std::unordered_map<std::size_t, std::size_t> const &output_cluster_df,
      TA::TiledRange by_cluster_trange) {

    auto const &pmap = *(A.get_pmap());
    std::unordered_map<std::size_t, TA::TensorD> tiles;
    for (auto ord = pmap.begin(); ord != pmap.end(); ++ord) {
        if (!A.is_zero(*ord)) {
            auto idx_in = A.trange().tiles().idx(*ord);
            auto idx_out = idx_in;
            for (auto i = 0; i < A.range().rank(); ++i) {
                if (i == 0) {
                    idx_out[i] = output_cluster_df.find(idx_in[i])->second;
                } else {
                    idx_out[i] = output_cluster_obs.find(idx_in[i])->second;
                }
            }


            auto in_tile = A.find(*ord).get();
            auto in_range = in_tile.range();

            auto by_cluster_ord = by_cluster_trange.tiles().ordinal(idx_out);
            auto out_range = by_cluster_trange.make_tile_range(by_cluster_ord);

            auto write_inner_tile =
                  [&](TA::TensorD const &inner, TA::TensorD &out) {
                auto in_start = inner.range().lobound();
                auto in_end = inner.range().upbound();

                out.block(in_start, in_end) = inner;
            };

            auto out_pair = tiles.find(by_cluster_ord);
            if (out_pair == tiles.end()) {
                auto out_tensor = TA::TensorD(out_range, 0.0);
                write_inner_tile(in_tile, out_tensor);
                tiles[by_cluster_ord] = std::move(out_tensor);
            } else {
                write_inner_tile(in_tile, out_pair->second);
            }
        }
    }

    return array_from_tile_map(A.get_world(), by_cluster_trange, tiles);
}


} // namespace scf
} // namespace mpqc
#endif //  MPQC_SCF_CADFHELPERFUNCTIONS_H
