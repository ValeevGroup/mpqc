
#ifndef MPQC_SCF_CADFHELPERFUNCTIONS_H
#define MPQC_SCF_CADFHELPERFUNCTIONS_H

#include <tiledarray.h>


#include <unordered_map>

namespace mpqc {
namespace scf {

inline TA::DistArray<TA::TensorD, TA::SparsePolicy>
array_from_tile_map(madness::World &world, TA::TiledRange const &trange,
                    std::unordered_map<std::size_t, TA::TensorD> const &tiles) {
    auto shape_data = TA::TensorF(trange.tiles_range(), 0.0);
    for (auto &pair : tiles) {
        *(shape_data.data() + pair.first) = pair.second.norm();
    }

    TA::SparseShape<float> shape(shape_data, trange);
    TA::DistArray<TA::TensorD, TA::SparsePolicy> array(world, trange, shape);

    for (auto &pair : tiles) {
        if(!array.is_zero(pair.first) && array.is_local(pair.first)){
            array.set(pair.first, pair.second);
        }
    }

    array.truncate();
    return array;
}

inline TA::DistArray<TA::TensorD, TA::SparsePolicy> reblock_from_atoms(
      TA::DistArray<TA::TensorD, TA::SparsePolicy> const &A,
      std::unordered_map<std::size_t, std::size_t> const &output_cluster_obs,
      std::unordered_map<std::size_t, std::size_t> const &output_cluster_df,
      TA::TiledRange by_cluster_trange) {

    // auto const &pmap = *(A.pmap());
    std::unordered_map<std::size_t, TA::TensorD> tiles;
    // for (auto ord = pmap.begin(); ord != pmap.end(); ++ord) {
    for (auto ord = 0ul; ord != A.trange().tiles_range().volume(); ++ord) {
        if (!A.is_zero(ord)) {
            auto idx_in = A.trange().tiles_range().idx(ord);
            auto idx_out = idx_in;
            for (auto i = 0ul; i < A.range().rank(); ++i) {
                if (i == 0) {
                    idx_out[i] = output_cluster_df.find(idx_in[i])->second;
                } else {
                    idx_out[i] = output_cluster_obs.find(idx_in[i])->second;
                }
            }

            auto in_tile = A.find(ord).get();
            auto in_range = in_tile.range();

            auto by_cluster_ord = by_cluster_trange.tiles_range().ordinal(idx_out);
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

    return array_from_tile_map(A.world(), by_cluster_trange, tiles);
}


} // namespace scf
} // namespace mpqc
#endif //  MPQC_SCF_CADFHELPERFUNCTIONS_H
