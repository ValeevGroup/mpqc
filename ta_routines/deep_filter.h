#pragma once
#ifndef TCC_ARRAYOPS_DEEPFILTER_H
#define TCC_ARRAYOPS_DEEPFILTER_H

#include "../include/tiledarray.h"
#include "../common/namespaces.h"
#include <vector>
#include <array>
#include <utility>

namespace tcc {
namespace array_ops {

template <typename T, unsigned int DIM, typename Tile>
TA::Array<T, DIM, Tile, TA::SparsePolicy>
deep_filter(TA::Array<T, DIM, Tile, TA::SparsePolicy> const &t,
            std::array<std::pair<unsigned long, unsigned long>, DIM> ranges) {

    std::vector<TA::TiledRange1> trange1s;
    std::vector<unsigned int> tile_starts;
    for (auto i = 0u; i < DIM; ++i) {
        auto const &trange1 = t.trange().data()[i];

        auto first_tile_it = std::find_if(
              trange1.begin(), trange1.end(),
              [=](std::pair<unsigned long, unsigned long> const &a) {
                  return ranges[i].first == a.first;
              });

        auto last_tile_it = std::find_if(
              trange1.begin(), trange1.end(),
              [=](std::pair<unsigned long, unsigned long> const &a) {
                  return ranges[i].second == a.first;
              });
        
        auto first_tile = trange1.element2tile(first_tile_it->first);
        tile_starts.push_back(first_tile);

        std::vector<unsigned long> blocking;
        for (; first_tile_it != last_tile_it; ++first_tile_it) {
            blocking.push_back(first_tile_it->first - ranges[i].first);
        }
        if(last_tile_it != trange1.end()){
            blocking.push_back(first_tile_it->first - ranges[i].first);
        } else {
            --first_tile_it;
            blocking.push_back(first_tile_it->second - ranges[i].first);
        }

        trange1s.emplace_back(blocking.begin(), blocking.end());
    }

    TA::TiledRange new_trange(trange1s.begin(), trange1s.end());

    TA::Tensor<float> tile_norms(new_trange.tiles(),
                                 new_trange.elements().volume()
                                 / new_trange.tiles().volume());

    TA::SparseShape<float> shape(t.get_world(), tile_norms, new_trange);

    TA::Array<T, DIM, Tile, TA::SparsePolicy> new_array(t.get_world(),
                                                        new_trange, shape);


    auto pmap = new_array.get_pmap();
    const auto end = pmap->end();
    for (auto it = pmap->begin(); it != end; ++it) {
        auto range = new_trange.make_tile_range(*it);
        auto idx = new_trange.tiles().idx(*it);
        for (auto i = 0u; i < DIM; ++i) {
            idx[i] += tile_starts[i];
        }

        auto old_ord = t.trange().tiles().ordinal(idx);
        if (!t.is_zero(old_ord)) {
            auto old_tile = t.find(old_ord).get();
            new_array.set(*it, Tile(range, old_tile.data()));
        } else {
            new_array.set(*it, Tile(range, 0.0));
        }
    }

    new_array.truncate();

    return new_array;
}


} // namespace array_ops
} // namespace tcc

#endif // TCC_ARRAYOPS_DEEPFILTER_H
