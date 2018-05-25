//
// Created by Chong Peng on 4/5/16.
//

#include "trange1_engine.h"

namespace mpqc {
namespace utility {

TA::TiledRange1 compute_trange1(std::size_t range,
                                 std::size_t block_size) {

  std::vector<std::size_t> blocks;
  blocks.push_back(0);
  for (std::size_t i = block_size; i < range; i += block_size) {
    blocks.push_back(i);
  }
  blocks.push_back(range);
  return TA::TiledRange1(blocks.begin(), blocks.end());
}

TA::TiledRange1 join_trange1(const TA::TiledRange1& range_a, const TA::TiledRange1& range_b) {
  std::vector<std::size_t> hashmarks;
  hashmarks.reserve(range_a.tile_extent() + range_b.tile_extent() + 1);
  // append hashmarks of from the first range
  for(const auto& t: range_a)
    hashmarks.push_back(t.first);
  // shift second's tiles to start at the end of first
  const auto offset = range_a.elements_range().second - range_b.elements_range().first;
  for(const auto& t: range_b)
    hashmarks.push_back(t.first + offset);
  // add the fence hashmark
  hashmarks.push_back(range_b.elements_range().second + offset);
  // done
  return TA::TiledRange1(hashmarks.begin(), hashmarks.end());
}

TA::TiledRange1 TRange1Engine::tr_occupied() {
  std::size_t occ = get_occ();

  return compute_trange1(occ, occ_block_size_);
}

TA::TiledRange1 TRange1Engine::tr_active_occupied() {
  std::size_t active_occ = get_active_occ();

  return compute_trange1(active_occ, occ_block_size_);
}

TA::TiledRange1 TRange1Engine::tr_virtual() {
  return compute_trange1(vir_, vir_block_size_);
}

TA::TiledRange1 TRange1Engine::tr_all() {
  return compute_trange1(all_, vir_block_size_);
}

void TRange1Engine::init() {
  tr_active_occupied_ = tr_active_occupied();
  tr_occupied_ = tr_occupied();
  tr_virtual_ = tr_virtual();
  tr_all_ = tr_all();
  active_occ_blocks_ = tr_active_occupied_.tiles_range().second;
  vir_blocks_ = tr_virtual_.tiles_range().second;
  // std::cout << tr_active_occupied_ << std::endl;
  // std::cout << active_occ_blocks_ << " " << vir_blocks_ << std::endl;
}

}  // namespace utility
}  // namespace mpqc
