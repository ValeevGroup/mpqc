//
// Created by Chong Peng on 4/5/16.
//

#include "trange1_engine.h"

namespace mpqc {

TA::TiledRange1 TRange1Engine::compute_range(std::size_t range,
                                             std::size_t block_size) {
  std::vector<std::size_t> blocks;
  blocks.push_back(0);
  for (std::size_t i = block_size; i < range; i += block_size) {
    blocks.push_back(i);
  }

  // if the boundary is less than 2/3 of the block size
  // include it to previous block
  if (3 * (range - blocks.back()) <= 2 * block_size && blocks.size() > 1) {
    blocks.back() = range;
  }
  // if not, add a new block
  else {
    blocks.push_back(range);
  }

  return TA::TiledRange1(blocks.begin(), blocks.end());
}

TA::TiledRange1 TRange1Engine::tr_occupied() {
  std::size_t occ = get_occ();

  return compute_range(occ, occ_block_size_);
}

TA::TiledRange1 TRange1Engine::tr_active_occupied() {
  std::size_t active_occ = get_active_occ();

  return compute_range(active_occ, occ_block_size_);
}

TA::TiledRange1 TRange1Engine::tr_virtual() {
  return compute_range(vir_, vir_block_size_);
}

TA::TiledRange1 TRange1Engine::tr_all() {
  return compute_range(all_, vir_block_size_);
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

}  // namespace mpqc