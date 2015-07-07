//
// Created by Chong Peng on 6/25/15.
//

#ifndef TILECLUSTERCHEM_BLOCK_SIZE_ENGINE_H
#define TILECLUSTERCHEM_BLOCK_SIZE_ENGINE_H

#include "../include/tiledarray.h"
#include "../common/namespaces.h"

#include <algorithm>

namespace tcc{
  class TRange1Engine{

  public:

    TRange1Engine() : occ_(0ul), all_(0ul), vir_(0ul), guess_(0ul),
                      tr_occupied_(), tr_virtual_(), tr_all_()
    {}

    TRange1Engine(const std::size_t occ, const std::size_t all,
                  const std::size_t guess) : occ_(occ), all_(all), vir_(all-occ), guess_(guess)
    {
      init();
    }

    void init();

    TA::TiledRange1 get_occ_tr1() const {return tr_occupied_;}
    TA::TiledRange1 get_vir_tr1() const  {return tr_virtual_;}
    TA::TiledRange1 get_all_tr1() const {return tr_all_;}
    std::size_t get_occ() const {return occ_;}
    std::size_t get_vir() const {return vir_;}
    std::size_t get_all() const {return all_;}
    std::size_t get_occ_blocks() const {return occ_blocks_;}
    std::size_t get_vir_blocks() const {return vir_blocks_;}
    std::size_t get_all_blocks() const {return vir_blocks_+occ_blocks_;}

  private:

    TA::TiledRange1 tr_occupied();
    TA::TiledRange1 tr_virtual();
    TA::TiledRange1 tr_all();


  private:
    std::size_t occ_;
    std::size_t all_;
    std::size_t vir_;
    std::size_t guess_;

    TA::TiledRange1 tr_occupied_;
    TA::TiledRange1 tr_virtual_;
    TA::TiledRange1 tr_all_;

    std::size_t occ_blocks_;
    std::size_t vir_blocks_;

  };


  TA::TiledRange1 TRange1Engine::tr_occupied() {
    std::size_t nblocks = (guess_ < occ_) ? guess_ : occ_;
    std::size_t block_size = std::max(occ_ / nblocks, 1ul);
    std::vector<std::size_t> blocks;
    blocks.reserve(nblocks + 1);
    blocks.push_back(0);
    for (std::size_t i = block_size; i < occ_; i += block_size) {
      blocks.push_back(i);
    }
    blocks.push_back(occ_);
    return TA::TiledRange1(blocks.begin(), blocks.end());
  }

  TA::TiledRange1 TRange1Engine::tr_virtual() {
    std::size_t nblocks = (guess_ < vir_) ? guess_ : vir_;
    std::size_t block_size = std::max(vir_/ nblocks, 1ul);
    std::vector<std::size_t> blocks;
    blocks.reserve(nblocks + 1);
    blocks.push_back(0);
    for (std::size_t i = block_size; i < vir_; i += block_size) {
      blocks.push_back(i);
    }
    blocks.push_back(vir_);
    return TA::TiledRange1(blocks.begin(), blocks.end());
  }

  TA::TiledRange1 TRange1Engine::tr_all() {

    // occ part
    std::size_t nblocks = (guess_ < occ_) ? guess_ : occ_;
    std::size_t block_size = std::max(occ_ / nblocks, 1ul);
    std::vector<std::size_t> blocks;
    blocks.push_back(0);
    for (std::size_t i = block_size; i < occ_; i += block_size) {
      blocks.push_back(i);
    }
    blocks.push_back(occ_);

    //vir part
    nblocks = (guess_ < vir_) ? guess_ : vir_;
    block_size = std::max(vir_/ nblocks, 1ul);
    for (std::size_t i = occ_ + block_size; i < all_; i += block_size) {
      blocks.push_back(i);
    }
    blocks.push_back(all_);
    return TA::TiledRange1(blocks.begin(), blocks.end());
  }

  void TRange1Engine::init() {
    tr_occupied_ = tr_occupied();
    tr_virtual_ = tr_virtual();
    tr_all_ = tr_all();
    occ_blocks_ = tr_occupied_.tiles().second;
    vir_blocks_ = tr_virtual_.tiles().second;
    //std::cout << tr_occupied_ << std::endl;
    //std::cout << occ_blocks_ << " " << vir_blocks_ << std::endl;
  }
}

#endif //TILECLUSTERCHEM_BLOCK_SIZE_ENGINE_H