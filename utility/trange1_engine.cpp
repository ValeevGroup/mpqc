//
// Created by Chong Peng on 4/5/16.
//

#include "trange1_engine.h"

namespace mpqc{

TA::TiledRange1 TRange1Engine::tr_occupied() {
    std::size_t block_size = occ_block_size_;
    std::size_t actual_occ = get_actual_occ();
    std::vector<std::size_t> blocks;
    blocks.push_back(0);
    for (std::size_t i = block_size; i < actual_occ; i += block_size) {
        blocks.push_back(i);
    }

    // if the boundary is less than 2/3 of the block size
    // include it to previous block
    if(3*(actual_occ - blocks.back()) <= 2*block_size && blocks.size() > 1){
        blocks.back() = actual_occ;
    }
        // if not, add a new block
    else{
        blocks.push_back(actual_occ);
    }

    return TA::TiledRange1(blocks.begin(), blocks.end());
}

TA::TiledRange1 TRange1Engine::tr_virtual() {
    std::size_t block_size = vir_block_size_;
    std::vector<std::size_t> blocks;
    blocks.push_back(0);
    for (std::size_t i = block_size; i < vir_; i += block_size) {
        blocks.push_back(i);
    }

    // if the boundary is less than 2/3 of the block size
    // include it to previous block
    if(3*(vir_ - blocks.back()) <= 2*block_size && blocks.size() > 1){
        blocks.back() = vir_;
    }
        // if not, add a new block
    else{
        blocks.push_back(vir_);
    }

    return TA::TiledRange1(blocks.begin(), blocks.end());
}

TA::TiledRange1 TRange1Engine::tr_all() {

    // occ part
    std::size_t block_size = occ_block_size_;
    std::size_t actual_occ = get_actual_occ();
    std::vector<std::size_t> blocks;
    blocks.push_back(0);
    for (std::size_t i = block_size; i < actual_occ; i += block_size) {
        blocks.push_back(i);
    }

    if(3*(actual_occ - blocks.back()) <= 2*block_size && blocks.size() > 1){
        blocks.back() = actual_occ;
    }
        // if not, add a new block
    else{
        blocks.push_back(actual_occ);
    }

    //vir part
    block_size = vir_block_size_;
    std::size_t actual_all = get_actual_all();
    for (std::size_t i = actual_occ + block_size; i < actual_all; i += block_size) {
        blocks.push_back(i);
    }

    // if the boundary is less than 2/3 of the block size
    // include it to previous block
    if(3*(actual_all - blocks.back()) <= 2*block_size){
        blocks.back() = actual_all;
    }
        // if not, add a new block
    else{
        blocks.push_back(actual_all);
    }

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


}// end of namespace mpqc