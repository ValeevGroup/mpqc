/*
 * array_info.cpp
 *
 *  Created on: Nov 1, 2016
 *      Author: evaleev
 */

#include "mpqc/math/external/tiledarray/array_info.h"

namespace mpqc {
namespace detail{

// average block size
std::size_t average_blocksize(TA::TiledRange1 tr1){
    std::vector<std::size_t> block_sizes;
    for (auto block = tr1.begin(); block != tr1.end(); ++block){
        auto block_size = block->second - block->first;
        block_sizes.push_back(block_size);
    }
    std::size_t n = block_sizes.size();
    std::size_t total = std::accumulate(block_sizes.cbegin(),block_sizes.cend(),0ul);

    std::size_t average = total / n;
    return average;
}

// compute the min and max block size in TiledRange1
std::pair<std::size_t, std::size_t>
minmax_blocksize(TiledArray::TiledRange1 tr1){

    std::vector<std::size_t> block_sizes;
    for (auto block = tr1.begin(); block != tr1.end(); ++block){
        auto block_size = block->second - block->first;
        block_sizes.push_back(block_size);
    }
    auto minmax_block_size = std::minmax_element(block_sizes.begin(),block_sizes.end());

    std::size_t min_block_size = *(minmax_block_size.first);
    std::size_t max_block_size = *(minmax_block_size.second);
    auto result = std::make_pair(min_block_size,max_block_size);
    return result;
};

}  // namespace detail
}  // namespace mpqc
