//
// Created by Chong Peng on 4/5/16.
//

#include "cc_utility.h"
#include <cmath>

namespace mpqc{
namespace cc{

std::vector<std::vector<libint2::Shell>>
reblock_basis(std::vector<libint2::Shell> shells, std::size_t blocksize){

    std::vector<std::vector<libint2::Shell>> result;

    std::vector<libint2::Shell> tmp;
    std::size_t tmp_size = 0;
    for (auto& shell : shells){

        std::size_t shell_size = shell.size();
        tmp_size += shell_size;

        if(4*tmp_size < 5*blocksize){
            tmp.push_back(shell);
        }
        else{
            result.push_back(tmp);
            tmp = std::vector<libint2::Shell>();
            tmp.push_back(shell);
            tmp_size = shell_size;
        }
    }

    // handle the boundary condition
    if(4*tmp_size < 5*blocksize){

        // if boundary is less than 2/3 of the block size
        // include it to previous block
        if(3*tmp_size < 2*blocksize){
            result.back().insert(result.back().end(),tmp.begin(),tmp.end());
        }
        else{
            result.push_back(tmp);
        }
    }
    return result;
}

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


// print progress
void print_progress(int lowprogress, int upprogress, int total){
    int divide = 10;
    if(total < 10){
        divide = total;
        int percent = 100*double(upprogress)/divide;
        std::cout << percent << "% done." << std::endl;
    }else{
        int increase = std::round(double(total)/divide);
        std::vector<int> progress_points;
        for(int i = 0; i < total; i += increase){
            progress_points.push_back(i);
        }

        for (int i = lowprogress; i < upprogress; i++){
            if (std::find(progress_points.begin(),progress_points.end(),i) != progress_points.end()){
                int percent = 100*double(i)/total;
                std::cout << percent << "% done." << std::endl;
            }
        }
    }

}


} // end of namespace cc
} // end of namespace mpqc
