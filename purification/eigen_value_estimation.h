#pragma once
#ifndef TCC_PURIFICATION_EIGENVALUEESTIMATION_H
#define TCC_PURIFICATION_EIGENVALUEESTIMATION_H

#include <vector>
#include <limits>
#include <cmath>
#include <iostream>

template <typename Array>
Array create_eval_scaled_guess(Array const &a) {
    Array guess(a.get_world(), a.trange());

    auto extent = a.trange().tiles().extent();
    auto ncols = extent[1];

    auto max = 0.0;
    auto min = std::numeric_limits<double>::max();

    // Loop over rows of tiles. 
    for (auto row_it = a.begin(); row_it != a.end();) {

        auto tile_nrows = row_it->get().range().extent()[0];
        std::vector<double> block_row_sums(tile_nrows,0);
        std::vector<double> diags(tile_nrows,0);

        auto row_it_copy = row_it;
        std::advance(row_it_copy, ncols);

        // For loop over each tile in the row of tiles. 
        for (auto tile_it = row_it; tile_it != row_it_copy; ++tile_it) {

            auto matrix = tile_it->get().tile().matrix();
            auto tile_ncols = matrix.cols(); // Width of current tile. 
            auto idx = tile_it.index();

            // accumulate sum across rows. 
            for (auto i = 0ul; i < tile_nrows; ++i) {
                for (auto j = 0ul; j < tile_ncols; ++j) {
                    if(idx[0] == idx[1] && i == j){
                        diags[i] = matrix(i,j);
                    }
                    block_row_sums[i] += std::abs(matrix(i, j));
                }
            }

        }

        for(auto i = 0ul; i < diags.size(); ++i){
            min = std::min(diags[i] - block_row_sums[i], min);
            max = std::max(diags[i] + block_row_sums[i], max);
        }

        // Save vector for this row of tiles. 
        row_it = row_it_copy; // advance the iterator
    }

    std::cout << "Max row sum was " << max << std::endl;
    std::cout << "Min row sum was " << min << std::endl;
    double scale = 1.0/(static_cast<double>(max - min));
    std::cout << "The guess scale factor was " << scale << std::endl;

    guess("i,j") = a("i,j");

    for(auto it = guess.begin(); it != guess.end(); ++it){
        auto index = it.index();
        if(index[0] == index[1]){
            // Scale the diagonal by the max value; 
        }
    }

    return guess;
}


#endif /* end of include guard: TCC_PURIFICATION_EIGENVALUEESTIMATION_H */
