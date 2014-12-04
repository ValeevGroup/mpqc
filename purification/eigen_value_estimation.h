#pragma once
#ifndef TCC_PURIFICATION_EIGENVALUEESTIMATION_H
#define TCC_PURIFICATION_EIGENVALUEESTIMATION_H

#include <vector>
#include <limits>
#include <cmath>
#include <iostream>

#include "../tensor/full_rank_tile.h"
#include "../tensor/low_rank_tile.h"
#include "../tensor/tile_variant.h"

// Compute the min eval guess for a row of tiles.
template <typename It>
double min_eval_guess(It first, It second) {
    auto nrows = first->get().range().extent()[0];
    std::vector<double> row_min_evals(nrows, 0);

    for (; first != second; ++first) {
        // Will need to workout something better for low rank tiles.
        auto matrix = first->get().tile().matrix();
        auto idx = first.index();

        if (idx[0] == idx[1]) {
            for (auto i = 0l; i < matrix.rows(); ++i) {
                for (auto j = 0l; j < matrix.cols(); ++j) {
                    row_min_evals[i] += (i == j) ? matrix(i, i)
                                                 : -std::abs(matrix(i, j));
                }
            }
        } else {
            for (auto i = 0l; i < matrix.rows(); ++i) {
                for (auto j = 0l; j < matrix.cols(); ++j) {
                    row_min_evals[i] -= std::abs(matrix(i, j));
                }
            }
        }
    }

    return *std::min_element(row_min_evals.begin(), row_min_evals.end());
}


template <typename Array>
double min_eval_est(Array const &H) {
    auto min = std::numeric_limits<double>::max();
    auto extent = H.trange().tiles().extent();

    // Check every row of tiles for its minimum eval guess
    for (auto tile_row_it = H.begin(); tile_row_it != H.end();) {
        auto next_tile_row = tile_row_it;
        std::advance(next_tile_row, extent[1]);

        min = std::min(min, min_eval_guess(tile_row_it, next_tile_row));

        tile_row_it = next_tile_row; // advance iterator
    }

    return min;
}


// Computing guess based on Niklasson, A. M. N., Weber, V., & Challacombe, M.
// (2005). Nonorthogonal density-matrix perturbation theory. The Journal of
// Chemical Physics, 123(4), 044107. doi:10.1063/1.1944725
template <typename Array>
Array create_eval_scaled_guess(Array const &H, Array const &S) {

    auto min_eval = min_eval_est(H);
    auto beta = min_eval - 1;

    // Compute X_0 guess
    Array HbetaS{};
    HbetaS("i,j") = H("i,j") - beta * S("i,j");

    Array X = HbetaS;

    auto matrix = X.begin()->get().tile().matrix();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(matrix);
    std::cout << "X_0 matrix = \n" << matrix << std::endl;
    std::cout << "Corret inverse = \n" << matrix.inverse() << std::endl;
    std::cout << "X_0 matrix eig_vals = \n" << es.eigenvalues().transpose()
              << std::endl;

    // Iterator till inverted
    auto iter = 0;
    while (iter < 10) { // 10 is for testing
        X("i,j") = 2 * X("i,j") - X("i,k") * HbetaS("k,l") * X("l,j");
        matrix = X.begin()->get().tile().matrix();
        es.compute(matrix);
        std::cout << "X_0 matrix = \n" << matrix << std::endl;
        std::cout << "X_0 matrix eig_vals = \n" << es.eigenvalues().transpose()
                  << std::endl;
        ++iter;
    }

    X.get_world().gop.fence();
    for (auto it = X.begin(); it != X.end(); ++it) {
        auto matrix = it->get().tile().matrix();
        std::cout << "X matrix = \n" << matrix << std::endl;
    }

    return X;
}


#endif /* end of include guard: TCC_PURIFICATION_EIGENVALUEESTIMATION_H */
