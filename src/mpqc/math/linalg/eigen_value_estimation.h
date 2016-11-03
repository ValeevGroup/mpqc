
#ifndef TCC_PURIFICATION_EIGENVALUEESTIMATION_H
#define TCC_PURIFICATION_EIGENVALUEESTIMATION_H

#include <vector>
#include <limits>
#include <cmath>
#include <iostream>

#include "mpqc/math/external/eigen/eigen.h"
#include <tiledarray.h>

// Compute the min eval guess for a row of tiles.
template <typename It>
double min_eval_guess(It first, It second) {
    auto nrows = first->get().range().extent()[0];
    std::vector<double> row_min_evals(nrows, 0);

    for (; first != second; ++first) {
        auto tensor = first->get();
        std::vector<std::size_t> extent = tensor.range().extent();
        auto matrix = TiledArray::eigen_map(tensor, extent[0], extent[1]);

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
std::pair<double, double> symmetric_min_max_evals(Array const &S) {
    auto const &array_extent = S.trange().tiles_range().extent();
    std::vector<Eigen::VectorXd> row_norms(array_extent[0]);

    for (auto it = S.begin(); it != S.end(); ++it) {
        auto const &tile = it->get();
        auto const extent = tile.range().extent();
        auto matrix_map = TiledArray::eigen_map(tile, extent[0], extent[1]);

        Eigen::VectorXd tile_row_sums(extent[0]);

        for (auto i = 0l; i < matrix_map.rows(); ++i) {
            auto row_sum = 0.0;
            for (auto j = 0l; j < matrix_map.cols(); ++j) {
                row_sum += std::abs(matrix_map(i, j));
            }
            tile_row_sums[i] = row_sum;
        }

        auto tile_row = it.index()[0];
        if (row_norms[tile_row].size() == 0) {
            row_norms[tile_row] = tile_row_sums;
        } else {
            row_norms[it.index()[0]] += tile_row_sums;
        }
    }

    double max = 0.0;
    double min = std::numeric_limits<double>::max();
    for (auto const &row_vector : row_norms) {
        max = std::max(max, row_vector.maxCoeff());
        min = std::min(min, row_vector.minCoeff());
    }
    return std::make_pair(min, max);
}


template <typename Array>
Array init_X(Array const &S) {
    Array X;

    auto evals = symmetric_min_max_evals(S);
    X("i,j") = 1 / (evals.first * evals.second) * S("i,j");
    return X;
}

template <typename Array>
Array invert(Array const &S) {

    Array X = init_X(S);
    Array product;
    product("i,j") = X("i,k") * S("k,j");

    auto iter = 0;
    double trace_ideal = S.trange().tiles_range().extent()[0];
    double trace_real = 0.0;
    while (iter < 1000 && std::abs(trace_real - trace_ideal) >= 1e-10) {
        X("i,j") = 2 * X("i,j") - X("i,k") * S("k,l") * X("l,j");
        product("i,j") = X("i,k") * S("k,j");
        trace_real = product("i,j").trace();
        ++iter;
    }
    X.world().gop.fence();
    return X;
}


template <typename Array>
double min_eval_est(Array const &H, Array const &S) {
    Array HSinv;
    HSinv("i,j") = invert(S)("i,k") * H("k,j");

    auto min = std::numeric_limits<double>::max();
    auto extent = HSinv.trange().tiles_range().extent();

    // Check every row of tiles for its minimum eval guess
    for (auto tile_row_it = HSinv.begin(); tile_row_it != HSinv.end();) {
        auto next_tile_row = tile_row_it;
        std::advance(next_tile_row, extent[1]);

        min = std::min(min, min_eval_guess(tile_row_it, next_tile_row));

        tile_row_it = next_tile_row; // advance iterator
    }

    return min;
}

/* // Computing guess based on Niklasson, A. M. N., Weber, V., & */
/* // Challacombe, M. */
/* // (2005). Nonorthogonal density-matrix perturbation theory. The */
/* // Journal of */
/* // Chemical Physics, 123(4), 044107. doi:10.1063/1.1944725 */
/* template <typename Array> */
/* Array create_eval_scaled_guess(Array const &H, Array const &S) { */

/*     auto min_eval = min_eval_est(H, S); */
/*     std::cout << "eval min guess = " << min_eval << std::endl; */
/*     auto eig_H = TiledArray::array_to_eigen(H); */
/*     auto eig_S = TiledArray::array_to_eigen(S); */
/*     Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd>
 * es(eig_H,eig_S); */
/*     std::cout << "Real min eval = " << es.eigenvalues()[0] << std::endl; */
/*     auto beta = min_eval - 1; */

/*     // Compute X_0 guess */
/*     Array HbetaS; */
/*     HbetaS("i,j") = H("i,j") - beta * S("i,j"); */

/*     return invert(HbetaS); */
/* } */


#endif /* end of include guard: TCC_PURIFICATION_EIGENVALUEESTIMATION_H */
