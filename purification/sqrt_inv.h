#pragma once
#ifndef TCC_PUREIFICATION_SQRTINV_H
#define TCC_PUREIFICATION_SQRTINV_H

#include "../include/tiledarray.h"
#include "../include/eigen.h"

namespace tcc {
namespace pure {

template <typename Array>
Array create_diagonal_matrix(Array const &model, double val) {
    Array diag(model.get_world(), model.trange());
    diag.set_all_local(0.0);
    auto end = diag.end();
    for (auto it = diag.begin(); it != end; ++it) {

        auto idx = it.index();
        auto diagonal_tile = std::all_of(
            idx.begin(), idx.end(), [&](typename Array::size_type const &x) {
                return x == idx.front();
            });

        if (diagonal_tile) {
            auto &tile = it->get();
            auto const &extent = tile.range().size();
            auto map = TiledArray::eigen_map(tile, extent[0], extent[1]);
            for (auto i = 0ul; i < extent[0]; ++i) {
                map(i, i) = val;
            }
        }
    }

    return diag;
}

template <typename Array>
double max_eval_est(Array const &S) {
    auto const &array_extent = S.trange().tiles().extent();
    std::vector<Eigen::VectorXd> row_norms(array_extent[0]);

    for (auto it = S.begin(); it != S.end(); ++it) {
        auto const &tile = it->get();
        auto const &extent = tile.range().size();
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
    for (auto const &row_vector : row_norms) {
        max = std::max(max, row_vector.maxCoeff());
    }
    return max;
}


template <typename Array>
void add_to_diag(Array &A, double val) {
    auto end = A.end();
    for (auto it = A.begin(); it != end; ++it) {
        auto idx = it.index();
        auto diagonal_tile = std::all_of(
            idx.begin(), idx.end(), [&](typename Array::size_type const &x) {
                return x == idx.front();
            });

        if (diagonal_tile) {
            auto &tile = it->get();
            auto const &extent = tile.range().size();
            auto map = TiledArray::eigen_map(tile, extent[0], extent[1]);
            for (auto i = 0ul; i < extent[0]; ++i) {
                map(i, i) += val;
            }
        }
    }
}

template <typename Array>
void third_order_update(Array const &S, Array &Z) {

    auto emax = max_eval_est(S);
    auto S_scale = 1 / emax;
    auto Ss_sqrt = std::sqrt(S_scale);
    Array Y = S;
    Array T;
    Array X;
    Array Sdiff;
    auto Tscale = 1.0 / 8.0;

    auto iter = 0;
    while (100 > iter++) {
        // Xn = \lambda*Yn*Zn
        X("i,j") = S_scale * Y("i,k") * Z("k,j");

        // Third order update
        T("i,j") = -10 * X("i,j") + 3 * X("i,k") * X("k,j");
        add_to_diag(T, 15);
        T("i,j") = Tscale * T("i,j");

        // Zn+1 = Zn*Tn
        Z("i,j") = Z("i,k") * T("k,j");
        // Yn+1 = Tn*Yn
        Y("i,j") = T("i,k") * Y("k,j");
        std::cout << "Iteration " << iter << std::endl;
    }

    Z("i,j") = std::sqrt(S_scale) * Z("i,j");
}


// Taken from J. Chem. Phys. 126. 124104 (2007)
// Uses the third order function, because I didn't feel like typing the longer
// ones. --Drew
template <typename Array>
Array inverse_sqrt(Array const &S) {
    Array Z = create_diagonal_matrix(S, 1.0);
    third_order_update(S, Z);
    return Z;
}


} // namespace pure
} // namespace tcc


#endif /* end of include guard: TCC_PUREIFICATION_SQRTINV_H */
