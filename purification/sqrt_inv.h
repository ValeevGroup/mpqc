#pragma once
#ifndef TCC_PUREIFICATION_SQRTINV_H
#define TCC_PUREIFICATION_SQRTINV_H

#include "../include/tiledarray.h"
#include "../include/eigen.h"

#include "diagonal_array.h"

#include <limits>
#include <type_traits>
#include <cmath>

namespace tcc {
namespace pure {

template <typename T, typename AT>
std::array<TiledArray::Tensor<T, AT>, 2>
    eigen_estimator(std::array<TiledArray::Tensor<T, AT>, 2> &result,
                    TiledArray::Tensor<T, AT> const &tile) {

    typedef typename TiledArray::Tensor<T, AT>::size_type size_type;

    if (result[0].empty()) {
        // Construct result tensors
        const std::array<size_type, 1> range_start
            = {{tile.range().start()[0]}};
        const std::array<size_type, 1> range_finish
            = {{tile.range().finish()[0]}};
        TiledArray::Range range(range_start, range_finish);
        result[0] = TiledArray::Tensor<T, AT>{range, 0.0};
        result[1] = TiledArray::Tensor<T, AT>{range, 0.0};
    }

    // Sum the rows of tile into result
    auto reduce_op =
        [](T &restrict result, const T arg) { result += std::abs(arg); };

    TiledArray::math::row_reduce(tile.range().size()[0], tile.range().size()[1],
                                 tile.data(), result[0].data(), reduce_op);


    auto const &start = tile.range().start();
    auto const &finish = tile.range().finish();
    auto const &weight = tile.range().weight();

    size_type const start_max = *std::max_element(start.begin(), start.end());
    size_type const finish_min
        = *std::min_element(finish.begin(), finish.end());

    const size_type n = tile.range().dim();
    if (start_max < finish_min) {
        // Compute the first and last ordinal index
        size_type tile_first = 0ul, tile_last = 0ul, tile_stride = 0ul;
        for (size_type i = 0ul; i < n; ++i) {
            const size_type start_i = start[i];
            const size_type weight_i = weight[i];

            tile_first += (start_max - start_i) * weight_i;
            tile_last += (finish_min - start_i) * weight_i;
            tile_stride += weight_i;
        }
        size_type result_first = tile_first / weight[0];

        // Compute the trace
        const T *restrict const tile_data = tile.data();
        T *restrict const result_data = result[1].data();
        for (; tile_first < tile_last;
             tile_first += tile_stride, ++result_first) {
            result_data[result_first] = tile_data[tile_first];
        }
    }

    return result;
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
class pair_accumulator {
  public:
    using result_type = std::array<typename Array::value_type, 2>;
    using argument_type = typename Array::value_type;
    using value_type = argument_type;

    pair_accumulator(TiledArray::Range const &range) : range_{range} {}
    pair_accumulator(pair_accumulator const &other) : range_(other.range_) {}
    pair_accumulator &operator=(pair_accumulator const &other) {
        range_ = other.range_;
        return *this;
    }

    result_type operator()() const {
        return result_type{{value_type{}, value_type{}}};
    }

    result_type operator()(result_type result) {
        if (result[0].empty()) {
            result[0] = value_type{range_, 0.0};
            result[1] = value_type{range_, 0.0};
        }

        return result;
    }

    void operator()(result_type &result, result_type const &arg) const {
        if (result[0].empty()) {
            result[0] = value_type{range_, 0.0};
            result[1] = value_type{range_, 0.0};
        }
        assert(!result.empty() && !arg.empty());
        result[0] += arg[0];
        result[1] += arg[1];
    }

    void operator()(result_type &result, argument_type const &arg) const {
        eigen_estimator(result, arg);
    }

  private:
    TiledArray::Range range_;
};

class pair_smasher {
  public:
    using tensor_type = TiledArray::Tensor<double>;
    using argument_type = std::array<tensor_type, 2>;
    using result_type = std::array<double, 2>;

    result_type operator()() const {
        return result_type{{std::numeric_limits<double>::max(),
                            std::numeric_limits<double>::min()}};
    }

    const result_type &operator()(const result_type &result) const {
        return result;
    }

    void operator()(result_type &result, result_type const &arg) const {
        result[0] = std::min(arg[0], result[0]);
        result[1] = std::max(arg[1], result[1]);
    }

    void operator()(result_type &result, argument_type const &arg) {
        const auto n = arg[0].size();
        const auto &row_sums = arg[0];
        const auto &diag_elems = arg[1];

        for (auto i = 0ul; i < n; ++i) {
            auto dval = diag_elems[i];
            auto corrected_sum = row_sums[i] - std::abs(dval);
            result[0] = std::min(result[0], dval - corrected_sum);
            result[1] = std::max(result[1], dval + corrected_sum);
        }
    }
};


template <typename Array>
std::array<typename Array::value_type::numeric_type, 2>
eval_guess(Array const &A) {

    A.get_world().gop.fence();

    using global_accumlator = pair_accumulator<Array>;
    std::vector<TiledArray::detail::ReduceTask<global_accumlator>> tasks;
    tasks.reserve(A.trange().tiles().size()[0]);
    for (auto i = 0ul; i < A.trange().tiles().size()[0]; ++i) {
        auto row_ranges = A.trange().data()[0].tile(i);
        std::array<typename Array::size_type, 1> start = {{row_ranges.first}};
        std::array<typename Array::size_type, 1> finish = {{row_ranges.second}};
        TiledArray::Range range(start, finish);
        tasks.emplace_back(A.get_world(), global_accumlator(range));
    }

    auto end = A.end();
    for (auto it = A.begin(); it != end; ++it) {
        tasks[it.index()[0]].add(it->get());
    }

    TiledArray::detail::ReduceTask<pair_smasher> local_reduce(A.get_world(),
                                                              pair_smasher{});

    auto counter = 0;
    for (auto i = 0ul; i < A.trange().tiles().size()[0]; ++i) {
        auto pair = tasks[i].submit();
        auto row_ranges = A.trange().data()[0].tile(i);
        std::array<typename Array::size_type, 1> start = {{row_ranges.first}};
        std::array<typename Array::size_type, 1> finish = {{row_ranges.second}};
        TiledArray::Range range(start, finish);
        madness::DistributedID key(A.id(), counter++);
        local_reduce.add(
            A.get_world().gop.all_reduce(key, pair, global_accumlator{range}));
    }

    return local_reduce.submit();
}


/* template <typename Array> */
/* std::pair<double, double> correct_scale(Array const &S) { */
/*     const auto dim = S.elements().size()[0]; */
/*     Eigen::MatrixXd eig_S = Eigen::MatrixXd::Zero(dim, dim); */
/*     for (auto it = S.begin(); it != S.end(); ++it) { */
/*         auto const &tile = (*it).get(); */
/*         const auto i = tile.range().start()[0]; */
/*         const auto j = tile.range().start()[1]; */
/*         const auto m = tile.range().size()[0]; */
/*         const auto n = tile.range().size()[1]; */
/*         eig_S.block(i, j, m, n) = TiledArray::eigen_map(tile, m, n); */
/*     } */
/*     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(eig_S); */
/*     auto max = eig.eigenvalues().maxCoeff(); */
/*     auto min = eig.eigenvalues().minCoeff(); */
/*     return std::make_pair(min, max); */
/* } */

template <typename Array>
void third_order_update(Array const &S, Array &Z) {

    // Calculate the lambda parameter, since we are currently only
    // trying to
    // invert Posative defininate matrices, we will cap the smallest
    // eigen value
    // guess at 0.
    auto spectral_range = eval_guess(S);
    const auto max_eval = spectral_range[1];
    const auto min_eval = std::max(0.0, spectral_range[0]);
    auto S_scale = 2.0 / (max_eval + min_eval);

    Array Y = S;
    Array T;
    Array X;
    auto Tscale = 1.0 / 8.0;

    auto ident = create_diagonal_matrix(Z, 1.0);
    Array approx_zero;
    auto iter = 0;
    auto norm_diff = std::numeric_limits<double>::max();
    while (norm_diff > 1.0e-13 && iter < 25) {
        // Xn = \lambda*Yn*Z
        X("i,j") = S_scale * Y("i,k") * Z("k,j");
        X.truncate();


        // Third order update
        T("i,j") = -10 * X("i,j") + 3 * X("i,k") * X("k,j");
        add_to_diag(T, 15);
        T("i,j") = Tscale * T("i,j");


        // Updating Z and Y
        Z("i,j") = Z("i,k") * T("k,j"); // Zn+1 = Zn*Tn
        Y("i,j") = T("i,k") * Y("k,j"); // Yn+1 = Tn*Yn

        Z.truncate();
        Y.truncate();

        approx_zero("i,j") = X("i,j") - T("i,j");

        const auto current_norm = approx_zero("i,j").norm().get();

        if (current_norm >= norm_diff) { // Once norm is increasing exit!
            Z("i,j") = std::sqrt(S_scale) * Z("i,j");
            return;
        }
        norm_diff = current_norm;
        ++iter;
    }

    Z("i,j") = std::sqrt(S_scale) * Z("i,j");
}


// Taken from J. Chem. Phys. 126. 124104 (2007)
// Uses the third order function, because I didn't feel like typing the
// longer
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