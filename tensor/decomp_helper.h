#pragma once

#include "tensor_fwd.h"
#include <tuple>
#include <iostream>

namespace tcc {
namespace tensor {
namespace detail {


TATensor to_full(TARange const &range, std::vector<TATensor> const &ts) {
    if (ts.size() == 1) {
        return ts[0].clone();
    } else if (ts.size() == 2) { // There must be at least 2
        math::GemmHelper h(madness::cblas::CBLAS_TRANSPOSE::NoTrans,
                           madness::cblas::CBLAS_TRANSPOSE::NoTrans,
                           range.dim(), ts[0].range().dim(),
                           ts[1].range().dim());

        return ts[0].gemm(ts[1], 1.0, std::move(h));
    }
}

template <int LeftDivisions, int RightDivisions>
std::tuple<TARange, std::vector<TATensor>, double>
add(TARange const &Lrange, std::vector<TATensor> const &left_t, double lcut,
    TARange const &Rrange, std::vector<TATensor> const &right_t, double rcut,
    double factor);

template <>
std::tuple<TARange, std::vector<TATensor>, double>
add<1, 1>(TARange const &Lrange, std::vector<TATensor> const &left_t,
          double lcut, TARange const &Rrange,
          std::vector<TATensor> const &right_t, double rcut, double factor) {
    std::vector<TATensor> vec = {left_t[0].add(right_t[0], factor)};
    return std::make_tuple(Lrange, std::move(vec), std::max(lcut, rcut));
}

template <>
std::tuple<TARange, std::vector<TATensor>, double>
add<1, 2>(TARange const &Lrange, std::vector<TATensor> const &left_t,
          double lcut, TARange const &Rrange,
          std::vector<TATensor> const &right_t, double rcut, double factor) {

    auto full_right = to_full(Rrange, right_t);
    std::vector<TATensor> vec = {left_t[0].add(full_right, factor)};
    return std::make_tuple(Lrange, std::move(vec), std::max(lcut, rcut));
}

template <>
std::tuple<TARange, std::vector<TATensor>, double>
add<2, 1>(TARange const &Lrange, std::vector<TATensor> const &left_t,
          double lcut, TARange const &Rrange,
          std::vector<TATensor> const &right_t, double rcut, double factor) {

    auto full_left = to_full(Lrange, left_t);
    std::vector<TATensor> vec = {full_left.add(right_t[0], factor)};
    return std::make_tuple(Lrange, std::move(vec), std::max(lcut, rcut));
}

template <>
std::tuple<TARange, std::vector<TATensor>, double>
add<2, 2>(TARange const &Lrange, std::vector<TATensor> const &left_t,
          double lcut, TARange const &Rrange,
          std::vector<TATensor> const &right_t, double rcut, double factor) {

    auto const &l0_range = left_t[0].range();
    auto const &r0_range = right_t[0].range();

    auto const &l1_range = left_t[1].range();
    auto const &r1_range = right_t[1].range();


    // Before doing to much work see if rounded addition makes sense
    {
        const auto Lout_size = l0_range.volume() + r0_range.volume();
        const auto Rout_size = l1_range.volume() + r1_range.volume();
        const auto Full_size = Lrange.volume();

        if (Full_size <= Lout_size + Rout_size) {
            std::vector<TATensor> vec = {
                to_full(Lrange, left_t).add(to_full(Rrange, right_t), factor)};
            return std::make_tuple(Lrange, std::move(vec),
                                   std::max(lcut, rcut));
        }
    }


    assert(l0_range.dim() == r0_range.dim());
    assert(l1_range.dim() == r1_range.dim());

    TARange Lrange_out;
    TARange Rrange_out;

    // Compute range for the left tensor
    if (l0_range.dim() == 2) {
        Lrange_out = TARange(l0_range.size()[0],
                             l0_range.size()[1] + r0_range.size()[1]);
    }
    // If adding order 3 tensors grow the rank dim
    if (l0_range.dim() == 3) {
        Lrange_out = TARange(l0_range.size()[0], l0_range.size()[1],
                             l0_range.size()[2] + r0_range.size()[2]);
    }

    // Compute range for the right tensor
    if (l1_range.dim() == 2) {
        Rrange_out = TARange(l1_range.size()[0] + r1_range.size()[0],
                             l1_range.size()[1]);
    }
    // If adding order 3 tensors grow the rank dim
    if (l1_range.dim() == 3) {
        Rrange_out = TARange(l1_range.size()[0] + r1_range.size()[0],
                             l1_range.size()[1], l1_range.size()[2]);
    }

    TATensor Left(Lrange_out);
    TATensor Right(Rrange_out);

    // Rounded addition into the Left tensor
    if (l0_range.dim() == 2) {
        auto Lmap
            = TA::eigen_map(Left, Lrange_out.size()[0], Lrange_out.size()[1]);

        auto l0_map
            = TA::eigen_map(left_t[0], l0_range.size()[0], l0_range.size()[1]);

        auto r0_map
            = TA::eigen_map(right_t[0], r0_range.size()[0], r0_range.size()[1]);

        Lmap.leftCols(l0_map.cols()) = factor * l0_map;
        Lmap.rightCols(r0_map.cols()) = r0_map;
    }
    if (l0_range.dim() == 3) {
        auto Lmap = TA::eigen_map(Left, Lrange_out.size()[0],
                                  Lrange_out.size()[1] * Lrange_out.size()[2]);

        auto l0_map = TA::eigen_map(left_t[0], l0_range.size()[0],
                                    l0_range.size()[1] * l0_range.size()[2]);

        auto r0_map = TA::eigen_map(right_t[0], r0_range.size()[0],
                                    r0_range.size()[1] * r0_range.size()[2]);

        Lmap.leftCols(l0_map.cols()) = factor * l0_map;
        Lmap.rightCols(r0_map.cols()) = r0_map;
    }

    // Rounded addition into the Right tensor
    if (l0_range.dim() == 2) {
        auto Rmap
            = TA::eigen_map(Right, Rrange_out.size()[0], Rrange_out.size()[1]);

        auto l1_map
            = TA::eigen_map(left_t[1], l1_range.size()[0], l1_range.size()[1]);

        auto r1_map
            = TA::eigen_map(right_t[1], r1_range.size()[0], r1_range.size()[1]);

        Rmap.topRows(l1_map.rows()) = l1_map;
        Rmap.bottomRows(r1_map.rows()) = r1_map;
    }
    if (l0_range.dim() == 3) {
        auto Rmap
            = TA::eigen_map(Right, Rrange_out.size()[0] * Rrange_out.size()[1],
                            Rrange_out.size()[2]);

        auto l1_map
            = TA::eigen_map(left_t[1], l1_range.size()[0] * l1_range.size()[1],
                            l1_range.size()[2]);

        auto r1_map
            = TA::eigen_map(right_t[1], r1_range.size()[0] * r1_range.size()[1],
                            r1_range.size()[2]);

        Rmap.topRows(l1_map.rows()) = l1_map;
        Rmap.bottomRows(r1_map.rows()) = r1_map;
    }

    std::vector<TATensor> vec = {std::move(Left), std::move(Right)};
    return std::make_tuple(Lrange, std::move(vec), std::max(lcut, rcut));
}

template <int LeftDivisions, int RightDivisions>
std::tuple<TARange, std::vector<TATensor>, double>
gemm(TARange const &out_range, TARange const &Lrange,
     std::vector<TATensor> const &left_t, double lcut, TARange const &Rrange,
     std::vector<TATensor> const &right_t, double rcut, double factor,
     math::GemmHelper const &gh);

template <>
std::tuple<TARange, std::vector<TATensor>, double>
gemm<1, 1>(TARange const &out_range, TARange const &Lrange,
           std::vector<TATensor> const &left_t, double lcut,
           TARange const &Rrange, std::vector<TATensor> const &right_t,
           double rcut, double factor, math::GemmHelper const &gh) {

    std::vector<TATensor> vec = {left_t[0].gemm(right_t[0], factor, gh)};
    return std::make_tuple(out_range, std::move(vec), std::max(lcut, rcut));
}

template <>
std::tuple<TARange, std::vector<TATensor>, double>
gemm<1, 2>(TARange const &out_range, TARange const &Lrange,
           std::vector<TATensor> const &left_t, double lcut,
           TARange const &Rrange, std::vector<TATensor> const &right_t,
           double rcut, double factor, math::GemmHelper const &gh) {

    // Calculate how many dims to contract over. For now assume this
    // is always possible
    const auto dims_contracted = gh.num_contract_ranks();
    const auto ranks_from_right = right_t[0].range().dim() - dims_contracted;
    const auto Lout_dim = Lrange.dim() - dims_contracted + ranks_from_right;

    if (gh.right_op() == madness::cblas::CBLAS_TRANSPOSE::NoTrans) {
        auto left_gh = math::GemmHelper(gh.left_op(), gh.right_op(), Lout_dim,
                                        left_t[0].range().dim(),
                                        right_t[0].range().dim());
        std::cout << "Lout_dim = " << Lout_dim << std::endl;
        std::cout << "Left = " << left_t[0].range().dim() << std::endl;
        std::cout << "right = " << right_t[0].range().dim() << std::endl;
        std::cout << std::endl;

        std::vector<TATensor> vec
            = {left_t[0].gemm(right_t[0], factor, left_gh), right_t[1]};
        return std::make_tuple(out_range, std::move(vec), std::min(lcut, rcut));

    } else if (gh.right_op() == madness::cblas::CBLAS_TRANSPOSE::Trans) {
        assert(false); // I don't have this case yet or know exactly how to
                       // handle it.
    }
}

template <>
std::tuple<TARange, std::vector<TATensor>, double>
gemm<2, 1>(TARange const &out_range, TARange const &Lrange,
           std::vector<TATensor> const &left_t, double lcut,
           TARange const &Rrange, std::vector<TATensor> const &right_t,
           double rcut, double factor, math::GemmHelper const &gh) {

    // Calculate how many dims to contract over. For now assume this
    // is always possible
    const auto dims_contracted = gh.num_contract_ranks();
    const auto ranks_from_left = left_t[1].range().dim() - dims_contracted;
    const auto Rout_dim = Rrange.dim() - dims_contracted + ranks_from_left;

    if (gh.left_op() == madness::cblas::CBLAS_TRANSPOSE::NoTrans) {
        auto left_gh = math::GemmHelper(gh.left_op(), gh.right_op(), Rout_dim,
                                        left_t[1].range().dim(),
                                        right_t[0].range().dim());

        std::vector<TATensor> vec = {
            left_t[0], left_t[1].gemm(right_t[0], factor, left_gh), right_t[1]};
        return std::make_tuple(out_range, std::move(vec), std::min(lcut, rcut));

    } else if (gh.left_op() == madness::cblas::CBLAS_TRANSPOSE::Trans) {
        assert(false); // I don't have this case yet or know exactly how to
                       // handle it.
    }
}

template <>
std::tuple<TARange, std::vector<TATensor>, double>
gemm<2, 2>(TARange const &out_range, TARange const &Lrange,
           std::vector<TATensor> const &left_t, double lcut,
           TARange const &Rrange, std::vector<TATensor> const &right_t,
           double rcut, double factor, math::GemmHelper const &gh) {

    // Calculate how many dims to contract over. For now assume this
    // is always possible
    const auto dims_contracted = gh.num_contract_ranks();
    const auto ranks_from_right = right_t[0].range().dim() - dims_contracted;
    const auto ranks_from_left = left_t[1].range().dim() - dims_contracted;
    const auto dims_mid = ranks_from_right + ranks_from_left;

    if (gh.right_op() == madness::cblas::CBLAS_TRANSPOSE::NoTrans
        && gh.left_op() == gh.right_op()) {

        auto left_gh = math::GemmHelper(gh.left_op(), gh.right_op(), dims_mid,
                                        left_t[1].range().dim(),
                                        right_t[0].range().dim());

        auto mid = left_t[1].gemm(right_t[0], 1.0, left_gh);

        if (left_t[0].range().volume() <= right_t[1].range().volume()) {
            auto mid_gh = math::GemmHelper(
                gh.left_op(), gh.right_op(), right_t[1].range().dim(),
                mid.range().dim(), right_t[1].range().dim());

            std::vector<TATensor> vec
                = {left_t[0], mid.gemm(right_t[1], factor, mid_gh)};
            return std::make_tuple(out_range, std::move(vec),
                                   std::min(lcut, rcut));
        } else { // Left_t[0] was larger than right_t[1]
            auto mid_gh = math::GemmHelper(
                gh.left_op(), gh.right_op(), left_t[0].range().dim(),
                left_t[0].range().dim(), mid.range().dim());

            std::vector<TATensor> vec
                = {left_t[0].gemm(mid, factor, mid_gh), right_t[1]};
            return std::make_tuple(out_range, std::move(vec),
                                   std::min(lcut, rcut));
        }
    } else if (gh.right_op() == madness::cblas::CBLAS_TRANSPOSE::Trans) {
        assert(false); // I don't have this case yet or know exactly how to
                       // handle it.
    }
}

} // namespace detail
} // namespace tensor
} // namespace tcc
