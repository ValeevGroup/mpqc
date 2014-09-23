#ifndef TileClusterChem_TILE_OPS_H
#define TileClusterChem_TILE_OPS_H

#include "low_rank_tile.h"
#include "full_rank_tile.h"
#include "tile_algebra.h"
#include "tile_variant.h"

// Gemm_AB functions
namespace tile_ops {

template <typename T>
FullRankTile<T>
gemm(const FullRankTile<T> &left, const FullRankTile<T> &right, double alpha) {
    return FullRankTile<T>{
        algebra::cblas_gemm(left.matrix(), right.matrix(), alpha)};
}

template <typename T>
LowRankTile<T>
gemm(const LowRankTile<T> &left, const FullRankTile<T> &right, double alpha) {
    return LowRankTile<T>{
        (alpha * left.matrixL()).eval(),
        algebra::cblas_gemm(left.matrixR(), right.matrix(), 1.0)};
}

template <typename T>
LowRankTile<T>
gemm(const FullRankTile<T> &left, const LowRankTile<T> &right, double alpha) {
    return LowRankTile<T>{
        algebra::cblas_gemm(left.matrix(), right.matrixL(), alpha),
        right.matrixR()};
}

template <typename T>
LowRankTile<T>
gemm(const LowRankTile<T> &left, const LowRankTile<T> &right, double alpha) {
    assert(left.Cols() == right.Rows()); // Check k index

    auto mid = algebra::cblas_gemm(left.matrixR(), right.matrixL(), 1.0);
    if (left.rank() >= right.rank()) {
        typename LowRankTile<T>::template Matrix<T> R(right.rank(),
                                                      right.Cols());

        // jumping through copy hoops because it's faster.
        const auto start = right.matrixR().data();
        const auto end = right.matrixR().data() + right.matrixR().size();
        std::copy(start, end, R.data());
        return LowRankTile<T>{algebra::cblas_gemm(left.matrixL(), mid, alpha),
                              std::move(R)};
    } else {
        typename LowRankTile<T>::template Matrix<T> L(left.Rows(), left.rank());
        const auto start = left.matrixL().data();
        const auto end = left.matrixL().data() + left.matrixL().size();
        std::transform(start, end, L.data(),
                       [=](const T &x) { return alpha * x; });
        return LowRankTile<T>{std::move(L),
                              algebra::cblas_gemm(mid, right.matrixR(), 1.0)};
    }
}

struct gemm_AB {
  private:
    double alpha = 1.0;

  public:
    gemm_AB(double a) : alpha(a) {}

    template <typename Left, typename Right>
    TileVariant<typename Left::scaler_type>
    operator()(Left const &left, Right const &right) const {
        return TileVariant<typename Left::scaler_type>{gemm(left, right, alpha)};
    }
};

} // namespace tile_ops ending gemm

// Gemm_inplace_functions
namespace tile_ops {
template <typename T>
FullRankTile<T> gemm(FullRankTile<T> result, const FullRankTile<T> &left,
                     const FullRankTile<T> &right, double alpha, double beta) {
    algebra::cblas_gemm_inplace(left.matrix(), right.matrix(), result.matrix(),
                                alpha, beta);
    return result;
}

template <typename T>
FullRankTile<T> gemm(FullRankTile<T> result, const FullRankTile<T> &left,
                     const LowRankTile<T> &right, double alpha, double beta) {
    auto temp = algebra::cblas_gemm(left.matrix(), right.matrixL(), 1.0);
    algebra::cblas_gemm_inplace(temp, right.matrixR(), result.matrix(), alpha,
                                beta);
    return result;
}

template <typename T>
FullRankTile<T> gemm(FullRankTile<T> result, const LowRankTile<T> &left,
                     const FullRankTile<T> &right, double alpha, double beta) {
    auto temp = algebra::cblas_gemm(left.matrixR(), right.matrix(), 1.0);
    algebra::cblas_gemm_inplace(left.matrixL(), temp, result.matrix(), alpha,
                                beta);
    return result;
}

template <typename T> // Tricky case assume this is never called.
FullRankTile<T> gemm(LowRankTile<T> result, const FullRankTile<T> &left,
                     const FullRankTile<T> &right, double alpha, double beta) {
    auto out = algebra::cblas_gemm(result.matrixL(), result.matrixR(), beta);
    algebra::cblas_gemm_inplace(left.matrix(), right.matrix(), out, alpha, 1.0);
    return out;
}

template <typename T>
FullRankTile<T> gemm(FullRankTile<T> result, const LowRankTile<T> &left,
                     const LowRankTile<T> &right, double alpha, double beta) {

    auto mid = algebra::cblas_gemm(left.matrixR(), right.matrixL(), 1.0);

    // Doing some math to determine which direction to take.
    const auto n = left.Cols();
    const auto m = right.Rows();
    const auto ra = left.rank();
    const auto rb = right.rank();

    const auto left_cost = n * (ra * ra + m * rb);  // Figure out which
    const auto right_cost = m * (rb * rb + n * ra); // path uses fewer flops

    if (left_cost < right_cost) { // Go left
        auto temp = algebra::cblas_gemm(left.matrixL(), mid, 1.0);
        algebra::cblas_gemm_inplace(temp, right.matrixR(), result.matrix(),
                                    alpha, beta);
    } else { // Go right
        auto temp = algebra::cblas_gemm(mid, right.matrixR(), 1.0);
        algebra::cblas_gemm_inplace(left.matrixL(), temp, result.matrix(),
                                    alpha, beta);
    }

    return result;
}

template <typename T>
LowRankTile<T> gemm(LowRankTile<T> result, const LowRankTile<T> &left,
                    const FullRankTile<T> &right, double alpha, double beta) {
    const auto rows = result.Rows();
    const auto cols = result.Cols();

    const auto rank_A = left.rank();
    const auto rank_C = result.rank();
    const auto rank_out = rank_C + rank_A;

    using matrix = typename LowRankTile<T>::template Matrix<T>;

    auto L = matrix{rows, rank_out};
    L.leftCols(rank_A) = alpha * left.matrixL();
    L.rightCols(rank_C) = beta * result.matrixL();

    auto R = matrix{rank_out, cols};
    R.topRows(rank_A)
        = algebra::cblas_gemm(left.matrixR(), right.matrix(), 1.0);
    R.bottomRows(rank_C) = result.matrixR();

    result = LowRankTile<T>{std::move(L), std::move(R)};

    return result;
}

template <typename T>
LowRankTile<T> gemm(LowRankTile<T> result, const FullRankTile<T> &left,
                    const LowRankTile<T> &right, double alpha, double beta) {
    const auto rows = result.Rows();
    const auto cols = result.Cols();

    const auto rank_A = right.rank();
    const auto rank_C = result.rank();
    const auto rank_out = rank_C + rank_A;

    using matrix = typename LowRankTile<T>::template Matrix<T>;

    auto L = matrix{rows, rank_out};
    L.leftCols(rank_A)
        = algebra::cblas_gemm(left.matrix(), right.matrixL(), alpha);
    L.rightCols(rank_C) = beta * result.matrixL();

    auto R = matrix{rank_out, cols};
    R.topRows(rank_A) = right.matrixR();
    R.bottomRows(rank_C) = result.matrixR();

    result = LowRankTile<T>{std::move(L), std::move(R)};

    return result;
}

template <typename T>
LowRankTile<T> gemm(LowRankTile<T> result, const LowRankTile<T> &left,
                    const LowRankTile<T> &right, double alpha, double beta) {
    assert(left.Cols() == right.Rows()); // Check k index

    const auto rows = result.Rows();
    const auto cols = result.Cols();

    const auto rank_AB = std::min(left.rank(), right.rank());
    const auto rank_C = result.rank();
    const auto rank_out = rank_C + rank_AB;

    using matrix = typename LowRankTile<T>::template Matrix<T>;

    const auto mid = algebra::cblas_gemm(left.matrixR(), right.matrixL(), 1.0);

    auto L = matrix{rows, rank_out};
    L.rightCols(rank_C) = beta * result.matrixL();

    auto R = matrix{rank_out, cols};
    R.bottomRows(rank_C) = result.matrixR();

    if (left.rank() > right.rank()) {
        L.leftCols(rank_AB) = algebra::cblas_gemm(left.matrixL(), mid, alpha);
        R.topRows(rank_AB) = right.matrixR();
    } else {
        L.leftCols(rank_AB) = alpha * left.matrixL();
        R.topRows(rank_AB) = algebra::cblas_gemm(mid, right.matrixR(), 1.0);
    }

    result = LowRankTile<T>{std::move(L), std::move(R)};

    return result;
}


struct gemm_inplace {
  private:
    double alpha = 1.0, beta = 1.0;

  public:
    gemm_inplace(double a) : alpha(a) {}
    gemm_inplace(double a, double b) : alpha(a), beta(b) {}

    template <typename Result, typename Left, typename Right>
    TileVariant<typename Result::scaler_type>
    operator()(Result result, Left const &left, Right const &right) const {
        return TileVariant<typename Result::scaler_type>{
            gemm(result, left, right, alpha, beta)};
    }
};

} // namespace tile_ops ending gemm_inplace


//
namespace tile_ops {

template <typename T>
FullRankTile<T>
add(const FullRankTile<T> &left, const FullRankTile<T> &right, double beta) {
    return FullRankTile<T>{beta * left.matrix() + right.matrix()};
}

template <typename T>
LowRankTile<T>
add(const LowRankTile<T> &left, const LowRankTile<T> &right, double beta) {
    assert(left.Rows() == right.Rows());
    assert(left.Cols() == right.Cols());

    const auto rows = left.Rows();
    const auto cols = left.Cols();
    const auto rank_out = left.rank() + right.rank();

    using matrix = typename LowRankTile<T>::template Matrix<T>;

    auto L = matrix{rows, rank_out};
    L.leftCols(left.rank()) = beta * left.matrixL();
    L.rightCols(right.rank()) = right.matrixL();

    auto R = matrix{rank_out, cols};
    R.topRows(left.rank()) = left.matrixR();
    R.bottomRows(right.rank()) = right.matrixR();

    return LowRankTile<T>{std::move(L), std::move(R)};
}

template <typename T>
LowRankTile<T> &compress(LowRankTile<T> &tile, double cut) {
    assert(tile.L_.cols() == tile.R_.rows());
    if (tile.matrixL().size() <= tile.matrixR().size()) {
        algebra::CompressLeft(tile.L_, tile.R_, cut);
    } else {
        algebra::CompressRight(tile.L_, tile.R_, cut);
    }
    assert(tile.L_.cols() == tile.R_.rows());
    tile.rank_ = tile.L_.cols();
    return tile;
}

} // namespace tile_ops

#endif // TileClusterChem_TILE_OPS_H
