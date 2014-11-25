#ifndef TileClusterChem_TILE_OPS_H
#define TileClusterChem_TILE_OPS_H

#include "low_rank_tile.h"
#include "full_rank_tile.h"
#include "tile_algebra.h"
#include "tile_variant.h"

namespace unary_ops {
template <typename T>
FullRankTile<T> scale(FullRankTile<T> const &t, T factor) {
    return FullRankTile<T>{factor * t.matrix()};
}

template <typename T>
LowRankTile<T> scale(LowRankTile<T> const &t, T factor) {
    return LowRankTile<T>{factor * t.matrixL(), t.matrixR()};
}

struct scale_functor {
    double factor = 1.0;
    scale_functor(double f) : factor{f} {}

    template <typename T>
    TileVariant<T> operator()(FullRankTile<T> const &t) {
        return TileVariant<T>{scale(t, factor)};
    }

    template <typename T>
    TileVariant<T> operator()(LowRankTile<T> const &t) {
        return TileVariant<T>{scale(t, factor)};
    }
};

} // namespace unary_ops

namespace binary_ops {

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

struct gemm_functor {
  private:
    double alpha = 1.0;

  public:
    gemm_functor(double a) : alpha(a) {}

    template <typename Left, typename Right>
    TileVariant<typename Left::scaler_type>
    operator()(Left const &left, Right const &right) const {
        return TileVariant<typename Left::scaler_type>{
            gemm(left, right, alpha)};
    }
};

template <typename T>
FullRankTile<T>
add(FullRankTile<T> const &left, FullRankTile<T> const &right, double beta) {
    return FullRankTile<T>{beta * left.matrix() + right.matrix()};
}

template <typename T>
FullRankTile<T>
add(FullRankTile<T> const &left, LowRankTile<T> const &right, double beta) {
    return FullRankTile<T>{beta * left.matrix() + right.matrix()};
}

template <typename T>
FullRankTile<T>
add(LowRankTile<T> const &left, FullRankTile<T> const &right, double beta) {
    return FullRankTile<T>{beta * left.matrix() + right.matrix()};
}

template <typename T>
LowRankTile<T>
add(LowRankTile<T> const &left, LowRankTile<T> const &right, double beta) {
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

struct add_functor {
    double beta = 1.0;
    add_functor(double b) : beta(b) {}

    template <typename Left, typename Right>
    TileVariant<typename Left::scaler_type>
    operator()(Left const &left, Right const &right) const {
        return TileVariant<typename Left::scaler_type>{add(left, right, beta)};
    }
};

template <typename T>
FullRankTile<T>
subt(FullRankTile<T> const &left, FullRankTile<T> const &right, double beta) {
    return FullRankTile<T>{beta * left.matrix() - right.matrix()};
}

template <typename T>
FullRankTile<T>
subt(FullRankTile<T> const &left, LowRankTile<T> const &right, double beta) {
    return FullRankTile<T>{beta * left.matrix() - right.matrix()};
}

template <typename T>
FullRankTile<T>
subt(LowRankTile<T> const &left, FullRankTile<T> const &right, double beta) {
    return FullRankTile<T>{beta * left.matrix() - right.matrix()};
}

template <typename T>
LowRankTile<T>
subt(LowRankTile<T> const &left, LowRankTile<T> const &right, double beta) {
    assert(left.Rows() == right.Rows());
    assert(left.Cols() == right.Cols());

    const auto rows = left.Rows();
    const auto cols = left.Cols();
    const auto rank_out = left.rank() + right.rank();

    using matrix = typename LowRankTile<T>::template Matrix<T>;

    auto L = matrix{rows, rank_out};
    L.leftCols(left.rank()) = beta * left.matrixL();
    L.rightCols(right.rank()) = -right.matrixL();

    auto R = matrix{rank_out, cols};
    R.topRows(left.rank()) = left.matrixR();
    R.bottomRows(right.rank()) = right.matrixR();

    return LowRankTile<T>{std::move(L), std::move(R)};
}

struct subt_functor {
    double beta = 1.0;
    subt_functor(double b) : beta(b) {}

    template <typename Left, typename Right>
    TileVariant<typename Left::scaler_type>
    operator()(Left const &left, Right const &right) const {
        return TileVariant<typename Left::scaler_type>{subt(left, right, beta)};
    }
};

} // namespace binary_ops

#endif // TileClusterChem_TILE_OPS_H
