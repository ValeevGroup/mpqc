#ifndef TCC_MATRIX_TILE_MUTATIONS_H
#define TCC_MATRIX_TILE_MUTATIONS_H

#include "tile_algebra.h"
#include "tile_variant.h"


namespace unary_mutations {

template <typename T>
FullRankTile<T> scale(FullRankTile<T> &&t, T factor) {
    t.matrix() *= factor;
    return t;
}

template <typename T>
LowRankTile<T> scale(LowRankTile<T> &&t, T factor) {
    t.matrixL() *= factor;
    return t;
}

struct scale_functor {
    double factor = 1.0;
    scale_functor(double s) : factor(s) {}

    template <typename T>
    TileVariant<T> operator()(FullRankTile<T> &&t) {
        return TileVariant<T>{scale(std::move(t), factor)};
    }

    template <typename T>
    TileVariant<T> operator()(LowRankTile<T> &&t) {
        return TileVariant<T>{scale(std::move(t), factor)};
    }
};


template <typename T>
void compress(LowRankTile<T> &t, double cut) {
    if (t.matrixL().size() < t.matrixR().size()) {
        algebra::CompressLeft(t.matrixL(), t.matrixR(), cut);
    } else {
        algebra::CompressRight(t.matrixL(), t.matrixR(), cut);
    }
}

struct compress_functor {
    double cut = 1e-7;
    compress_functor() = default;
    compress_functor(double c) : cut(c) {}

    template <typename T>
    TileVariant<T> operator()(FullRankTile<T> &&t) const {
        typename FullRankTile<T>::template Matrix<T> L, R;

        if (!algebra::Decompose_Matrix(t.matrix(), L, R, cut)) {
            return TileVariant<T>{LowRankTile<T>{std::move(L), std::move(R)}};
        } else {
            return TileVariant<T>{std::move(t)};
        }
    }

    template <typename T>
    TileVariant<T> operator()(LowRankTile<T> &&t) const {
        compress(t, cut);
        return TileVariant<T>{std::move(t)};
    }
};

} // namespace unary mutations

namespace binary_mutations {

template <typename T>
FullRankTile<T>
subt(FullRankTile<T> &&left, FullRankTile<T> const &right, double factor) {
    left.matrix() -= factor * right.matrix();
    return left;
}

template <typename T>
FullRankTile<T>
subt(FullRankTile<T> &&left, LowRankTile<T> const &right, double factor) {
    left.matrix() -= factor * right.matrix();
    return left;
}

template <typename T>
FullRankTile<T>
subt(LowRankTile<T> &&left, FullRankTile<T> const &right, double factor) {
    return FullRankTile<T>{left.matrix() - factor * right.matrix()};
}

template <typename T>
LowRankTile<T>
subt(LowRankTile<T> &&left, LowRankTile<T> const &right, double factor) {
    const auto rows = left.Rows();
    const auto cols = left.Cols();
    const auto rank_out = left.rank() + right.rank();

    using matrix = typename LowRankTile<T>::template Matrix<T>;

    auto L = matrix{rows, rank_out};
    L.leftCols(left.rank()) = left.matrixL();
    L.rightCols(right.rank()) = -factor * right.matrixL();

    auto R = matrix{rank_out, cols};
    R.topRows(left.rank()) = left.matrixR();
    R.bottomRows(right.rank()) = right.matrixR();

    left.matrixL().swap(L);
    left.matrixR().swap(R);

    return left;
}

struct subt_functor {
    double factor = 1.0;
    subt_functor(double f) : factor{f} {}

    template <typename Result, typename Right>
    TileVariant<typename Result::scaler_type>
    operator()(Result &&r, // Warning forwarding ref, make sure works
               Right const &right) {
        return TileVariant<typename Result::scaler_type>{
            subt(std::move(r), right, factor)};
    }
};

} // namespace binary_mutations

namespace ternary_mutations {

template <typename T>
FullRankTile<T> gemm(FullRankTile<T> &&result, FullRankTile<T> const &left,
                     FullRankTile<T> const &right, double alpha, double beta) {
    algebra::cblas_gemm_inplace(left.matrix(), right.matrix(), result.matrix(),
                                alpha, beta);
    return result;
}

template <typename T>
FullRankTile<T> gemm(FullRankTile<T> &&result, FullRankTile<T> const &left,
                     LowRankTile<T> const &right, double alpha, double beta) {
    auto temp = algebra::cblas_gemm(left.matrix(), right.matrixL(), 1.0);
    algebra::cblas_gemm_inplace(temp, right.matrixR(), result.matrix(), alpha,
                                beta);
    return result;
}

template <typename T>
FullRankTile<T> gemm(FullRankTile<T> &&result, LowRankTile<T> const &left,
                     FullRankTile<T> const &right, double alpha, double beta) {
    auto temp = algebra::cblas_gemm(left.matrixR(), right.matrix(), 1.0);
    algebra::cblas_gemm_inplace(left.matrixL(), temp, result.matrix(), alpha,
                                beta);
    return result;
}

template <typename T> // Tricky case assume this is never called.
FullRankTile<T> gemm(LowRankTile<T> &&result, FullRankTile<T> const &left,
                     FullRankTile<T> const &right, double alpha, double beta) {
    auto out = algebra::cblas_gemm(result.matrixL(), result.matrixR(), beta);
    algebra::cblas_gemm_inplace(left.matrix(), right.matrix(), out, alpha, 1.0);
    return out;
}

template <typename T>
FullRankTile<T> gemm(FullRankTile<T> &&result, LowRankTile<T> const &left,
                     LowRankTile<T> const &right, double alpha, double beta) {

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
LowRankTile<T> gemm(LowRankTile<T> &&result, LowRankTile<T> const &left,
                    FullRankTile<T> const &right, double alpha, double beta) {
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
LowRankTile<T> gemm(LowRankTile<T> &&result, FullRankTile<T> const &left,
                    LowRankTile<T> const &right, double alpha, double beta) {
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
LowRankTile<T> gemm(LowRankTile<T> &&result, LowRankTile<T> const &left,
                    LowRankTile<T> const &right, double alpha, double beta) {
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

    return LowRankTile<T>{std::move(L), std::move(R)};
}


struct gemm_functor {
  private:
    double alpha = 1.0, beta = 1.0;

  public:
    gemm_functor(double a) : alpha(a) {}
    gemm_functor(double a, double b) : alpha(a), beta(b) {}

    template <typename Result, typename Left, typename Right>
    TileVariant<typename Result::scaler_type>
    operator()(Result &&result, // Fowarding Ref Check that works.
               Left const &left, Right const &right) const {
        return TileVariant<typename Result::scaler_type>{
            gemm(std::move(result), left, right, alpha, beta)};
    }
};

} // namespace ternary_muations

#endif // TCC_MATRIX_TILE_MUTATIONS_H
