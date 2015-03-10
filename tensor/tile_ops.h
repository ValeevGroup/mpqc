#ifndef TileClusterChem_TILE_OPS_H
#define TileClusterChem_TILE_OPS_H

#include "low_rank_tile.h"
#include "full_rank_tile.h"
#include "tile_algebra.h"
#include "tile_variant.h"
#include <iostream>

namespace tcc {
namespace tensor {

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

// D_gemm is for doing the DF exchange intermediate contraction with the
// density.
template <typename T>
FullRankTile<T> D_gemm(const FullRankTile<T> &left,
                       const FullRankTile<T> &right, double alpha) {
    // Reshape for contraction
    Eigen::MatrixXd L = left.matrix();
    L.resize(L.rows() * L.cols() / right.matrix().rows(),
             right.matrix().rows());

    // Contract and reshape to expected shape
    auto outL = algebra::cblas_gemm(L, right.matrix(), alpha);
    L.resize(left.matrix().rows(), left.matrix().cols());
    return FullRankTile<T>{std::move(outL)};
}

template <typename T>
LowRankTile<T>
D_gemm(const LowRankTile<T> &left, const FullRankTile<T> &right, double alpha) {
    // Compute contraction dim and reshape
    const auto K = right.matrix().rows();
    Eigen::MatrixXd leftR = left.matrixR();
    leftR.resize(leftR.rows() + leftR.cols() / K, K);

    auto out_leftR = algebra::cblas_gemm(leftR, right.matrix(), alpha);
    out_leftR.resize(left.matrixR().rows(), left.matrixR().cols());

    return LowRankTile<T>{left.matrixL(), std::move(out_leftR)};
}

template <typename T>
FullRankTile<T>
D_gemm(const FullRankTile<T> &left, const LowRankTile<T> &right, double alpha) {
    const auto K = right.matrixL().rows();
    Eigen::MatrixXd L = left.matrix();
    L.resize(L.rows() * L.cols() / K, K);

    auto tmpL = algebra::cblas_gemm(L, right.matrixL(), alpha);
    L = algebra::cblas_gemm(tmpL, right.matrixR(), 1.0);
    L.resize(left.matrix().rows(), left.matrix().cols());

    return FullRankTile<T>{std::move(L)};
}

template <typename T>
LowRankTile<T>
D_gemm(const LowRankTile<T> &left, const LowRankTile<T> &right, double alpha) {

    const auto K = right.matrixL().rows();
    Eigen::MatrixXd leftR = left.matrixR();
    leftR.resize(leftR.rows() * leftR.cols() / K, K);

    auto mid = algebra::cblas_gemm(leftR, right.matrixL(), 1.0);
    Eigen::MatrixXd R = algebra::cblas_gemm(mid, right.matrixR(), 1.0);
    R.resize(left.rank(), right.matrixR().cols());

    return LowRankTile<T>{(alpha * left.matrixL()).eval(), std::move(R)};
}

struct D_functor {
  private:
    double alpha = 1.0;

  public:
    D_functor(double a) : alpha(a) {}

    template <typename Left, typename Right>
    TileVariant<typename Left::scaler_type>
    operator()(Left const &left, Right const &right) const {
        return TileVariant<typename Left::scaler_type>{
            D_gemm(left, right, alpha)};
    }
};

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
        typename LowRankTile<T>::Matrix R(right.rank(), right.Cols());

        // jumping through copy hoops because it's faster.
        const auto start = right.matrixR().data();
        const auto end = right.matrixR().data() + right.matrixR().size();
        std::copy(start, end, R.data());
        return LowRankTile<T>{algebra::cblas_gemm(left.matrixL(), mid, alpha),
                              std::move(R)};
    } else {
        typename LowRankTile<T>::Matrix L(left.Rows(), left.rank());
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
Dgemm(const FullRankTile<T> &left, const FullRankTile<T> &right, double alpha) {
    RowMatrixXd L = left.matrix();

    // Reshape L for contraction
    const auto K = right.matrix().rows();
    const auto other = L.cols() / K;
    const auto L_rows_out = L.rows();
    const auto L_cols_out = other * right.matrix().cols();
    L.resize(L.rows() * other, K);

    L = algebra::cblas_gemm(L, right.matrix(), alpha);
    L.resize(L_rows_out, L_cols_out);

    return FullRankTile<T>{std::move(L)};
}

template <typename T>
LowRankTile<T>
Dgemm(const LowRankTile<T> &left, const FullRankTile<T> &right, double alpha) {

    RowMatrixXd Lr = left.matrixR();

    // Reshape Lr for contraction
    const auto K = right.matrix().rows();
    const auto other = Lr.cols() / K;
    const auto Lr_rows_out = Lr.rows();
    const auto Lr_cols_out = other * right.matrix().cols();
    Lr.resize(Lr.rows() * other, K);

    // Contract and reshape back to original
    Lr = algebra::cblas_gemm(Lr, right.matrix(), 1.0);
    Lr.resize(Lr_rows_out, Lr_cols_out);

    return LowRankTile<T>{(alpha * left.matrixL()).eval(), std::move(Lr)};
}

template <typename T>
FullRankTile<T>
Dgemm(const FullRankTile<T> &left, const LowRankTile<T> &right, double alpha) {
    RowMatrixXd L = left.matrix();

    // Reshape L for contraction
    const auto K = right.matrixL().rows();
    const auto other = L.cols() / K;
    const auto L_rows_out = L.rows();
    const auto L_cols_out = other * right.matrix().cols();
    L.resize(L.rows() * other, K);

    // Contract and reshape back to original
    L = algebra::cblas_gemm(L, right.matrixL(), alpha);

    // This is the odd guy out in the D contraction because it is
    // (Xi, j) * (j,r2) * (r2, a) So I don't have an easy way to make
    // the output (X,ia) which for now I depend on so make it full rank
    L = algebra::cblas_gemm(L, right.matrixR(), 1.0);
    L.resize(L_rows_out, L_cols_out);

    return FullRankTile<T>{std::move(L)};
}

template <typename T>
LowRankTile<T>
Dgemm(const LowRankTile<T> &left, const LowRankTile<T> &right, double alpha) {

    RowMatrixXd Lr = left.matrixR();
    // Reshape Lr for contraction
    const auto K = right.matrixR().cols();
    const auto other = Lr.cols() / K;
    const auto Lr_rows_out = Lr.rows();
    const auto Lr_cols_out = other * right.matrixR().cols();
    Lr.resize(Lr.rows() * other, K);

    auto mid = algebra::cblas_gemm(Lr, right.matrixL(), 1.0);
    Lr = algebra::cblas_gemm(mid, right.matrixR(), 1.0);
    Lr.resize(Lr_rows_out, Lr_cols_out);

    return LowRankTile<T>{(alpha * left.matrixL()).eval(), std::move(Lr)};
}

struct Dgemm_functor {
  private:
    double alpha = 1.0;

  public:
    Dgemm_functor(double a) : alpha(a) {}

    template <typename Left, typename Right>
    TileVariant<typename Left::scaler_type>
    operator()(Left const &left, Right const &right) const {
        return TileVariant<typename Left::scaler_type>{
            Dgemm(left, right, alpha)};
    }
};

template <typename T>
FullRankTile<T>
Xgemm(const FullRankTile<T> &left, const FullRankTile<T> &right, double alpha) {
    
    RowMatrixXd L = left.matrix();
    RowMatrixXd R = right.matrix();

    // (X,ai) => (Xa,i) => (i, Xa) * (Xa,j) = (i,j)
    const int bs_dim = std::sqrt(R.cols());
    const auto K = R.size() / bs_dim;
    L.resize(K, bs_dim);
    RowMatrixXd Lt = L.transpose();
    R.resize(K, bs_dim);

    return FullRankTile<T>{algebra::cblas_gemm(Lt, R, alpha)};
}

template <typename T>
FullRankTile<T>
Xgemm(const LowRankTile<T> &left, const FullRankTile<T> &right, double alpha) {

    RowMatrixXd Lr = left.matrixR();
    RowMatrixXd LlT = left.matrixL().transpose();

    // Reshape for contraction
    const int bs_dim = std::sqrt(right.matrix().cols());

    // (r1,X) * (X,aj) => (r1, aj)
    auto mid = algebra::cblas_gemm(LlT, right.matrix(), 1.0);

    // (r1, ai) => (r1a, i) => (i, r1a)
    Lr.resize(Lr.rows() * bs_dim, bs_dim);
    RowMatrixXd LrT = Lr.transpose();

    // (r1, aj) => (r1a, j)
    mid.resize(mid.rows() * bs_dim, bs_dim);

    // (i, r1a) * (r1a, j)
    return FullRankTile<T>{algebra::cblas_gemm(LrT, mid, alpha)};
}

template <typename T>
FullRankTile<T>
Xgemm(const FullRankTile<T> &left, const LowRankTile<T> &right, double alpha) {

    RowMatrixXd Rr = right.matrixR();

    // (ai, X) * (X,r2) => (ai,r2)
    RowMatrixXd L = left.matrix().transpose();
    auto mid
        = algebra::cblas_gemm(L, right.matrixL(), 1.0);

    // (ai, r2) => (a, ir2) => (ir2, a) => (i, r2a)
    const int bs_dim = std::sqrt(Rr.cols());
    mid.resize(bs_dim, mid.cols() * bs_dim);
    RowMatrixXd midT = mid.transpose();
    midT.resize(bs_dim, mid.size() / bs_dim);

    // (r2, aj) => (r2a,j)
    Rr.resize(Rr.size() / bs_dim, bs_dim);

    // (i, r2a) * (r2a, j)
    return FullRankTile<T>{algebra::cblas_gemm(midT, Rr, alpha)};
}

template <typename T>
FullRankTile<T>
Xgemm(const LowRankTile<T> &left, const LowRankTile<T> &right, double alpha) {

    // (r1, X) * (X,r2) => (r1,r2)
    RowMatrixXd Ll = left.matrixL().transpose();
    auto mid
        = algebra::cblas_gemm(Ll, right.matrixL(), 1.0);

    // (r1,r2) * (r2, aj) => (r1, aj) => (r1a,j)
    mid = algebra::cblas_gemm(mid, right.matrixR(), 1.0);
    const auto bs_dim = std::sqrt(right.matrixR().cols());
    mid.resize(mid.size()/bs_dim, bs_dim);

    //(r1, ai) => (r1a, i) => (i, r1a)
    RowMatrixXd Lr = left.matrixR();
    Lr.resize(bs_dim, Lr.size()/bs_dim);
    RowMatrixXd LrT = Lr.transpose();

    // (i, r1a) * r1a, j)
    return FullRankTile<T>{algebra::cblas_gemm(LrT, mid, alpha)};
}

struct Xgemm_functor {
  private:
    double alpha = 1.0;

  public:
    Xgemm_functor(double a) : alpha(a) {}

    template <typename Left, typename Right>
    TileVariant<typename Left::scaler_type>
    operator()(Left const &left, Right const &right) const {
        return TileVariant<typename Left::scaler_type>{
            Xgemm(left, right, alpha)};
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

    using matrix = typename LowRankTile<T>::Matrix;

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

    using matrix = typename LowRankTile<T>::Matrix;

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

} // namespace tensor
} // namespace tcc

#endif // TileClusterChem_TILE_OPS_H
