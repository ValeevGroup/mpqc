#ifndef LOW_RANK_TILE_H
#define LOW_RANK_TILE_H

#include "../include/eigen.h"
#include "tile_algebra.h"

template <typename T>
class LowRankTile {
  public:
    template <typename U>
    using Matrix = Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>;

    using scaler_type = T;

    /*
     * Constructors
     */
  public:
    LowRankTile() = default;
    LowRankTile(LowRankTile const &) = default;
    LowRankTile &operator=(LowRankTile const &t) {
      L_.resize(t.L_.rows(),t.L_.cols());
      R_.resize(t.R_.rows(),t.R_.cols());
      eigen_copy(t.L_, t.R_);
    }

    LowRankTile(LowRankTile &&t) noexcept : L_(), R_() {
        L_.swap(t.L_);
        R_.swap(t.R_);
    }
    LowRankTile &operator=(LowRankTile &&t) noexcept {
        L_.swap(t.L_);
        R_.swap(t.R_);
        return *this;
    }

    LowRankTile(const Matrix<T> &L, const Matrix<T> &R)
        : L_(L.rows(), L.cols()), R_(R.rows(), R.cols()) {
      eigen_copy(L,R);
    }
    LowRankTile(Matrix<T> &&L, Matrix<T> &&R) noexcept : L_(), R_() {
        L_.swap(L);
        R_.swap(R);
    }

    /*
     * Functions
     */
  public:
    inline unsigned long rank() const { return L_.cols(); }
    inline unsigned long Rows() const { return L_.rows(); }
    inline unsigned long Cols() const { return R_.cols(); }
    inline unsigned long size() const { return R_.size() + L_.size(); }

    inline Matrix<T> const &matrixL() const { return L_; }
    inline Matrix<T> const &matrixR() const { return R_; }
    inline Matrix<T> &matrixL() { return L_; }
    inline Matrix<T> &matrixR() { return R_; }
    inline Matrix<T> matrix() const { return algebra::cblas_gemm(L_, R_, 1.0); }

  private:
    Matrix<T> L_;
    Matrix<T> R_;

    void inline eigen_copy(Matrix<T> const &L, Matrix<T> const &R){
      std::copy(L.data(), L.data()+L.size(), L_.data());
      std::copy(R.data(), R.data()+R.size(), R_.data());
    }
};

#endif // LOW_RANK_TILE_H
