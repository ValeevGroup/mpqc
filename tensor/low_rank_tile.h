#ifndef LOW_RANK_TILE_H
#define LOW_RANK_TILE_H

#include "../include/eigen.h"
#include "tile_algebra.h"
#include "../include/tiledarray.h"

namespace tcc {
namespace tensor {

template <typename T>
class LowRankTile {
  public:
    using Matrix = RowMatrix<T>;
    using scaler_type = T;

    /*
     * Constructors
     */
  public:
    LowRankTile() = default;
    LowRankTile(LowRankTile const &) = default;
    LowRankTile &operator=(LowRankTile const &t) = default;

    LowRankTile(LowRankTile &&t) noexcept : L_(), R_(), zero_{t.zero_} {
        L_.swap(t.L_);
        R_.swap(t.R_);
    }
    LowRankTile &operator=(LowRankTile &&t) noexcept {
        L_.swap(t.L_);
        R_.swap(t.R_);
        zero_ = t.zero_;
        return *this;
    }


    explicit LowRankTile(bool zero) : L_(), R_(), zero_{zero} {
        assert(zero == true);
    }

    LowRankTile(const Matrix &L, const Matrix &R)
        : L_(L.rows(), L.cols()), R_(R.rows(), R.cols()) {
        assert(L_.cols() == R_.rows());
        assert(L_.size() != 0 && R_.size() != 0);
        const auto Lsize = L.size();
        const auto Rsize = R.size();
        std::copy(L.data(), L.data() + Lsize, L_.data());
        std::copy(R.data(), R.data() + Rsize, R_.data());
    }
    LowRankTile(Matrix &&L, Matrix &&R) noexcept : L_(), R_() {
        L_.swap(L);
        R_.swap(R);
        assert(L_.cols() == R_.rows());
        assert(L_.size() != 0 && R_.size() != 0);
    }

    /*
     * Functions
     */
  public:
    inline unsigned long rank() const {
        assert(L_.cols() == R_.rows());
        return L_.cols();
    }
    inline unsigned long Rows() const { return L_.rows(); }
    inline unsigned long Cols() const { return R_.cols(); }
    inline unsigned long size() const { return R_.size() + L_.size(); }
    inline bool iszero() const { return zero_; }

    inline Matrix const &matrixL() const { return L_; }
    inline Matrix const &matrixR() const { return R_; }
    inline Matrix &matrixL() { return L_; }
    inline Matrix &matrixR() { return R_; }
    inline Matrix matrix() const { return algebra::cblas_gemm(L_, R_, 1.0); }

    template <typename Archive>
    typename madness::enable_if<madness::archive::is_output_archive<Archive>>::
        type
        serialize(Archive &ar) {
        ar &zero_ &L_.rows() & L_.cols() & R_.rows() & R_.cols()
            & madness::archive::wrap(L_.data(), L_.size())
            & madness::archive::wrap(R_.data(), R_.size());
    }

    template <typename Archive>
    typename madness::enable_if<madness::archive::is_input_archive<Archive>>::
        type
        serialize(Archive &ar) {
        bool zero;
        ar &zero;
        decltype(L_.rows()) Lrows, Lcols, Rrows, Rcols;
        ar &Lrows &Lcols &Rrows &Rcols;
        Matrix L(Lrows, Lcols), R(Rrows, Rcols);
        ar &madness::archive::wrap(L.data(), L.size())
            & madness::archive::wrap(R.data(), R.size());

        zero_ = zero;
        L_.swap(L);
        R_.swap(R);
    }

  private:
    Matrix L_;
    Matrix R_;
    bool zero_ = false;
};

} // namespace tensor
} // namespace tcc

#endif // LOW_RANK_TILE_H
