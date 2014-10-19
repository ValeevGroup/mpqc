#ifndef FULL_RANK_TILE_H
#define FULL_RANK_TILE_H

#include "../include/tiledarray.h"
#include "../include/eigen.h"
#include "tile_algebra.h"
#include <type_traits>

template <typename T, typename = typename std::enable_if
                      <std::is_arithmetic<T>::value, T>::type>
class FullRankTile {
  public:
    template <typename U>
    using Matrix = Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>;

    using scaler_type = T;

    FullRankTile() = default;
    FullRankTile(FullRankTile const &t) = default;
    FullRankTile &operator=(FullRankTile<T> const &t) {
        tile_ = t.tile_;
        rank_ = t.rank_;
        return *this;
    }
    FullRankTile(FullRankTile<T> &&t) noexcept : tile_(),
                                                 rank_(std::move(t.rank_)) {
        tile_.swap(t.tile_);
    }
    FullRankTile &operator=(FullRankTile<T> &&t) noexcept {
        rank_ = std::move(t.rank_);
        tile_.swap(t.tile_);
        return *this;
    }

    FullRankTile(Matrix<T> t) : tile_(), rank_(std::min(t.rows(), t.cols())) {
        tile_.swap(t);
    }

    inline Matrix<T> const &matrix() const { return tile_; }
    inline Matrix<T> &matrix() { return tile_; }

    inline T const *data() const { return tile_.data(); }
    inline T *data() { return tile_.data(); }

    inline unsigned long rank() const { return rank_; }
    inline unsigned long Rows() const { return tile_.rows(); }
    inline unsigned long Cols() const { return tile_.cols(); }
    inline unsigned long size() const { return tile_.size(); }
    inline bool iszero() const {return false;}

  private:
    Matrix<T> tile_;
    unsigned long rank_ = 0;
};

#endif // FULL_RANK_TILE_H
