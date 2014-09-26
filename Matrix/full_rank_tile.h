#ifndef FULL_RANK_TILE_H
#define FULL_RANK_TILE_H

#include "../include/tiledarray.h"
#include "../include/eigen.h"
#include "tile_algebra.h"
#include <type_traits>

template <typename T>
class FullRankTile {
  public:
    template <typename U>
    using Matrix = Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>;

    using scaler_type = T;

    FullRankTile() = default;
    FullRankTile(FullRankTile<T> const &t) = default;
    FullRankTile &operator=(FullRankTile<T> const &t) {
        tile_.resize(t.Rows(), t.Cols());
        std::copy(t.data(), t.data() + t.size(), tile_.data());
        return *this;
    }
    FullRankTile(FullRankTile<T> &&t) noexcept : tile_() {
        tile_.swap(t.tile_);
    }
    FullRankTile &operator=(FullRankTile<T> &&t) noexcept {
        tile_.swap(t.tile_);
        return *this;
    }

    FullRankTile(Matrix<T> const &t) : tile_(t.rows(), t.cols()) {
        std::copy(t.data(), t.data() + t.size(), tile_.data());
    }

    FullRankTile(Matrix<T> &&t) : tile_() { tile_.swap(t); }

    inline Matrix<T> const &matrix() const { return tile_; }
    inline Matrix<T> &matrix() { return tile_; }

    inline T const *data() const { return tile_.data(); }
    inline T *data() { return tile_.data(); }

    inline unsigned long rank() const {
        return std::min(tile_.rows(), tile_.cols());
    }
    inline unsigned long Rows() const { return tile_.rows(); }
    inline unsigned long Cols() const { return tile_.cols(); }
    inline unsigned long size() const { return tile_.size(); }

  private:
    Matrix<T> tile_;
};

#endif // FULL_RANK_TILE_H
