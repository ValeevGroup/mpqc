#ifndef FULL_RANK_TILE_H
#define FULL_RANK_TILE_H

#include "../include/tiledarray.h"
#include "../include/eigen.h"
#include <type_traits>

template <typename T, typename = typename std::enable_if
                      <std::is_arithmetic<T>::value, T>::type>
class FullRankTile {
  public:
    template <typename U>
    using Matrix = Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>;

    using numeric_type = T;

    FullRankTile() = default;
    FullRankTile(FullRankTile<T> const &t) = default;
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
    }

    FullRankTile(Matrix<T> t) : tile_(), rank_(std::min(t.rows(), t.cols())) {
        tile_.swap(t);
    }

    inline Matrix<T> const &data() const { return tile_; }
    inline unsigned long rank() const { return rank_; }

  private:
    Matrix<T> tile_;
    unsigned long rank_ = 0;
};

#endif // FULL_RANK_TILE_H
