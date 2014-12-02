#ifndef FULL_RANK_TILE_H
#define FULL_RANK_TILE_H

#include "../include/tiledarray.h"
#include "../include/eigen.h"
#include "tile_algebra.h"
#include <type_traits>

namespace tcc {
    namespace tensor {

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
        return *this;
    }
    FullRankTile(FullRankTile<T> &&t) noexcept : tile_() {
        tile_.swap(t.tile_);
    }
    FullRankTile &operator=(FullRankTile<T> &&t) noexcept {
        tile_.swap(t.tile_);
        return *this;
    }

    FullRankTile(Matrix<T> t) : tile_() {
        tile_.swap(t);
    }

    inline Matrix<T> const &matrix() const { return tile_; }
    inline Matrix<T> &matrix() { return tile_; }

    inline T const *data() const { return tile_.data(); }
    inline T *data() { return tile_.data(); }

    inline unsigned long Rows() const { return tile_.rows(); }
    inline unsigned long Cols() const { return tile_.cols(); }
    inline unsigned long size() const { return tile_.size(); }
    inline unsigned long rank() const {
        return std::min(tile_.rows(), tile_.cols());
    }
    inline bool iszero() const {return false;}


    template <typename Archive>
    typename madness::enable_if<madness::archive::is_output_archive<Archive> >::type
    serialize(Archive& ar) {
        ar & tile_.rows() & tile_.cols() 
           & madness::archive::wrap(tile_.data(), tile_.size());
    }

    template <typename Archive>
    typename madness::enable_if<madness::archive::is_input_archive<Archive> >::
        type
        serialize(Archive &ar) {
        decltype(tile_.rows())  rows, cols;
        ar & rows & cols;
        Matrix<T> mat(rows, cols); 
        ar & madness::archive::wrap(mat.data(), mat.size());
        tile_.swap(mat);
    }

  private:
    Matrix<T> tile_;
};

} // namespace tensor
} // namespace tcc 

#endif // FULL_RANK_TILE_H
