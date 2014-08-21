#ifndef LR_TILE_H
#define LR_TILE_H

#include "../include/eigen.h"
#include "../include/tiledarray.h"
#include <utility>
#include <iostream>

namespace detail {

template <typename T>
using EigMat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T> using LR_pair = std::pair<EigMat<T>, EigMat<T>>;

/**
* qr_decomp returns the low rank Q and R matrices as a pair.
*/
template <typename T> LR_pair<T> qr_decomp(const EigMat<T> &input, double cut) {
  Eigen::ColPivHouseholderQR<EigMat<T>> qr(input);
  qr.setThreshold(cut);

  return std::make_pair(EigMat<T>(qr.householderQ()).leftCols(qr.rank()),
                        EigMat<T>(qr.matrixR()
                                      .topLeftCorner(qr.rank(), input.cols())
                                      .template triangularView<Eigen::Upper>())
                        * qr.colsPermutation().transpose());
}

} // namespace detail

/**
 * LRTile is a class which can be used as a low rank tile for TiledArray
 * matrices
 */
template <typename T> class LRTile {
public:
  typedef LRTile eval_type; /// The data type of the tile, typically double.
  typedef T value_type;     /// The data type of the tile, typically double.
  typedef TiledArray::Range range_type; /// Tensor range type
  typedef double numeric_type;
  typedef std::size_t size_type; /// The data type of the tile, typically
                                 /// double.

  template <typename U> using EigMat = detail::EigMat<U>;

  LRTile() = default;
  ~LRTile() = default;
  LRTile(const LRTile &rhs) = default;
  // LRTile(LRTile &&rhs) = default; Must write by hand :P.

  /**
   * @brief LRTile constructor builds a low rank representation of a matrix
   * given an input matrix.
   * @param input is an eigen matrix which with matching data type to the class.
   */
  explicit LRTile(const EigMat<T> &input, double cut = 1e-08) : L_(), R_() {
    auto QR_pair = detail::qr_decomp(input, cut);
    L_ = std::get<0>(QR_pair);
    R_ = std::get<1>(QR_pair);
    rank_ = L_.cols();
  }

  /**
   * @brief LRTile constructor to assign the low rank matrix directly
   * @param left_mat the left hand matrix
   * @param right_mat the right hand matrix
   */
  explicit LRTile(const EigMat<T> &left_mat, const EigMat<T> &right_mat)
      : L_(left_mat), R_(right_mat), rank_(right_mat.cols()) {}

  /**
   * @brief mult multiples two tiles together.
   * @param a Low Rank tile to multiple the current tile with.
   * @return a new Low Rank tile.
   */
  LRTile mult(const LRTile<T> &right) const {
    // Check which rank is smaller
    bool use_left_rank = (rank() < right.rank());
    EigMat<T> L;
    EigMat<T> R;
    if (use_left_rank) {
      L = matrixL();
      R = (matrixR() * right.matrixL()) * right.matrixR();
    } else {
      L = matrixL() * (matrixR() * right.matrixL());
      R = right.matrixR();
    }
    return LRTile(L, R);
  }

  /**
   * @brief add adds two low rank tiles togather
   * @param right the right side low rank being added to this
   * @return new low rank tile with rank this->rank() + right.rank()
   * @warning the output tile will have higher rank than the input tiles.
   */
  LRTile add(const LRTile<T> right) const {
    // output matrices
    EigMat<T> L(matrixL().rows(), rank() + right.rank());
    EigMat<T> R(rank() + right.rank(), right.matrixR().cols());

    // Fill the matrices
    L.leftCols(rank()) = matrixL();
    L.rightCols(right.rank()) = right.matrixL();

    R.topRows(rank()) = matrixR();
    R.bottomRows(right.rank()) = right.matrixR();

    LRTile out(L, R);

    // TODO need heuristic to decide when to compress.
    out.compress();

    return out;
  }

  bool empty() const { return false; }

  LRTile add(const LRTile &right, const TiledArray::Permutation &perm) const {}
  LRTile add(const LRTile &right,
             const numeric_type factor,
             const TiledArray::Permutation &perm) const {}
  LRTile
  add(const value_type &value, const TiledArray::Permutation &perm) const {}


  LRTile clone() const {}
  LRTile permute(const TiledArray::Permutation &perm) const {}
  LRTile scale(const numeric_type factor) const {}
  LRTile
  scale(const numeric_type factor, const TiledArray::Permutation &perm) const {}
  LRTile &scale_to(const numeric_type factor) const {}
  LRTile &add_to(const LRTile &right) const {}
  LRTile &add_to(const LRTile &right, const numeric_type factor) const {}
  LRTile &add_to(const value_type &value) const {}

  LRTile mult(const LRTile &right, const numeric_type factor) const {}
  LRTile mult(const LRTile &right, const TiledArray::Permutation &perm) const {}
  LRTile mult(const LRTile &right,
              const numeric_type factor,
              const TiledArray::Permutation &perm) const {}

  LRTile &mult_to(const LRTile &right) {}
  LRTile &mult_to(const LRTile &right, const numeric_type factor) {}

  /**
   * @brief compress attempts to reduce the rank of the tile
   * @param cut is the threshold parameter for Eigen::ColPivHouseholder.
   */
  void compress(double cut = 1e-8) {
    // TODO_PAR make these tasks which can run seperately.
    auto qrL = detail::qr_decomp(matrixL(), cut);
    auto qrR = detail::qr_decomp(matrixR(), cut);
    // TODO_PAR don't forget to fence since matrixR/L are refs

    bool use_left_rank = (std::get<0>(qrL).cols() < std::get<0>(qrR).cols());
    if (use_left_rank) {
      L_ = std::get<0>(qrL);
      R_ = (std::get<1>(qrL) * std::get<0>(qrR)) * std::get<1>(qrR);
    } else {
      L_ = std::get<0>(qrL) * (std::get<1>(qrL) * std::get<0>(qrR));
      R_ = std::get<1>(qrR);
    }
  }

  LRTile gemm(const LRTile &arg1,
              const LRTile::numeric_type factor,
              const TiledArray::math::GemmHelper &gemm_config) {}

  // GEMM operation with fused indices as defined by gemm_config; multiply arg1
  // by arg2, return the result
  LRTile gemm(const LRTile &arg1,
              const LRTile &arg2,
              const LRTile::numeric_type factor,
              const TiledArray::math::GemmHelper &gemm_config) {}

  // GEMM operation with fused indices as defined by gemm_config; multiply left
  // by right, store to result
  void gemm(LRTile &result,
            const LRTile &arg1,
            const LRTile &arg2,
            const LRTile::numeric_type factor,
            const TiledArray::math::GemmHelper &gemm_config) {}


  /**
   * @brief matrixL
   * @return the left hand side low rank matrix.
   */
  inline const EigMat<T> &matrixL() const { return L_; }

  /**
   * @brief matrixR
   * @return the right hand side low rank matrix.
   */
  inline const EigMat<T> &matrixR() const { return R_; }


  /**
   * @brief matrixLR
   * @return the full sized product of the low rank matrices.
   * @warning Can possibly increase memory usage by a large amount.
   */
  inline EigMat<T> matrixLR() const { return L_ * R_; }

  /**
   * @brief rank returns the rank of the tile by value.
   */
  inline std::size_t rank() const { return rank_; }

  /**
   * @brief size
   * @return The combined size of the low rank matrices.
   */
  inline std::size_t size() const { return L_.size() + R_.size(); }

  template <typename Archive> void serialize(Archive &ar) {}

private:
  EigMat<T> L_;
  EigMat<T> R_;
  std::size_t rank_ = 0;
};

template <typename T>
inline std::ostream &operator<<(std::ostream &os, const LRTile<T> &t) {
  os << "{\n";
  os << "L \n" << t.matrixL();
  os << "\n\nR \n" << t.matrixR();
  os << "\n}";
  return os;
}

#endif // LR_TILE_H
