#ifndef LR_TILE_H
#define LR_TILE_H

#include "../include/eigen.h"

/**
 * LRTile is a class which can be used as a low rank tile for TiledArray
 * matrices
 */
template <typename T> class LRTile {
public:
  typedef T value_type; /// The data type of the tile, typically double.

  template <typename U>
  using EigMat = Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>;


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
    Eigen::ColPivHouseholderQR<EigMat<T>> qr(input);
    qr.setThreshold(cut);
    rank_ = qr.rank();

    L_ = EigMat<T>(qr.householderQ()).leftCols(rank());
    R_ = qr.matrixR()
             .topLeftCorner(rank(), input.cols())
             .template triangularView<Eigen::Upper>();
    R_ *= qr.colsPermutation().transpose();
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
  LRTile mult(const LRTile<T> right) const {
    // Check which rank is smaller
    bool use_left_rank = (rank() < right.rank());
    EigMat<T> L;
    EigMat<T> R;
    if(use_left_rank){
      L = matrixL();
      R = (matrixR() * right.matrixL()) * right.matrixR();
    } else{
      L = matrixL() * (matrixR() * right.matrixL());
      R = right.matrixR();
    }
    return LRTile(L,R);
  }

  /**
   * @brief matrixL
   * @return the left hand side low rank matrix.
   */
  inline const EigMat<T> &matrixL() const {return L_;}

  /**
   * @brief matrixR
   * @return the right hand side low rank matrix.
   */
  inline const EigMat<T> &matrixR() const {return R_;}

  /**
   * @brief matrixLR
   * @return the full sized product of the low rank matrices.
   * @warning Can possibly increase memory use by a large amount.
   */
  inline EigMat<T> matrixLR() const {return L_ * R_;}

  /**
   * @brief rank returns the rank of the tile by value.
   */
  inline std::size_t rank() const { return rank_; }

private:
  EigMat<T> L_;
  EigMat<T> R_;
  std::size_t rank_ = 0;
};

#endif // LR_TILE_H
