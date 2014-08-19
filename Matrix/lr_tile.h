#ifndef LR_TILE_H
#define LR_TILE_H

#include "../include/eigen.h"
#include <iostream> //TODO remove printing from file

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
  explicit LRTile(const EigMat<T> &input) : Q_(), R_() {
    Eigen::ColPivHouseholderQR<EigMat<T>> qr(input);
    qr.setThreshold(cut_);
    rank_ = qr.rank();

    Q_ = EigMat<T>(qr.householderQ()).leftCols(rank());
    R_ = qr.matrixR()
             .topLeftCorner(rank(), input.cols())
             .template triangularView<Eigen::Upper>();
    // TODO figure out what todo with the cols
    EigMat<T> P = qr.colsPermutation();
    // TODO remove printing from file.
    std::cout << "rank of input = " << rank() << "\nP\n" << P << std::endl;
  }

  /**
   * @brief rank returns the rank of the tile by value.
   */
  inline std::size_t rank() const { return rank_; }

private:
  EigMat<T> Q_;
  EigMat<T> R_;
  double cut_ = 1e-08;
  std::size_t rank_ = 0;
};

#endif // LR_TILE_H
