#ifndef LR_TILE_H
#define LR_TILE_H

#include "../include/eigen.h"
#include "../include/tiledarray.h"
#include "../include/lapack.h"
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

  EigMat<T> L = EigMat<T>(qr.householderQ()).leftCols(qr.rank());
  EigMat<T> R = EigMat<T>(qr.matrixR()
                              .topLeftCorner(qr.rank(), input.cols())
                              .template triangularView<Eigen::Upper>())
                * qr.colsPermutation().transpose();

  return std::make_pair(L, R);
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
  LRTile &operator=(const LRTile &rhs) {
    L_ = rhs.L_;
    R_ = rhs.R_;
    rank_ = rhs.rank_;
    return *this;
  }

  // LRTile(LRTile &&rhs) = default; Must write by hand :P.
  /**
   * @brief LRTile constructor builds a low rank representation of a matrix
   * given an input matrix.
   * @param input is an eigen matrix which with matching data type to the class.
   */
  explicit LRTile(const EigMat<T> &input, double cut = 1e-08)
      : L_(), R_(), range_() {
    auto QR_pair = detail::qr_decomp(input, cut);
    L_ = std::get<0>(QR_pair);
    R_ = std::get<1>(QR_pair);
    rank_ = L_.cols();
  }

  /**
   * @brief LRTile constructor builds a low rank representation of a matrix
   * given an input matrix.
   * @param input is an eigen matrix which with matching data type to the class.
   */
  explicit LRTile(TiledArray::Range range,
                  const EigMat<T> &input,
                  double cut = 1e-08)
      : L_(), R_(), range_(range) {
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
      : L_(left_mat), R_(right_mat), rank_(right_mat.rows()) {}

  /**
   * @brief mult multiples two tiles together.
   * @param a Low Rank tile to multiple the current tile with.
   * @return a new Low Rank tile.
   */
  LRTile mult(const LRTile<T> &right) const {
    assert(false); // TODO
    return LRTile();
  }

  /**
   * @brief add adds two low rank tiles togather
   * @param right the right side low rank being added to this
   * @return new low rank tile with rank this->rank() + right.rank()
   * @warning the output tile will have higher rank than the input tiles.
   */
  LRTile add(const LRTile<T> right) const {
    const auto new_rank = rank_ + right.rank_;
    EigMat<T> L(L_.rows(), new_rank);
    EigMat<T> R(new_rank, right.R_.cols());

    L.leftCols(rank()) = L_;
    L.rightCols(right.rank()) = right.L_;

    R.topRows(rank()) = R_;
    R.bottomRows(right.rank()) = right.R_;

    compress(L, R, 1e-08);

    return LRTile(L, R);
  }

  bool empty() const { return false; }

  LRTile add(const LRTile &right, const TiledArray::Permutation &perm) const {
    assert(false);
    return LRTile();
  }
  LRTile add(const LRTile &right,
             const numeric_type factor,
             const TiledArray::Permutation &perm) const {}
  LRTile
  add(const value_type &value, const TiledArray::Permutation &perm) const {}


  LRTile clone() const {
    assert(false); // TODO
    return LRTile();
  }
  LRTile permute(const TiledArray::Permutation &perm) const {
    assert(false); // TODO
    return LRTile();
  }
  LRTile scale(const numeric_type factor) const {
    assert(false); // TODO
    return LRTile();
  }
  LRTile
  scale(const numeric_type factor, const TiledArray::Permutation &perm) const {
    assert(false); // TODO
    return LRTile();
  }
  LRTile &scale_to(const numeric_type factor) {
    assert(false); // TODO
    return LRTile();
  }
  LRTile &add_to(const LRTile &right) {
    assert(false); // TODO
    return *this;
  }
  LRTile &add_to(const LRTile &right, const numeric_type factor) {
    assert(false); // TODO
    return *this;
  }
  LRTile &add_to(const value_type &value) {
    assert(false); // TODO
    return *this;
  }

  LRTile mult(const LRTile &right, const numeric_type factor) const {}
  LRTile mult(const LRTile &right, const TiledArray::Permutation &perm) const {
    assert(false); // TODO
    return LRTile();
  }
  LRTile mult(const LRTile &right,
              const numeric_type factor,
              const TiledArray::Permutation &perm) const {
    assert(false); // TODO
    return LRTile();
  }

  LRTile &mult_to(const LRTile &right) {
    assert(false); // TODO
    return *this;
  }
  LRTile &mult_to(const LRTile &right, const numeric_type factor) {
    assert(false); // TODO
    return *this;
  }

  /**
   * @brief compress attempts to reduce the rank of the tile
   * @param cut is the threshold parameter for Eigen::ColPivHouseholder.
   */
  void compress(EigMat<T> &L, EigMat<T> &R, double cut = 1e-8) const {
    const auto Lrows = L.rows();
    const auto Lcols = L.cols();

    double first_qr = madness::wall_time();
    Eigen::ColPivHouseholderQR<EigMat<T>> qr(L);
    first_qr = madness::wall_time() - first_qr;
    std::cout << "First Qr time = " << first_qr << "s\n";
    qr.setThreshold(cut);
    const auto rankL = qr.rank();

    double Form_Ll = madness::wall_time();
    L.resize(Lrows,rankL);
    L.setIdentity(Lrows, rankL);
    qr.householderQ().applyThisOnTheLeft(L);
    Form_Ll = madness::wall_time() - Form_Ll;
    std::cout << "Form ll time = " << Form_Ll << "s\n";

    double Form_Rl = madness::wall_time();
    EigMat<T> Rl = EigMat<T>(qr.matrixR()
                                 .topLeftCorner(rankL, Lcols)
                                 .template triangularView<Eigen::Upper>())
                   * qr.colsPermutation().transpose();
    Form_Rl = madness::wall_time() - Form_Rl;
    std::cout << "Form Rl time = " << Form_Rl << "s\n";


    double sec_qr = madness::wall_time();
    qr.compute(R);
    const auto rankR = qr.rank();
    sec_qr = madness::wall_time() - sec_qr;
    std::cout << "Second Qr time = " << sec_qr
              << "s\n\tTotal Qr time = " << sec_qr + first_qr << "s\n";

    double contract_middle = madness::wall_time();
    Rl = Rl * EigMat<T>(qr.householderQ()).leftCols(rankR);
    contract_middle = madness::wall_time() - contract_middle;
    std::cout << "Make middle time = " << contract_middle << "s\n";
    double finalize = 0;

    if (Rl.rows() >= Rl.cols()) {
      finalize = madness::wall_time();
      L *= Rl;
      R = EigMat<T>(qr.matrixR()
                        .topLeftCorner(qr.rank(), R.cols())
                        .template triangularView<Eigen::Upper>())
          * qr.colsPermutation().transpose();
      finalize = madness::wall_time() - finalize;
      std::cout << "finalize time = " << finalize << "s\n";
    } else {
      finalize = madness::wall_time();
      R = Rl * (qr.matrixR()
                             .topLeftCorner(qr.rank(), R.cols())
                             .template triangularView<Eigen::Upper>())
          * qr.colsPermutation().transpose();
      finalize = madness::wall_time() - finalize;
      std::cout << "finalize time = " << finalize << "s\n";
    }
    std::cout << "Total time = " << first_qr + sec_qr + contract_middle
                                    + finalize << "s\n";
  }


  LRTile operator*(const LRTile &right) {
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

  LRTile gemm(const LRTile &right,
              const LRTile::numeric_type factor,
              const TiledArray::math::GemmHelper &gemm_config) const {

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
    assert(L.cols() == R.rows());
    return LRTile(L, R);
  }

  // GEMM operation with fused indices as defined by gemm_config; multiply arg1
  // by arg2, return the result
  LRTile &gemm(const LRTile &left,
               const LRTile &right,
               const LRTile::numeric_type factor,
               const TiledArray::math::GemmHelper &gemm_config) {
    // TODO This is wasted copying so fix it.
    auto AB = left.gemm(right, factor, gemm_config);
    if (rank() != 0) {
      auto temp = AB.add(*this);
      *this = temp;
    } else {
      *this = AB;
    }
    return *this;
  }

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
  TiledArray::Range range_;
};

template <typename T>
inline std::ostream &operator<<(std::ostream &os, const LRTile<T> &t) {
  os << "{\n";
  os << "L^T\n" << t.matrixL().transpose();
  os << "\nR\n" << t.matrixR();
  os << "\n}\n";
  return os;
}

#endif // LR_TILE_H
