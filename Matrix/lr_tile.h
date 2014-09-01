#ifndef LR_TILE_H
#define LR_TILE_H

#include "../include/eigen.h"
#include "../include/tiledarray.h"
#include "lin_algebra.h"
#include <utility>
#include <iostream>

namespace detail {

template <typename T>
using EigMat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T>
using LR_pair = std::pair<EigMat<T>, EigMat<T>>;

/**
* qr_decomp returns the low rank Q and R matrices as a pair.
*/
template <typename T>
LR_pair<T> qr_decomp(const EigMat<T> &input, double cut) {
    double lapack_time = madness::wall_time();
    EigMat<T> L_la;
    EigMat<T> R_la;
    ColPivQR(input, L_la, R_la, cut);
    lapack_time = madness::wall_time() - lapack_time;

    double eigen_time = madness::wall_time();
    Eigen::ColPivHouseholderQR<EigMat<T>> qr(input);
    qr.setThreshold(cut);

    EigMat<T> L = EigMat<T>(qr.householderQ()).leftCols(qr.rank());
    EigMat<T> R = EigMat<T>(qr.matrixR()
                                .topLeftCorner(qr.rank(), input.cols())
                                .template triangularView<Eigen::Upper>())
                  * qr.colsPermutation().transpose();
    eigen_time = madness::wall_time() - eigen_time;

    // std::cout << "Eigen takes " << eigen_time << "s\n";
    // ststd::cout << "Lapack takes " << lapack_time << "s\n\n";

    return std::make_pair(L, R);
}

} // namespace detail

/**
 * LRTile is a class which can be used as a low rank tile for TiledArray
 * matrices
 */
template <typename T>
class LRTile {
  public:
    typedef LRTile eval_type; /// The data type of the tile
    typedef T value_type;     /// The data type of the tile
    typedef TiledArray::Range range_type; /// Tensor range type
    typedef T numeric_type; /// The scalar type of tile.
    typedef std::size_t size_type;

    template <typename U>
    using EigMat = detail::EigMat<U>;

    LRTile() = default;
    ~LRTile() = default;
    LRTile(const LRTile &rhs) = default;

    LRTile(TiledArray::Range range)
        : L_(), R_(), rank_(), range_(std::move(range)) {}

    LRTile(LRTile &&rhs) noexcept : L_(),
                                    R_(),
                                    rank_(std::move(rhs.rank_)),
                                    range_(std::move(rhs.range_)) {
        L_.swap(rhs.L_);
        R_.swap(rhs.R_);
    }

    LRTile &operator=(const LRTile &rhs) {
        L_ = rhs.L_;
        R_ = rhs.R_;
        rank_ = rhs.rank_;
        return *this;
    }

    LRTile &operator=(LRTile &&rhs) noexcept {
        L_.swap(rhs.L_);
        R_.swap(rhs.R_);
        rank_ = std::move(rhs.rank_);
        range_.swap(rhs.range_);
        return *this;
    }

    /**
     * @brief LRTile constructor builds a low rank representation of a matrix
     * given an input matrix.
     * @param input is an eigen matrix which with matching data type to the
     * class.
     */
    explicit LRTile(const EigMat<T> &input, double cut = 1e-16)
        : L_(), R_(), range_() {
        auto QR_pair = detail::qr_decomp(input, cut);
        L_ = std::get<0>(QR_pair);
        R_ = std::get<1>(QR_pair);
        rank_ = L_.cols();
    }

    /**
     * @brief LRTile constructor builds a low rank representation of a matrix
     * given an input matrix.
     * @param input is an eigen matrix which with matching data type to the
     * class.
     */
    explicit LRTile(TiledArray::Range range, const EigMat<T> &input,
                    double cut = 1e-16)
        : L_(), R_(), range_(range) {
        auto QR_pair = detail::qr_decomp(input, cut);
        L_ = std::get<0>(QR_pair);
        R_ = std::get<1>(QR_pair);
        rank_ = L_.cols();
    }

    /**
     * @brief LRTile constructor to assign the low rank matrices directly
     * @param left_mat the left hand matrix
     * @param right_mat the right hand matrix
     */
    explicit LRTile(const EigMat<T> &left_mat, const EigMat<T> &right_mat)
        : L_(left_mat), R_(right_mat), rank_(right_mat.rows()), range_() {}

    /**
     * @brief LRTile constructor which takes values for each member.
     * @param range a TiledArray::Range object.
     * @param left_mat an Eigen::Matrix for the left matrix.
     * @param right_mat an Eigen::Matrix for the right matrix.
     */
    explicit LRTile(const TiledArray::Range &range, const EigMat<T> &left_mat,
                    const EigMat<T> &right_mat)
        : L_(left_mat), R_(right_mat), rank_(right_mat.rows()), range_(range) {}

    /**
     * @brief LRTile move constructor to assign the low rank matrices directly
     * @param range A TiledArray::Range object for the tile.
     * @param left_mat the left hand matrix
     * @param right_mat the right hand matrix
     */
    explicit LRTile(TiledArray::Range &&range, EigMat<T> &&left_mat,
                    EigMat<T> &&right_mat) noexcept : L_(),
                                                      R_(),
                                                      rank_(right_mat.rows()),
                                                      range_() {
        swap(range_, range);
        L_.swap(left_mat);
        R_.swap(right_mat);
    }

    /**
     * @brief empty is the tile empty?
     */
    bool empty() const { return (rank_ == 0); }

    /**
     * @brief clone returns a copy of the current tile.
     */
    LRTile clone() const { return LRTile(range(), matrixL(), matrixR()); }

    /**
     * @brief compress attempts to reduce the rank of the tile
     * @param L,R are the matrice used to make the tile they will be modified.
     * @param cut is the threshold parameter for Eigen::ColPivHouseholder.
     * @note L and R should be moved into compress, but Eigen does not yet
     * support
     * moving.
     */
    LRTile compress(TiledArray::Range range, EigMat<T> &L, EigMat<T> &R,
                    double cut = 1e-16) const {
        const auto Lrows = L.rows();
        const auto Lcols = L.cols();
        const auto Rcols = R.cols();

        Eigen::ColPivHouseholderQR<EigMat<T>> qr(L);
        qr.setThreshold(cut);
        const auto rankL = qr.rank();

        // This seems to be the most efficent way to form the Q matrix,
        // Eigen's documentation is not very clear on the subject.
        L.resize(Lrows, rankL);
        L.setIdentity(Lrows, rankL);
        qr.householderQ().applyThisOnTheLeft(L);

        EigMat<T> Rl = EigMat<T>(qr.matrixR()
                                     .topLeftCorner(rankL, Lcols)
                                     .template triangularView<Eigen::Upper>())
                       * qr.colsPermutation().transpose();

        qr.compute(R);
        const auto rankR = qr.rank();

        Rl = Rl * EigMat<T>(qr.householderQ()).leftCols(rankR);

        R = EigMat
            <T>(qr.matrixR().topLeftCorner(rankR, Rcols).template triangularView
                <Eigen::Upper>()) * qr.colsPermutation().transpose();

        // Contract into smaller index.
        (rankL > rankR) ? L *= Rl : R = Rl * R;

        return LRTile(range, std::move(L), std::move(R));
    }

    /**
     * @brief operator * multiplies to tiles togather.
     * @param right the tile to multiply this by.
     */
    LRTile operator*(const LRTile &right) {
        // Check which rank is smaller
        bool use_left_rank = (rank() < right.rank());

        EigMat<T> L = matrixL();
        EigMat<T> R = right.matrixR();
        const auto mid = matrixR() * right.matrixL();

        (use_left_rank) ? R = mid *R : L *= mid;

        return LRTile(std::move(L), std::move(R));
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

    /**
     * @brief range returns the TiledArray::Range object associated with the
     * tile.
     */
    inline TiledArray::Range range() const { return range_; }

    /********************************************************
     * BEGIN SECTION FOR TILEDARRAY MATH FUNCTIONS
     ********************************************************/

    LRTile permute(const TiledArray::Permutation &perm) const {
        assert(false); // TODO
        return LRTile();
    }

    /**
     * @brief scale
     * @param factor
     * @return
     */
    LRTile scale(const numeric_type factor) const {
        return LRTile(std::move(range()), std::move(factor * matrixL()),
                      std::move(matrixR()));
    }

    LRTile scale(const numeric_type factor,
                 const TiledArray::Permutation &perm) const {
        assert(false); // TODO
        return LRTile();
    }

    /**
     * @brief scale_to scales the current tile by a constant factor.
     * @param factor a scalor type to scale the tile by.
     * @return this.
     */
    LRTile &scale_to(const numeric_type factor) {
        L_ = factor * L_;
        return *this;
    }


    /**
     * @brief add adds two low rank tiles togather
     * @param right the right side low rank being added to this
     * @return new low rank tile with rank this->rank() + right.rank()
     * @warning the output tile will have higher rank than the input tiles.
     */
    LRTile add(const LRTile<T> &right) const {
        // If right is empty just copy this.
        if (right.rank_ == 0) {
            return LRTile(*this);
        }

        const auto new_rank = rank_ + right.rank_;

        EigMat<T> L(L_.rows(), new_rank);
        EigMat<T> R(new_rank, right.R_.cols());

        L << L_, right.L_;
        R << R_, right.R_;

        return compress(range(), L, R, 1e-08);
    }

    LRTile add(const LRTile &right, const TiledArray::Permutation &perm) const {
        assert(false);
        return LRTile();
    }

    LRTile add(const LRTile &right, const numeric_type factor,
               const TiledArray::Permutation &perm) const {
        assert(false);
        return LRTile();
    }

    LRTile
    add(const value_type &value, const TiledArray::Permutation &perm) const {
        assert(false);
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

    /**
     * @brief subt
     * @param right tile to subtract from this
     */
    LRTile subt(const LRTile &right) const {
        // If right is empty just copy this.
        if (right.rank_ == 0) {
            return LRTile(*this);
        }

        const auto new_rank = rank_ + right.rank_;

        EigMat<T> L(L_.rows(), new_rank);
        EigMat<T> R(new_rank, right.R_.cols());

        L << L_, right.L_; // Just added a minus here.
        R << R_, -right.R_;

        return compress(range(), L, R, 1e-16);
    }

    /**
     * @brief subt
     * @param right
     * @param perm
     * @return
     */
    LRTile
    subt(const LRTile &right, const TiledArray::Permutation &perm) const {
        // TODO FIX
        return subt(right);
    }

    LRTile &neg_to() {
        assert(false);
        return *this;
    }

    /**
     * @brief subt_to
     * @param right
     * @return
     */
    LRTile &subt_to(const LRTile &right) {
        *this = subt(right);
        return *this;
    }

    /**
     * @brief subt_to subtracts a tile from the current tile.
     * @param right the tile to subtract from this
     * @param factor a scalor which scales the result of the subtraction
     * @return this.
     */
    LRTile &subt_to(const LRTile &right, numeric_type factor) {
        // assert(false);
        *this = subt(right).scale_to(factor);
        return *this;
    }

    LRTile neg() const {
        assert(false);
        return LRTile();
    }

    LRTile neg(const TiledArray::Permutation &perm) const {
        assert(false);
        return LRTile();
    }


    LRTile mult(const LRTile<T> &right) const {
        assert(false); // TODO
        return LRTile();
    }

    LRTile mult(const LRTile &right, const numeric_type factor) const {
        assert(false); // TODO
        return LRTile();
    }

    LRTile
    mult(const LRTile &right, const TiledArray::Permutation &perm) const {
        assert(false); // TODO
        return LRTile();
    }

    LRTile mult(const LRTile &right, const numeric_type factor,
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
     * @brief gemm
     * @param right
     * @param factor
     * @param gemm_config
     * @return
     */
    LRTile gemm(const LRTile &right, const LRTile::numeric_type factor,
                const TiledArray::math::GemmHelper &gemm_config) const {

        bool use_left_rank = (rank() < right.rank());

        EigMat<T> L = factor * matrixL();
        EigMat<T> R = right.matrixR();

        // auto mid = cblas_gemm(matrixR(), right.matrixL());
        const auto mid = matrixR() * right.matrixL();

        //(use_left_rank) ? R = cblas_gemm(mid, R) : L = cblas_gemm(L, mid);
        (use_left_rank) ? R = mid *R : L *= mid;

        return LRTile(range(), std::move(L), std::move(R));
    }

    /**
     * @brief gemm
     * @param left
     * @param right
     * @param factor
     * @param gemm_config
     * @return
     */
    LRTile &gemm(const LRTile &left, const LRTile &right,
                 const LRTile::numeric_type factor,
                 const TiledArray::math::GemmHelper &gemm_config) {
        *this = left.gemm(right, factor, gemm_config).add(*this);
        return *this;
    }


    template <typename Archive>
    void serialize(Archive &ar) {}

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
