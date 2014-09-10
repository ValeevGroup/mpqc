#ifndef LR_TILE_H
#define LR_TILE_H

#include "../include/eigen.h"
#include "../include/tiledarray.h"
#include "tile_algebra.h"
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
LR_pair<T> qr_decomp(const EigMat<T> &input, bool &is_full_rank, double cut) {

    EigMat<T> L;
    EigMat<T> R;
    is_full_rank = algebra::ColPivQR(input, L, R, cut);

    // L_la and R_la will be empty if is_full_rank returns true
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
    typedef LRTile eval_type;             /// The data type of the tile
    typedef T value_type;                 /// The data type of the tile
    typedef TiledArray::Range range_type; /// Tensor range type
    typedef T numeric_type;               /// The scalar type of tile.
    typedef std::size_t size_type;

    template <typename U>
    using EigMat = detail::EigMat<U>;

    LRTile() = default;
    ~LRTile() = default;
    LRTile(const LRTile &rhs) = default;
    LRTile &operator=(const LRTile &rhs) = default;

    /**
     * @brief LRTile move constructor, must write by hand since eigen objects
     * don't have move constructors yet.
     * @param rhs a rvalue ref to a LRTile
     */
    LRTile(LRTile &&rhs) noexcept : L_(),
                                    R_(),
                                    rank_(std::move(rhs.rank_)),
                                    is_full_rank_(std::move(rhs.is_full_rank_)),
                                    range_(),
                                    cut_(std::move(rhs.cut_)) {
        range_.swap(rhs.range_);
        L_.swap(rhs.L_);
        R_.swap(rhs.R_);
    }

    /**
     * @brief operator = move assignment operator, see move constructor for
     * issue with eigen objects
     * @param rhs the LRTile being assigned from.
     * @return this LRTile
     */
    LRTile &operator=(LRTile &&rhs) noexcept {
        L_.swap(rhs.L_);
        R_.swap(rhs.R_);
        rank_ = std::move(rhs.rank_);
        is_full_rank_ = std::move(rhs.is_full_rank_);
        range_.swap(rhs.range_);
        cut_ = std::move(rhs.cut_);
        return *this;
    }

    /**
     * @brief LRTile this constructor takes a range object, and default
     * initializes everything else.
     * @param range a TiledArray::Range object.
     */
    explicit LRTile(TiledArray::Range range, double cut = 1e-07)
        : L_(), R_(), range_(std::move(range)), cut_(cut) {}

    /**
     * @brief LRTile constructor builds a low rank representation of a matrix
     * given an input matrix.
     * @param range is a TiledArray::Range object
     * @param input is an eigen matrix which with matching data type to the
     * class.
     * @param cut is the threshold used to determine the rank of the tile.
     */
    LRTile(TiledArray::Range range, const EigMat<T> &input, bool decomp = true,
           double cut = 1e-7)
        : L_(), R_(), range_(range), cut_(cut) {

        if (decomp) {
            auto QR_pair = detail::qr_decomp(input, is_full_rank_, cut_);
            if (!is_full_rank_) {
                L_ = std::get<0>(QR_pair);
                R_ = std::get<1>(QR_pair);
            } else { // Tile was full rank so just use original matrix.
                L_ = input;
            }
        } else {
            is_full_rank_ = true;
            L_ = input;
        }

        rank_ = (is_full_rank_) ? std::min(L_.cols(), L_.rows()) : L_.cols();
    }

    /**
     * @brief LRTile constructor builds a low rank representation of a
     * matrix
     * given an input matrix. If the matrix is not low rank then keep the
     * original input.
     * @param input is an eigen matrix which with matching data type to the
     * class.
     * @cut is the threshold passed to the decomposition for determining the
     * rank of the tile.
     */
    LRTile(const EigMat<T> &input, bool decomp = true, double cut = 1e-07)
        : LRTile(TiledArray::Range::Range(), input, decomp, cut) {}


    /**
     * @brief LRTile constructor which takes values for each member.
     * @param range a TiledArray::Range object.
     * @param left_mat an Eigen::Matrix for the left matrix.
     * @param right_mat an Eigen::Matrix for the right matrix.
     */
    LRTile(const TiledArray::Range &range, const EigMat<T> &left_mat,
           const EigMat<T> &right_mat, double cut = 1e-07)
        : L_(left_mat), R_(right_mat), rank_(left_mat.cols()), range_(range),
          cut_(cut) {}

    /**
     * @brief LRTile constructor to assign the low rank matrices directly
     * @param left_mat the left hand matrix
     * @param right_mat the right hand matrix
     */
    LRTile(const EigMat<T> &left_mat, const EigMat<T> &right_mat,
           double cut = 1e-07)
        : LRTile(TiledArray::Range::Range(), left_mat, right_mat, cut) {}

    /**
     * @brief LRTile move constructor to assign the low rank matrices
     * directly
     * @param range A TiledArray::Range object for the tile.
     * @param left_mat the left hand matrix
     * @param right_mat the right hand matrix
     */
    LRTile(TiledArray::Range &&range, EigMat<T> &&left_mat,
           EigMat<T> &&right_mat, double cut = 1e-07) noexcept
        : L_(),
          R_(),
          rank_(left_mat.cols()),
          range_(),
          cut_(std::move(cut)) {
        swap(range_, range);
        if (rank_ > double(std::min(left_mat.rows(), right_mat.cols())) / 2.0) {
            L_ = left_mat * right_mat;
            is_full_rank_ = true;
        } else {
            L_.swap(left_mat);
            R_.swap(right_mat);
        }
    }

    /**
     * @brief LRTile move constructor which takes a left and right matrix.
     * @param left_mat the left matrix being moved in.
     * @param right_mat the right matrix being moved in.
     */
    LRTile(EigMat<T> &&left_mat, EigMat<T> &&right_mat,
           double cut = 1e-07) noexcept
        : LRTile(TiledArray::Range::Range(), std::move(left_mat),
                 std::move(right_mat), cut) {}

    /**
     * @brief empty is the tile empty?
     */
    bool empty() const { return (rank_ == 0); }

    /**
     * @brief clone returns a copy of the current tile, user the copy ctor.
     */
    LRTile clone() const {
        return *this;
    } // TODO this will have to change when we switch to shallow copy

    /**
     * @brief compress attempts to reduce the rank of the tile
     * @param L,R are the matrice used to make the tile they will be
     * modified.
     * @param cut is the threshold parameter for Eigen::ColPivHouseholder.
     * @note L and R should be moved into compress, but Eigen does not yet
     * support moving.
     *
     * @Warning compress will modify L and R
     */
    LRTile compress(TiledArray::Range range, EigMat<T> &L, EigMat<T> &R,
                    double cut) const {
        // First L is taken by value second is modified.  Returns true if not
        // compressiable.
        if (!algebra::ColPivQR(L, L, R, cut)) {
            return LRTile(std::move(range), std::move(L), std::move(R), cut);
        } else {
            return LRTile(std::move(range), L * R, false, cut);
        }
    }

    void compress(EigMat<T> &L, EigMat<T> &R, double cut) {
        algebra::ColPivQR(L, L, R, cut);
    }

    /**
     * @brief operator * multiplies to tiles togather.
     * @param right the tile to multiply this by.
     */
    LRTile operator*(const LRTile &right) {
        // TODO determine if it is likely that two full_rank tiles can be
        // multiplied together into a low rank tile.
        // for now assume they cannot be.
        EigMat<T> L;
        EigMat<T> R;
        if (is_full() && right.is_full()) {
            return LRTile(range(), EigMat<T>(L_ * right.L_), false, cut_);
        } else if (right.is_full()) {
            L = matrixL();
            R = matrixR() * right.matrixL();
            return LRTile(std::move(L), std::move(R), cut_);
        } else if (is_full()) {
            L = matrixL() * right.matrixL();
            R = right.matrixR();
            return LRTile(std::move(L), std::move(R), cut_);
        }

        bool use_left_rank = (rank() < right.rank());
        L = matrixL();
        R = right.matrixR();
        const EigMat<T> mid = matrixR() * right.matrixL();

        (use_left_rank) ? R = mid *R : L *= mid;

        return LRTile(std::move(L), std::move(R), cut_);
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
     * @return the full sized matrix. If full rank just return L_, else
     * L_*R_
     * @warning Can possibly increase memory usage by a large amount.
     */
    inline EigMat<T> matrixLR() const {
        return (is_full()) ? L_ : algebra::cblas_gemm(L_, R_, 1.0);
    }

    /**
     * @brief rank returns the rank of the tile by value.
     */
    inline std::size_t rank() const { return rank_; }

    /**
     * @brief size
     * @return The combined size of the low rank matrices.
     */
    inline std::size_t size() const { return L_.size() + R_.size(); }

    inline bool is_full() const { return is_full_rank_; }

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
     * @param factor the factor which scales the array.
     * @warning scale assumes that the rank of the tile is not changed by
     * scaling, but this is not necessarily true.
     * @return a new LRTile which as been scaled by factor.
     */
    LRTile scale(const numeric_type factor) const {
        auto temp = clone();
        temp.L_ *= factor;
        return temp;
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
        } else if (rank() == 0) { // If this is empty copy right.
            return LRTile(right);
        } else if (is_full() && right.is_full()) { // If both full
            return LRTile(range(), EigMat<T>(L_ + right.L_), false, cut_);
        } else if (is_full() || right.is_full()) {
            return LRTile(range(), EigMat<T>(matrixLR() + right.matrixLR()),
                          false, cut_);
        }

        const auto new_rank = rank_ + right.rank_;

        if (new_rank > static_cast
            <unsigned long>(std::min(L_.rows(), R_.cols()))) {
            return LRTile(range(), EigMat<T>(matrixLR() + right.matrixLR()),
                          true, cut_);
        }

        EigMat<T> L(L_.rows(), new_rank);
        EigMat<T> R(new_rank, right.R_.cols());

        L << L_, right.L_;
        R << R_, right.R_;

        if (new_rank < double(std::min(L.rows(), R_.cols())) / 4.0) {
            return LRTile(range(), std::move(L), std::move(R), cut_);
        }

        return compress(range(), L, R, cut_);
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
        } else if (rank() == 0) {
            return LRTile(right);
        } else if (is_full() && right.is_full()) {
            return LRTile(range(), EigMat<T>(L_ - right.L_), false, cut_);
        } else if (is_full() || right.is_full()) {
            return LRTile(range(), EigMat<T>(matrixLR() - right.matrixLR()),
                          false, cut_);
        }

        const auto new_rank = rank_ + right.rank_;
        if (new_rank > static_cast
            <unsigned long>(std::min(L_.rows(), R_.cols()))) {
            return LRTile(range(), EigMat<T>(matrixLR() - right.matrixLR()),
                          true, cut_);
        }

        EigMat<T> L(L_.rows(), new_rank);
        EigMat<T> R(new_rank, right.R_.cols());

        L << L_, -right.L_; // Just added a minus here.
        R << R_, right.R_;

        if (new_rank < double(std::min(L.rows(), R_.cols())) / 4.0) {
            return LRTile(range(), std::move(L), std::move(R), cut_);
        }

        return compress(range(), L, R, cut_);
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
        assert(false);
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
     * @brief gemm perform a Gemm opertion on two tiles
     * @param right the left tile being multiplied by this one.
     * @param factor a scaling factor for the gemm.
     * @param gemm_config
     * @return
     */
    LRTile gemm(const LRTile &right, const LRTile::numeric_type factor,
                const TiledArray::math::GemmHelper &gemm_config) const {
        auto result_range = gemm_config.make_result_range
                            <range_type>(range(), right.range());
        EigMat<T> L;
        EigMat<T> R;
        // TODO decide whether is full or not later.
        if (is_full() && right.is_full()) {
            // return LRTile(std::move(result_range),
            // re              EigMat<T>(factor * L_ * right.L_), false);
            return LRTile(std::move(result_range),
                          algebra::cblas_gemm(EigMat<T>(factor * matrixL()),
                                              right.matrixL(), 1.0),
                          false, cut_);

        } else if (is_full()) {
            // L = matrixL() * right.matrixL();
            L = algebra::cblas_gemm(EigMat<T>(factor * matrixL()),
                                    right.matrixL(), 1.0);
            R = right.matrixR();
            return LRTile(std::move(result_range), std::move(L), std::move(R),
                          cut_);
        } else if (right.is_full()) {
            L = factor * matrixL();
            R = algebra::cblas_gemm(matrixR(), right.matrixL(), 1.0);
            return LRTile(std::move(result_range), std::move(L), std::move(R),
                          cut_);
        }

        bool use_left_rank = (rank() < right.rank());
        L = factor * matrixL();
        R = right.matrixR();


        const auto mid = algebra::cblas_gemm(matrixR(), right.matrixL(), 1.0);
        // const auto mid = matrixR() * right.matrixL();

        (use_left_rank) ? R = algebra::cblas_gemm(mid, R, 1.0) : L
            = algebra::cblas_gemm(L, mid, 1.0);
        //(use_left_rank) ? R = mid *R : L *= mid;

        return LRTile(std::move(result_range), std::move(L), std::move(R),
                      cut_);
    }

    void is_full_gemm(const LRTile &left, const LRTile &right,
                      const LRTile::numeric_type factor) {
        if (left.is_full()) {
            if (right.is_full()) {
                algebra::cblas_gemm_inplace(left.L_, right.L_, L_, factor);
            } else {
                EigMat<T> Lab = algebra::cblas_gemm(left.L_, right.L_, 1.0);
                algebra::cblas_gemm_inplace(Lab, right.R_, L_, factor);
            }
        } else if (right.is_full()) {
            EigMat<T> Rab = algebra::cblas_gemm(left.R_, right.L_, 1.0);
            algebra::cblas_gemm_inplace(left.L_, Rab, L_, factor);
        } else {
            EigMat<T> mid = algebra::cblas_gemm(left.R_, right.L_, 1.0);
            bool use_left_rank = (left.rank() < right.rank());
            if (use_left_rank) {
                EigMat<T> Rab = algebra::cblas_gemm(mid, right.R_, 1.0);
                algebra::cblas_gemm_inplace(left.L_, Rab, L_, factor);
            } else {
                EigMat<T> Lab = algebra::cblas_gemm(left.L_, mid, 1.0);
                algebra::cblas_gemm_inplace(Lab, right.R_, L_, factor);
            }
        }
    }

    void is_not_full_gemm(const LRTile &left, const LRTile &right,
                          const LRTile::numeric_type factor) {
        assert(!(left.is_full() && right.is_full()));

        auto mult_rank = std::min(left.rank(), right.rank());
        auto C_out_rank = mult_rank + rank();

        EigMat<T> L(L_.rows(), C_out_rank);
        EigMat<T> R(C_out_rank, R_.cols());

        if(mult_rank == 0){
          return;
        }

        L.rightCols(rank()) = L_;
        R.bottomRows(rank()) = R_;

        if (left.is_full()) {
            { // Write into the left side of L
                const int K = left.L_.cols();
                const int M = L.rows();
                const int N = mult_rank; // TODO fix things here.
                const int LDA = M, LDB = K, LDC = M;
                if (K == 0 || M == 0 || N == 0) {
                    L.leftCols(mult_rank) = factor * left.L_ * right.L_;
                } else {
                    madness::cblas::gemm(
                        madness::cblas::CBLAS_TRANSPOSE::NoTrans,
                        madness::cblas::CBLAS_TRANSPOSE::NoTrans, M, N, K,
                        factor, left.L_.data(), LDA, right.L_.data(), LDB, 0,
                        L.data(), LDC);
                }
            }
            { // Write into the top of R
                R.topRows(mult_rank) = right.R_;
            }
        } else if (right.is_full()) {
            { // Write into the top of R
                const int K = right.L_.rows();
                const int M = mult_rank;
                const int N = R.cols();
                const int LDA = M, LDB = K, LDC = R.rows();
                if (K == 0 || M == 0 || N == 0) {
                    R.topRows(mult_rank) = factor * left.R_ * right.L_;
                } else {
                    madness::cblas::gemm(
                        madness::cblas::CBLAS_TRANSPOSE::NoTrans,
                        madness::cblas::CBLAS_TRANSPOSE::NoTrans, M, N, K,
                        factor, left.R_.data(), LDA, right.L_.data(), LDB, 0,
                        R.data(), LDC);
                }
            }
            { // Write into the top of R
                L.leftCols(mult_rank) = left.L_;
            }

        } else {
            EigMat<T> mid = algebra::cblas_gemm(left.R_, right.L_, 1.0);
            bool use_left_rank = (left.rank() >= right.rank());
            if (use_left_rank) {
                // L.leftCols(mult_rank) = left.L_ * mid;
                { // Write into the left side of L
                    const int K = left.L_.cols();
                    const int M = L.rows();
                    const int N = mult_rank;
                    const int LDA = M, LDB = K, LDC = M;
                    if (K == 0 || M == 0 || N == 0) {
                        L.leftCols(mult_rank) = decltype(L)::Zero(M,N);
                    } else {
                        madness::cblas::gemm(
                            madness::cblas::CBLAS_TRANSPOSE::NoTrans,
                            madness::cblas::CBLAS_TRANSPOSE::NoTrans, M, N, K,
                            factor, left.L_.data(), LDA, mid.data(), LDB, 0,
                            L.data(), LDC);
                    }
                }
                { // Write into the top of R
                    R.topRows(mult_rank) = right.R_;
                }
            } else {
                { // Write into the top of R
                    //    R.topRows(mult_rank) = mid * right.R_;
                    const int K = right.R_.rows();
                    const int M = mult_rank;
                    const int N = R.cols();
                    const int LDA = M, LDB = K, LDC = R.rows();
                    if (K == 0 || M == 0 || N == 0) {
                        R.topRows(mult_rank) = factor * mid * right.R_;
                    } else {
                        madness::cblas::gemm(
                            madness::cblas::CBLAS_TRANSPOSE::NoTrans,
                            madness::cblas::CBLAS_TRANSPOSE::NoTrans, M, N, K,
                            factor, mid.data(), LDA, right.R_.data(), LDB, 0,
                            R.data(), LDC);
                    }
                }
                { // Write into the top of R
                    L.leftCols(mult_rank) = left.L_;
                }
            }
        }
        // compress(L,R,cut_);
        L_.swap(L);
        R_.swap(R);
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

        range_ = gemm_config.make_result_range
                 <range_type>(range(), right.range());

        if (is_full()) {
            is_full_gemm(left, right, factor);
        } else {
            auto mult_rank = std::min(left.rank(), right.rank());
            auto C_rank = mult_rank + rank();

            if ((C_rank >= 0.5 * std::min(L_.rows(), R_.cols()))
                || (left.is_full() && right.is_full())) {
                L_ = matrixLR();
                R_.resize(0, 0);
                rank_ = std::min(L_.rows(), L_.cols());
                is_full_rank_ = true;
                is_full_gemm(left, right, factor);
            } else {
                is_not_full_gemm(left, right, factor);
            }
        }

        rank_ = L_.cols();

        return *this;
    }


    template <typename Archive>
    void serialize(Archive &ar) {}

  private:
    EigMat<T> L_;
    EigMat<T> R_;
    std::size_t rank_ = 0;
    bool is_full_rank_ = false;
    TiledArray::Range range_;
    double cut_;
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
