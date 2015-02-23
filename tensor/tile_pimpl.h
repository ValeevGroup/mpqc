#pragma once
#ifndef TCC_MATRIX_TILE_PIMPLE_H
#define TCC_MATRIX_TILE_PIMPLE_H

#include "tile_variant.h"
#include "tile_ops.h"
#include "tile_mutations.h"
#include "../include/tiledarray.h"

#include <memory>

namespace tcc {
namespace tensor {

template <typename T>
class TilePimpl {
  public:
    typedef TilePimpl eval_type;
    typedef T value_type;
    typedef TiledArray::Range range_type;
    typedef T numeric_type;
    typedef std::size_t size_type;


  public:
    TilePimpl() = default;
    ~TilePimpl() = default;
    TilePimpl(TilePimpl const &t) = default;
    TilePimpl &operator=(TilePimpl const &t) = default;
    TilePimpl(TilePimpl &&t) = default;
    TilePimpl &operator=(TilePimpl &t) = default;

    /*
     * User defined constructors
     */

    // Just range ctor
    TilePimpl(TiledArray::Range r) : tile_(), range_(std::move(r)) {}


    TilePimpl(TiledArray::Range r, double cut)
        : tile_(), range_(std::move(r)), cut_(cut) {}

    // TileVariant ctors
    TilePimpl(TiledArray::Range r, TileVariant<T> t)
        : tile_(std::make_shared<TileVariant<T>>(std::move(t))),
          range_(std::move(r)) {}

    explicit TilePimpl(TiledArray::Range r, TileVariant<T> t, double cut)
        : tile_(std::make_shared<TileVariant<T>>(std::move(t))),
          range_(std::move(r)),
          cut_(cut) {}

    // Clone will make a deep copy
    TilePimpl clone() const { return TilePimpl(range_, *tile_, cut_); }

    /*
     * Tile information functions
     */
    bool isFull() const { return tile_->tag(); }
    unsigned long rank() const { return tile_->rank(); }
    TiledArray::Range const &range() const { return range_; }
    double cut() const { return cut_; }
    bool empty() const { return !tile_; }
    double norm() const { return tile_->norm(); }

    // maybe expensive
    void setCut(double cut) {
        const auto temp = cut_;
        cut_ = cut;
        if (temp <= cut_) {
            // TODO_TCC recompress may save some space.
        }
    }

    TileVariant<T> const &tile() const { return *tile_; }
    TileVariant<T> &tile() { return *tile_; }

  private:
    std::shared_ptr<TileVariant<T>> tile_;
    TiledArray::Range range_;
    double cut_ = 1e-7;


  public:
    /*
     * Tile operations
     */
    TilePimpl gemm(TilePimpl const &right, numeric_type factor,
                   TiledArray::math::GemmHelper const &gemm_config) const {

        auto result_range
            = gemm_config.make_result_range<range_type>(range(), right.range());

        return TilePimpl(std::move(result_range),
                         tile_->apply_binary_op(
                             right.tile(), binary_ops::gemm_functor(factor)),
                         std::max(cut(), right.cut()));
    }

    TilePimpl &gemm(TilePimpl const &left, TilePimpl const &right,
                    numeric_type factor,
                    TiledArray::math::GemmHelper const &gemm_config) {
        range_
            = gemm_config.make_result_range<range_type>(range(), right.range());

        // Will convert to full when gemm grows the rank to much.
        if (tile_->tag() == TileVariant<T>::LowRank
            && (left.tile().tag() == TileVariant<T>::LowRank
                || right.tile().tag() == TileVariant<T>::LowRank)) {
            convert_to_full(*tile_, left.tile(), right.tile());
        }

        tile_->apply_ternary_mutation(left.tile(), right.tile(),
                                      ternary_mutations::gemm_functor(factor));

        return *this;
    }

    TilePimpl &add_to(TilePimpl const &right) {
        tile_->apply_binary_mutation(right.tile(),
                                     binary_mutations::subt_functor(-1.0));
        return *this;
    }

    TilePimpl add(TilePimpl const &right) const {
        return TilePimpl{range_, tile_->apply_binary_op(
                                     right.tile(), binary_ops::add_functor(1)),
                         cut_};
    }

    TilePimpl
    add(TilePimpl const &right, TiledArray::Permutation const &perm) const {
        assert(false);
        return TilePimpl();
    }

    TilePimpl
    mult(TilePimpl const &right, TiledArray::Permutation const &perm) const {
        assert(false);
        return TilePimpl();
    }

    TilePimpl mult(TilePimpl const &right) const {
        assert(false);
        return TilePimpl();
    }

    TilePimpl &mult_to(TilePimpl const &right) {
        assert(false);
        return *this;
    }


    TilePimpl permute(TiledArray::Permutation const &perm) const {
        const auto dim = perm.dim();
        if (dim != 3) {
            assert(false);
        }

        // Get double ref so we have auto complete
        TileVariant<double> const &tile = *tile_;

        class Permute_ {
          private:
            TiledArray::Permutation p_;
            std::array<int, 3> sizes_;

          public:
            Permute_(TiledArray::Permutation const &p, int X, int i, int j)
                : p_(p), sizes_{{X, i, j}} {}

            TileVariant<T> operator()(FullRankTile<double> const &f) const {
                auto size_ij = sizes_[1] * sizes_[2];

                Eigen::MatrixXd out_mat(sizes_[0], size_ij);
                auto const &in_mat = f.matrix();

                std::cout << "In shape = " << in_mat.rows() << " "
                          << in_mat.cols() << std::endl;
                std::cout << "Out shape = " << out_mat.rows() << " "
                          << out_mat.cols() << std::endl;

                for (auto X = 0ul; X < sizes_[0]; ++X) {
                    for (auto i = 0ul; i < sizes_[1]; ++i) {
                        for (auto j = 0ul; j < sizes_[2]; ++j) {
                            out_mat(X, j * sizes_[2] + i)
                                = in_mat(X, i * sizes_[1] + j);
                        }
                    }
                }

                return TileVariant<T>{FullRankTile<double>{std::move(out_mat)}};
            }

            TileVariant<T> operator()(LowRankTile<double> const &lr) const {
                auto size_ij = sizes_[1] * sizes_[2];

                auto const &in_r = lr.matrixR();
                Eigen::MatrixXd out_r(lr.rank(), size_ij);

                for (auto X = 0ul; X < lr.rank(); ++X) {
                    for (auto i = 0ul; i < sizes_[1]; ++i) {
                        for (auto j = 0ul; j < sizes_[2]; ++j) {
                            out_r(X, j * sizes_[2] + i)
                                = in_r(X, i * sizes_[1] + j);
                        }
                    }
                }

                return TileVariant<T>{
                    LowRankTile<double>{lr.matrixL(), std::move(out_r)}};
            }
        };

        auto const &sizes = range_.size();
        return TilePimpl{range_, tile.apply_unary_op(Permute_(
                                     perm, sizes[0], sizes[1], sizes[2])),
                         cut()};
    }

    void compress() {
        if (!tile_->iszero()) {
            tile_->apply_unary_mutation(
                unary_mutations::compress_functor(cut()));
        }
    }

    void set_zero() { *tile_ = TileVariant<T>{LowRankTile<T>{true}}; }

    TilePimpl scale(const T factor, TiledArray::Permutation const &perm) const {
        assert(false);
    }

    TilePimpl scale(const T factor) const {
        if (!tile_->iszero()) {
            return TilePimpl{
                range(),
                tile_->apply_unary_op(unary_ops::scale_functor(factor)), cut()};
        } else {
            return TilePimpl{range(), TileVariant<T>{LowRankTile<T>{true}},
                             cut()};
        }
    }

    TilePimpl &scale_to(const T factor) {
        if (!tile_->iszero()) {
            tile_->apply_unary_mutation(unary_mutations::scale_functor(factor));
        }
        return *this;
    }

    TilePimpl &subt_to(TilePimpl const &right) {
        if (tile_->iszero()) { // We are zero
            if (!right.tile().iszero()) {
                *tile_ = right.tile().apply_unary_op(
                    unary_ops::scale_functor(-1.0));
            }
        } else if (!right.tile().iszero()) { // both aren't zero
            tile_->apply_binary_mutation(right.tile(),
                                         binary_mutations::subt_functor(1.0));
        }
        return *this;
    }

    // TODO_TCC figure out what is going on with needing factor second.
    TilePimpl &subt_to(TilePimpl const &right, T factor) {
        if (tile_->iszero()) {
            if (!right.tile().iszero()) {
                *tile_ = right.tile().apply_unary_op(
                    unary_ops::scale_functor(-1.0));
            }
        } else if (!right.tile().iszero()) {
            convert_to_full(*tile_, right.tile());
            tile_->apply_binary_mutation(right.tile(),
                                         binary_mutations::subt_functor(1.0));
            tile_->apply_unary_mutation(unary_mutations::scale_functor(factor));
        }
        return *this;
    }


    TilePimpl
    subt(TilePimpl const &right, TiledArray::Permutation const &perm) const {
        assert(false);
    }

    TilePimpl subt(TilePimpl const &right) const {
        if (tile_->iszero()) {
            if (right.tile().iszero()) {
                return clone();
            } else {
                return TilePimpl{range(), right.tile().apply_unary_op(
                                              unary_ops::scale_functor(-1.0)),
                                 cut()};
            }
        } else if (right.tile().iszero()) {
            return clone();
        } else {
            return TilePimpl{range(),
                             tile_->apply_binary_op(
                                 right.tile(), binary_ops::subt_functor(1.0)),
                             std::max(cut(), right.cut())};
        }
    }

    TilePimpl neg(TiledArray::Permutation const &perm) const { assert(false); }
    TilePimpl neg() const { assert(false); }

    TilePimpl &neg_to() { assert(false); }


    template <typename Archive>
    typename madness::enable_if<madness::archive::is_output_archive<Archive>>::
        type
        serialize(Archive &ar) {
        int tag = (tile_) ? tile_->tag() : -1;
        ar &tag &cut_ &range_;
        if (tag != -1) {
            if (tag == 0) {
                auto rows = static_cast<unsigned long>(
                    tile_->lrtile().matrixL().rows());
                auto cols = static_cast<unsigned long>(
                    tile_->lrtile().matrixL().cols());
                ar &rows &cols;
                rows = tile_->lrtile().matrixR().rows();
                cols = tile_->lrtile().matrixR().cols();
                ar &rows &cols;
                ar &madness::archive::wrap(tile_->lrtile().matrixL().data(),
                                           tile_->lrtile().matrixL().size());
                ar &madness::archive::wrap(tile_->lrtile().matrixR().data(),
                                           tile_->lrtile().matrixR().size());

            } else {
                auto rows = static_cast<unsigned long>(
                    tile_->ftile().matrix().rows());
                auto cols = static_cast<unsigned long>(
                    tile_->ftile().matrix().cols());
                ar &rows &cols;
                ar &madness::archive::wrap(tile_->ftile().matrix().data(),
                                           tile_->ftile().matrix().size());
            }
        }
    }

    template <typename Archive>
    typename madness::enable_if<madness::archive::is_input_archive<Archive>>::
        type
        serialize(Archive &ar) {
        int tag = 0;
        ar &tag &cut_ &range_;
        if (tag != -1) {
            if (tag == 0) {
                auto lrows = 0ul, lcols = 0ul, rrows = 0ul, rcols = 0ul;
                ar &lrows &lcols &rrows &rcols;
                typename LowRankTile<T>::template Matrix<T> L(lrows, lcols);
                typename LowRankTile<T>::template Matrix<T> R(rrows, rcols);
                ar &madness::archive::wrap(L.data(), L.size())
                    & madness::archive::wrap(R.data(), R.size());
                LowRankTile<double> l{std::move(L), std::move(R)};
                tile_ = std::make_shared<TileVariant<double>>(std::move(l));
            } else {
                auto rows = 0ul, cols = 0ul;
                ar &rows &cols;
                typename FullRankTile<T>::template Matrix<T> mat(rows, cols);
                ar &madness::archive::wrap(mat.data(), mat.size());
                FullRankTile<double> f{std::move(mat)};
                tile_ = std::make_shared<TileVariant<double>>(std::move(f));
            }
        }
    }


  private:
    /*
     *  Utility Functions
     */
    template <typename First, typename Second, typename... Rest>
    unsigned long
    mutation_rank(First first, Second second, Rest... rest) const {
        return first + std::min({second, rest...});
    }

    // May convert the tile to a full rank tile, but doesn't have to.
    template <typename... Rest>
    void convert_to_full(TileVariant<T> &result, Rest... rest) const {
        auto get_rank = [&](TileVariant<T> const &t) { return t.rank(); };
        const auto out_rank = mutation_rank(result.rank(), get_rank(rest)...);

        if (double(out_rank) > 0.5 * result.full_rank()) {
            result = TileVariant<T>{result.matrix()};
        }
    }
};

} // namespace tensor
} // namespace tcc

#endif // TCC_MATRIX_TILE_PIMPLE_H
