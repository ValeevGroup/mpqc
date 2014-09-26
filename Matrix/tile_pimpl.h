#ifndef TCC_MATRIX_TILE_PIMPLE_H
#define TCC_MATRIX_TILE_PIMPLE_H

#include "tile_variant.h"
#include "tile_ops.h"
#include "tile_mutations.h"
#include "../include/tiledarray.h"
#include <memory>

template <typename T>
class TilePimpl {
  public:
    /*
     * TiledArray typedefs
     */
    using eval_type = TilePimpl;
    using value_type = T;
    using range_type = TiledArray::Range;
    using numeric_type = T;
    using size_type = std::size_t;

  public:
    TilePimpl() = default;
    ~TilePimpl() = default;
    TilePimpl(TilePimpl const &t) = default;
    TilePimpl &operator=(TilePimpl const &t) = default;
    TilePimpl(TilePimpl &&t) = default;
    TilePimpl &operator=(TilePimpl &&t) = default;

    /*
     * User defined constructors
     */

    // Just range ctor
    explicit TilePimpl(TiledArray::Range r, double cut)
        : tile_(), range_(std::move(r)), cut_{cut} {}
    explicit TilePimpl(TiledArray::Range r) : TilePimpl{std::move(r), 1e-07} {}

    // TileVariant ctors
    explicit TilePimpl(TiledArray::Range r, TileVariant<T> t, double cut)
        : tile_(std::make_shared<TileVariant<T>>(std::move(t))),
          range_(std::move(r)), cut_{cut} {}
    explicit TilePimpl(TiledArray::Range r, TileVariant<T> t)
        : TilePimpl{std::move(r), std::move(t), 1e-07} {}

    // Clone will make a deep copy
    TilePimpl clone() const { return TilePimpl{range_, *tile_, cut_}; }

    /*
     * Tile information functions
     */
    bool isFull() const { return tile_->tag(); }
    unsigned long rank() const { return tile_->rank(); }
    TiledArray::Range const &range() const { return range_; }
    double cut() const { return cut_; }
    bool empty() const { return !tile_; }

    TileVariant<T> const &tile() const { return *tile_; }

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
        assert(false);
        return *this;
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
        assert(false);
        return TilePimpl();
    }

    void compress() {
        tile_->apply_unary_mutation(unary_mutations::compress_functor(cut()));
    }

    TilePimpl scale(const T factor, TiledArray::Permutation const &perm) const {
        assert(false);
    }

    TilePimpl scale(const T factor) const {
        return TilePimpl{
            range(), tile_->apply_unary_op(unary_ops::scale_functor(factor)),
            cut()};
    }

    TilePimpl &scale_to(const T factor) {
        tile_->apply_unary_mutation(unary_mutations::scale_functor(factor));
        return *this;
    }

    TilePimpl &subt_to(TilePimpl const &right) {
        tile_->apply_binary_mutation(right.tile(),
                                     binary_mutations::subt_functor(1.0));
        return *this;
    }

    // TODO_TCC figure out what is going on with needing factor second.
    TilePimpl &subt_to(TilePimpl const &right, T factor) {
        convert_to_full(*tile_, right.tile());
        tile_->apply_binary_mutation(right.tile(),
                                     binary_mutations::subt_functor(1.0));
        tile_->apply_unary_mutation(unary_mutations::scale_functor(factor));
        return *this;
    }


    TilePimpl
    subt(TilePimpl const &right, TiledArray::Permutation const &perm) const {
        assert(false);
    }

    TilePimpl subt(TilePimpl const &right) const {
        return TilePimpl{
            range(),
            tile_->apply_binary_op(right.tile(), binary_ops::subt_functor(1.0)),
            std::max(cut(), right.cut())};
    }

    TilePimpl neg(TiledArray::Permutation const &perm) const { assert(false); }
    TilePimpl neg() const { assert(false); }

    TilePimpl &neg_to() { assert(false); }


    template <typename Archive>
    void serialize(Archive &ar) {}

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
        if (result.tag() == TileVariant<T>::FullRank) {
            return;
        }
        auto get_rank = [&](TileVariant<T> const &t) { return t.rank(); };
        const auto out_rank = mutation_rank(result.rank(), get_rank(rest)...);

        if (double(out_rank) > 0.5 * result.full_rank()) {
            result = TileVariant<T>{FullRankTile<T>{result.lrtile().matrix()}};
        }
    }
};

#endif // TCC_MATRIX_TILE_PIMPLE_H
