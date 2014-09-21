#ifndef TCC_MATRIX_TILE_PIMPLE_H
#define TCC_MATRIX_TILE_PIMPLE_H

#include "tile_variant.h"
#include "tile_ops.h"
#include "../include/tiledarray.h"
#include <memory>

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

    TilePimpl(TiledArray::Range &&r) : tile_(), range_(std::move(r)) {}

    TilePimpl(TiledArray::Range r, double cut)
        : tile_(), range_(std::move(r)), cut_(cut) {}

    TilePimpl(TiledArray::Range &&r, double cut)
        : tile_(), range_(std::move(r)), cut_(cut) {}

    // TileVariant ctors
    TilePimpl(TiledArray::Range r, TileVariant<T> t)
        : tile_(std::make_shared<TileVariant<T>>(std::move(t))),
          range_(std::move(r)) {}

    explicit TilePimpl(TiledArray::Range r, TileVariant<T> t, double cut)
        : tile_(std::make_shared<TileVariant<T>>(std::move(t))),
          range_(std::move(r)), cut_(cut) {}

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

    // maybe expensive
    void setCut(double cut) {
        const auto temp = cut_;
        cut_ = cut;
        if (temp <= cut_) {
            // TODO_TCC recompress may save some space.
        }
    }

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

        // TODO complete this section.
        // auto result_range = gemm_config.make_result_range
        //                    <range_type>(range(), right.range());
        auto result_range = range();

        return TilePimpl(std::move(result_range),
                         std::move(tile_->apply_binary_op(
                             right.tile(), tile_ops::gemm_AB(factor))),
                         std::max(cut(), right.cut()));
    }

    TilePimpl &gemm(TilePimpl const &left, TilePimpl const &right,
                    numeric_type factor,
                    TiledArray::math::GemmHelper const &gemm_config) {

        // TODO complete this section.
        // auto result_range = gemm_config.make_result_range
        //                    <range_type>(range(), right.range());
        auto result_range = range();

        tile_->apply_binary_op_to(left.tile(), right.tile(),
                                  tile_ops::gemm_inplace(factor));

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

    TilePimpl
    mult(TilePimpl const &right) const {
        assert(false);
        return TilePimpl();
    }

    TilePimpl &
    mult_to(TilePimpl const &right) {
        assert(false);
        return *this;
    }


    TilePimpl permute(TiledArray::Permutation const &perm) const {
        assert(false);
        return TilePimpl();
    }


    template <typename Archive>
    void serialize(Archive &ar) {}
};

#endif // TCC_MATRIX_TILE_PIMPLE_H
