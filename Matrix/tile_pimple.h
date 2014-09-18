#ifndef TCC_MATRIX_TILE_PIMPLE_H
#define TCC_MATRIX_TILE_PIMPLE_H

#include "tile_variant.h"
#include "../include/tiledarray.h"
#include <memory>

template <typename T>
class TilePimple {
  public:
    typedef TilePimple eval_type;
    typedef T value_type;
    typedef TiledArray::Range range_type;
    typedef T numeric_type;
    typedef std::size_t size_type;


  public:
    TilePimple() = default;
    ~TilePimple() = default;
    TilePimple(TilePimple const &t) = default;
    TilePimple &operator=(TilePimple const &t) = default;
    TilePimple(TilePimple &&t) = default;
    TilePimple &operator=(TilePimple &t) = default;

    /*
     * User defined constructors
     */

    // Just range ctor
    TilePimple(TiledArray::Range r) : tile_(), range_(std::move(r)) {}

    TilePimple(TiledArray::Range &&r) : tile_(), range_(std::move(r)) {}

    TilePimple(TiledArray::Range r, double cut)
        : tile_(), range_(std::move(r)), cut_(cut) {}

    TilePimple(TiledArray::Range &&r, double cut)
        : tile_(), range_(std::move(r)), cut_(cut) {}

    // TileVariant ctors
    TilePimple(TiledArray::Range r, TileVariant<T> t)
        : tile_(std::make_shared<TileVariant<T>>(std::move(t))),
          range_(std::move(r)) {}

    TilePimple(TiledArray::Range r, TileVariant<T> t, double cut)
        : tile_(std::make_shared<TileVariant<T>>(std::move(t))),
          range_(std::move(r)), cut_(cut) {}

    TilePimple(TiledArray::Range &&r, TileVariant<T> &&t) noexcept
        : tile_(std::make_shared<TileVariant<T>>(std::move(t))),
          range_(std::move(r)) {}

    TilePimple(TiledArray::Range &&r, TileVariant<T> &&t, double cut) noexcept
        : tile_(std::make_shared<TileVariant<T>>(std::move(t))),
          range_(std::move(r)),
          cut_(cut) {}

    // Clone will make a deep copy
    TilePimple clone() const { return TilePimple(range_, *tile_, cut_); }

    /*
     * Tile information functions
     */
    bool isFull() const { return tile_->tag(); }
    TiledArray::Range const & range() const { return range_; }

    // maybe expensive
    void setCut(double cut) {
        const auto temp = cut_;
        cut_ = cut;
        if (temp <= cut_) {
            // TODO recompress may save some space.
        }
    }

    template <typename Archive>
    void serialize(Archive &ar) {}

  private:
    std::shared_ptr<TileVariant<T>> tile_;
    TiledArray::Range range_;
    double cut_ = 1e-7;
};

#endif // TCC_MATRIX_TILE_PIMPLE_H
