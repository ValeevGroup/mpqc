#pragma once
#ifndef TCC_TENSOR_TILEPIMPL_DEVEL_H
#define TCC_TENSOR_TILEPIMPL_DEVEL_H

#include "tile_variant_devel.h"
#include "../include/tiledarray.h"

#include <memory>

namespace tcc {
namespace tensor {

template <typename T, typename Range = TiledArray::Range>
class TilePimplDevel {
  public:
    typedef TilePimplDevel eval_type;
    typedef T value_type;
    typedef Range range_type;
    typedef T numeric_type;
    typedef std::size_t size_type;

  public:
    TilePimplDevel() = default;
    ~TilePimplDevel() = default;
    TilePimplDevel(TilePimplDevel const &t) = default;
    TilePimplDevel &operator=(TilePimplDevel const &t) = default;
    TilePimplDevel(TilePimplDevel &&t) = default;
    TilePimplDevel &operator=(TilePimplDevel &t) = default;

    TilePimplDevel(Range r) : tile_(), pimpl_range_(std::move(r)) {}

    TilePimplDevel(Range r, TileVariantDevel<T> t)
        : tile_(std::make_shared<TileVariantDevel<T>>(std::move(t))),
          pimpl_range_(std::move(r)) {}

    TilePimplDevel clone() const {
        return TilePimplDevel{pimpl_range_, *tile_};
    }

    TileVariantDevel<T> const &tile() const {
        assert(tile_);
        return *tile_;
    }

    bool empty() const { return !tile_; }
    Range const &range() const { return pimpl_range_; }

    template <typename Archive>
    void serialize(Archive &ar) {
        assert(false);
    }

  private:
    std::shared_ptr<TileVariantDevel<T>> tile_;
    Range pimpl_range_;
};

} // namespace tensor
} // namespace tcc


#endif /* end of include guard: TCC_TENSOR_TILEPIMPL_DEVEL_H */
