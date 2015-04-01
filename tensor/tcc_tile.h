#pragma once
#ifndef TCC_TENSOR_TCCTILE_H
#define TCC_TENSOR_TCCTILE_H

#include "../common/typedefs.h"
#include "../include/tiledarray.h"

#include <memory>

namespace tcc {
namespace tensor {

namespace detail {

// Helper class to facilitate ADL
template <typename T>
class TileModel {
    T tile_;

  public:
    using eval_type = T;
    using value_type = T;
    using numeric_type = typename T::numeric_type;
    using size_type = std::size_t;

    TileModel() = default;
    ~TileModel() = default;
    TileModel(TileModel const &) = default;
    TileModel(TileModel &&) = default;
    TileModel &operator=(TileModel &&) = default;
    TileModel &operator=(TileModel const &) = default;

    TileModel(T &&t) : tile_{std::move(t)} {}
    TileModel(T const &t) : tile_{t} {}

    T &tile() { return tile_; }
    T const &tile() const { return tile_; }

    bool empty_() const { return empty(tile_); }
    auto norm_() const -> decltype(norm(tile_)) { return norm(tile_); }

    T clone_() const { return clone(tile_); }

    template <typename Right>
    auto add_(Right const &r) const -> decltype(add(tile_, r)) {
        return add(tile_, r);
    }

    template <typename Right>
    auto add_(Right const &r,
              TA::Permutation const &p) const -> decltype(add(tile_, r, p)) {
        return add(tile_, r, p);
    }

    template <typename Right>
    T &add_to_(Right const &r) {
        add_to(tile_, r);
        return tile_;
    }

    template <typename Right>
    T &add_to_(Right const &r, TA::Permutation const &p) {
        add_to(tile_, r, p);
        return tile_;
    }

    T permute_(TA::Permutation const &p) const { return permute(tile_, p); }
};

} // namespace detail


template <typename T>
class Tile {
    TA::Range range_;
    std::shared_ptr<detail::TileModel<T>> tile_;

  public:
    using value_type = T;
    using range_type = TA::Range;
    using numeric_type = typename T::numeric_type;
    using size_type = std::size_t;

    Tile() = default;
    ~Tile() = default;
    Tile(Tile const &) = default;
    Tile(Tile &&) = default;
    Tile &operator=(Tile &&) = default;
    Tile &operator=(Tile const &) = default;

    Tile(TA::Range r) : range_{std::move(r)} {}
    Tile(TA::Range r, T t)
            : range_{std::move(r)},
              tile_{std::make_shared<detail::TileModel<T>>(
                    detail::TileModel<T>{std::move(t)})} {}

    template <typename Value>
    Tile(TA::Range r, Value v)
            : range_{std::move(r)},
              tile_{std::make_shared<detail::TileModel<T>>(
                    detail::TileModel<T>{detail::TileModel<T>{T{r, v}}})} {}

    T &tile() { return tile_->tile(); }
    T const &tile() const { return tile_->tile(); }

    Tile clone() const { return Tile{range_, tile_->clone_()}; }

    TA::Range const &range() const { return range_; }

    bool empty() const { return (!tile_ || tile_->empty_()); }

    auto norm() const -> decltype(tile_ -> norm_()) { return tile_->norm_(); }

    template <typename U>
    auto add(
          Tile<U> const &u) const -> Tile<decltype(tile_ -> add_(u.tile()))> {
        auto out = tile_->add_(u.tile());
        return Tile<decltype(out)>{range_, std::move(out)};
    }

    template <typename U>
    auto add(Tile<U> const &u, TA::Permutation const &p) const
          -> Tile<decltype(tile_ -> add_(u.tile(), p))> {
        auto out = tile_->add_(u.tile(), p);
        return Tile<decltype(out)>{range_, std::move(out)};
    }

    template <typename U>
    Tile &add_to(Tile<U> const &u) {
        tile_->add_to_(u.tile());
        return *this;
    }

    Tile permute(TA::Permutation const &p) const {
        return Tile{range_, tile_->permute_(p)};
    }

    template <typename Archive>
    void serialize(Archive &ar) {
        assert(false);
    }
};

template <typename T>
inline std::ostream &operator<<(std::ostream &os, Tile<T> const &tile) {
    os << tile.range() << ": " << tile.tile() << "\n";
    return os;
}

} // namespace tensor
} // namespace tcc

#endif // TCC_TENSOR_TCCTILE_H
