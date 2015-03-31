#pragma once
#ifndef TCC_TENSOR_TCCTILE_H
#define TCC_TENSOR_TCCTILE_H

#include "../common/typedefs.h"
#include "../include/tiledarray.h"

#include <memory>

namespace tcc {
namespace tensor {

template <typename T>
class Tile {
    TA::Range range_;
    std::shared_ptr<T> tile_;

  public:
    using eval_type = T;
    using value_type = Tile<T>;
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
        : range_{std::move(r)}, tile_{std::make_shared<T>(std::move(t))} {}
    Tile(TA::Range r, std::shared_ptr<T> st)
        : range_{std::move(r)}, tile_{std::move(st)} {}

    template<typename Value>
    Tile(TA::Range r, Value v)
        : range_{std::move(r)}, tile_{std::make_shared<T>(T{r,v})} {}

    std::shared_ptr<T> tile_ptr() { return tile_; }
    const std::shared_ptr<T> tile_ptr() const { return tile_; }
    T &tile() { return *tile_; }
    T const &tile() const { return *tile_; }
    Tile clone() const { return Tile{range_, *tile_}; }
    TA::Range const &range() const { return range; }
    bool empty() const { return !tile_; }
    double norm() const { return tile_->norm(); }

    template <typename Archive>
    void serialize(Archive &ar) {
        assert(false);
    }
};

template<typename T>
inline std::ostream& operator<<(std::ostream& os, Tile<T> const & tile){
    os << tile.range() << ": " << tile.tile() << "\n";
    return os;
}

} // namespace tensor
} // namespace tcc

#endif // TCC_TENSOR_TCCTILE_H
