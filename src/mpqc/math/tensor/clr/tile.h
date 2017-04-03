#ifndef MPQC4_SRC_MPQC_MATH_TENSOR_CLR_TILE_H_
#define MPQC4_SRC_MPQC_MATH_TENSOR_CLR_TILE_H_

#include <memory>

#include <tiledarray.h>

#include "mpqc/math/tensor/clr/decomposed_tensor_unary.h"
#include "mpqc/util/meta/get_type.h"

namespace mpqc {
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
  ~TileModel() { }
  TileModel(TileModel const &) = default;
  TileModel(TileModel &&) = default;
  TileModel &operator=(TileModel &&) = default;
  TileModel &operator=(TileModel const &) = default;

  TileModel(T &&t) : tile_{std::move(t)} {}
  TileModel(T const &t) : tile_{t} {}

  T &tile() { return tile_; }
  T const &tile() const { return tile_; }

  bool empty_() const { return empty(tile_); }
  auto norm_() const { return norm(tile_); }
  auto trace_() const { return trace(tile_); }
  auto sum_() const { return sum(tile_); }
  auto min_() const { return min(tile_); }
  auto max_() const { return max(tile_); }
  auto abs_min_() const { return abs_min(tile_); }
  auto abs_max_() const { return abs_max(tile_); }
  auto product_() const { return product(tile_); }
  auto squared_norm_() const {
    return squared_norm(tile_);
  }

  T clone_() const { return clone(tile_); }

  T permute_(TA::Permutation const &p) const { return permute(tile_, p); }

  /*
   * Add
   */
  template <typename... Args>
  auto add_(Args &&... args) const {
    return add(tile_, std::forward<Args>(args)...);
  }

  template <typename... Args>
  T &add_to_(Args &&... args) {
    add_to(tile_, std::forward<Args>(args)...);
    return tile_;
  }

  /*
   * Subtract
   */
  template <typename... Args>
  auto subt_(Args &&... args) const {
    return subt(tile_, std::forward<Args>(args)...);
  }

  template <typename... Args>
  T &subt_to_(Args &&... args) {
    subt_to(tile_, std::forward<Args>(args)...);
    return tile_;
  }

  /*
   * Multiply
   */
  template <typename... Args>
  auto mult_(Args &&... args) const {
    return mult(tile_, std::forward<Args>(args)...);
  }

  template <typename... Args>
  T &mult_to_(Args &&... args) {
    mult_to(tile_, std::forward<Args>(args)...);
    return tile_;
  }

  /*
   * Negate
   */
  template <typename... Args>
  auto neg_(Args &&... args) const {
    return neg(tile_, std::forward<Args>(args)...);
  }

  template <typename... Args>
  T &neg_to_(Args &&... args) {
    neg_to(tile_, std::forward<Args>(args)...);
    return tile_;
  }

  /*
   * Scale
   */
  template <typename... Args>
  auto scale_(Args &&... args) const {
    return scale(tile_, std::forward<Args>(args)...);
  }

  template <typename... Args>
  T &scale_to_(Args &&... args) {
    scale_to(tile_, std::forward<Args>(args)...);
    return tile_;
  }

  /*
   * Gemm
   */
  template <typename... Args>
  auto gemm_(Args &&... args) const {
    return gemm(tile_, std::forward<Args>(args)...);
  }

  template <typename... Args>
  T &gemm_to_(Args &&... args) {
    gemm(tile_, std::forward<Args>(args)...);
    return tile_;
  }
};

}  // namespace detail

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
  ~Tile() { }
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
      : range_{r},
        tile_{std::make_shared<detail::TileModel<T>>(
            detail::TileModel<T>(detail::TileModel<T>(T(r, v))))} {}

 public:  // Member functions for general use
  T &tile() { return tile_->tile(); }
  T const &tile() const { return tile_->tile(); }

  Tile clone() const { return Tile{range_, tile_->clone_()}; }

  TA::Range const &range() const { return range_; }

  bool empty() const { return (!tile_ || tile_->empty_()); }

  auto norm() const { return tile_->norm_(); }

  auto trace() const { return tile_->trace_(); }

  auto sum() const { return tile_->sum_(); }

  auto min() const { return tile_->min_(); }

  auto max() const { return tile_->max_(); }

  auto abs_min() const {
    return tile_->abs_min_();
  }

  auto abs_max() const {
    return tile_->abs_max_();
  }

  auto product() const {
    return tile_->product_();
  }

  auto squared_norm() {
    return tile_->squared_norm_();
  }

 private:
  template <typename... Args,
            typename std::enable_if<
                std::is_same<typename utility::meta::last_type<Args...>::type,
                             TA::Permutation const &>::value>::type * = nullptr>
  TA::Range create_new_range(TA::Range const &r, Args &&... args) const {
    auto const &perm = utility::meta::back(args...);
    return perm * r;
  }

  template <typename... Args,
            typename std::enable_if<!std::is_same<
                typename utility::meta::last_type<Args...>::type,
                TA::Permutation const &>::value>::type * = nullptr>
  TA::Range create_new_range(TA::Range const &r, Args &&...) const {
    return r;
  }

 public:  // TA math functions
  Tile permute(TA::Permutation const &p) const {
    return Tile{p * range_, tile_->permute_(p)};
  }

  template <typename Archive>
  typename std::enable_if<
      madness::archive::is_output_archive<Archive>::value>::type
  serialize(Archive &ar) {
    ar &range_;
    bool empty = !static_cast<bool>(tile_);
    ar &empty;
    if (!empty) {
      ar & tile_->tile();
    }
  }

  template <typename Archive>
  typename std::enable_if<
      madness::archive::is_input_archive<Archive>::value>::type
  serialize(Archive &ar) {
    ar &range_;
    bool empty = false;
    ar &empty;
    if (!empty) {
      T tile;
      ar &tile;
      tile_ = std::make_shared<detail::TileModel<T>>(
          detail::TileModel<T>{std::move(tile)});
    }
  }

  /*
   * Add functions
   */
  template <typename... Args>
  auto add(Args &&... args) const {
    auto range = create_new_range(range_, args...);
    auto out = tile_->add_(std::forward<Args>(args)...);
    return Tile<decltype(out)>{std::move(range), std::move(out)};
  }

  // Overload for tcc::Tile
  template <typename U, typename... Args>
  auto add(Tile<U> const &u, Args &&... args) const {
    auto range = create_new_range(range_, args...);
    auto out = tile_->add_(u.tile(), std::forward<Args>(args)...);
    return Tile<decltype(out)>{std::move(range), std::move(out)};
  }

  template <typename... Args>
  Tile &add_to(Args &&... args) {
    tile_->add_to_(std::forward<Args>(args)...);
    return *this;
  }

  template <typename U, typename... Args>
  Tile &add_to(Tile<U> const &u, Args &&... args) {
    tile_->add_to_(u.tile(), std::forward<Args>(args)...);
    return *this;
  }

  /*
   * Subtract functions
   */
  template <typename... Args>
  auto subt(Args &&... args) const {
    auto range = create_new_range(range_, args...);
    auto out = tile_->subt_(std::forward<Args>(args)...);
    return Tile<decltype(out)>{std::move(range), std::move(out)};
  }

  // Overload for tcc::Tile
  template <typename U, typename... Args>
  auto subt(Tile<U> const &u, Args &&... args) const {
    auto range = create_new_range(range_, args...);
    auto out = tile_->subt_(u.tile(), std::forward<Args>(args)...);
    return Tile<decltype(out)>{std::move(range), std::move(out)};
  }

  template <typename... Args>
  Tile &subt_to(Args &&... args) {
    tile_->subt_to_(std::forward<Args>(args)...);
    return *this;
  }

  template <typename U, typename... Args>
  Tile &subt_to(Tile<U> const &u, Args &&... args) {
    tile_->subt_to_(u.tile(), std::forward<Args>(args)...);
    return *this;
  }

  /*
   * Multiplication functions
   */
  template <typename... Args>
  auto mult(Args &&... args) const {
    auto range = create_new_range(range_, args...);
    auto out = tile_->mult_(std::forward<Args>(args)...);
    return Tile<decltype(out)>{std::move(range), std::move(out)};
  }

  // Overload for tcc::Tile
  template <typename U, typename... Args>
  auto mult(Tile<U> const &u, Args &&... args) const {
    auto range = create_new_range(range_, args...);
    auto out = tile_->mult_(u.tile(), std::forward<Args>(args)...);
    return Tile<decltype(out)>{std::move(range), std::move(out)};
  }

  template <typename... Args>
  Tile &mult_to(Args &&... args) {
    tile_->mult_to_(std::forward<Args>(args)...);
    return *this;
  }

  template <typename U, typename... Args>
  Tile &mult_to(Tile<U> const &u, Args &&... args) {
    tile_->mult_to_(u.tile(), std::forward<Args>(args)...);
    return *this;
  }

  /*
   * Neg functions
   */
  template <typename... Args>
  Tile neg(Args &&... args) const {
    TA::Range range = create_new_range(range_, args...);
    return Tile{std::move(range), tile_->neg_(std::forward<Args>(args)...)};
  };

  template <typename... Args>
  Tile scale(Args &&... args) const {
    TA::Range range = create_new_range(range_, args...);
    return Tile{std::move(range), tile_->scale_(std::forward<Args>(args)...)};
  };

  template <typename... Args>
  Tile &scale_to(Args &&... args) {
    tile_->scale_to_(std::forward<Args>(args)...);
    return *this;
  };

  Tile &neg_to() {
    tile_->neg_to_();
    return *this;
  };

  /*
   * Gemm functions
   */
  template <typename Other>
  auto gemm(Other const &o, const numeric_type factor,
            TA::math::GemmHelper const &gemm_helper) const {
    auto range = gemm_helper.make_result_range<TA::Range>(range_, o.range());

    auto out = tile_->gemm_(o, factor, gemm_helper);
    return Tile<decltype(out)>{std::move(range), std::move(out)};
  }

  template <typename U>
  auto gemm(Tile<U> const &u, const numeric_type factor,
            TA::math::GemmHelper const &gemm_helper) const {
    auto range = gemm_helper.make_result_range<TA::Range>(range_, u.range());
    auto out = tile_->gemm_(u.tile(), factor, gemm_helper);
    return Tile<decltype(out)>{std::move(range), std::move(out)};
  }

  template <typename Left, typename Right>
  Tile &gemm(Left const &l, Right const &r, const numeric_type factor,
             TA::math::GemmHelper const &gemm_helper) {
    tile_->gemm_to_(l, r, factor, gemm_helper);
    return *this;
  }

  template <typename U, typename V>
  Tile &gemm(Tile<U> const &l, Tile<V> const &r, const numeric_type factor,
             TA::math::GemmHelper const &gemm_helper) {
    tile_->gemm_to_(l.tile(), r.tile(), factor, gemm_helper);
    return *this;
  }
};

template <typename T>
inline std::ostream &operator<<(std::ostream &os, Tile<T> const &tile) {
  os << tile.range() << ": " << tile.tile();
  return os;
}

template <typename T>
unsigned long tile_clr_storage(Tile<DecomposedTensor<T>> const &tile) {
  auto size = 0ul;
  for (auto const &t : tile.tile().tensors()) {
    size += t.range().volume();
  }
  return size;
}

}  // namespace tensor
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_MATH_TENSOR_CLR_TILE_H_
