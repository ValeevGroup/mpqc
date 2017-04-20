
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_DIRECT_TILE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_DIRECT_TILE_H_

#include "mpqc/chemistry/qc/lcao/integrals/integral_builder.h"
#include "mpqc/chemistry/qc/lcao/integrals/task_integrals_common.h"

#include <madness/world/worldptr.h>

#include <functional>
#include <memory>
#include <vector>

#include <tiledarray.h>

namespace mpqc {
namespace lcao {
namespace gaussian {

/*! \brief A direct tile for integral construction
 *
 */
template <typename Tile, typename Builder = DirectIntegralBuilder<Tile, libint2::Engine> >
class DirectTile {
 private:
  std::vector<std::size_t> idx_;
  TA::Range range_;
  std::shared_ptr<Builder> builder_;

 public:
  using value_type = typename Tile::numeric_type;
  using eval_type = Tile;
  using range_type = TA::Range;

  DirectTile() = default;
  DirectTile(DirectTile const &) = default;
  DirectTile(DirectTile &&) = default;
  DirectTile &operator=(DirectTile const &) = default;
  DirectTile &operator=(DirectTile &&) = default;

  DirectTile(std::vector<std::size_t> index, TA::Range range,
             std::shared_ptr<Builder> builder)
      : idx_(std::move(index)),
        range_(std::move(range)),
        builder_(std::move(builder)) {}

  operator eval_type() const { return builder_->operator()(idx_, range_); }

  /*! \brief Allows for truncate to be called on direct tiles
   *
   * \note In reduced scaling code this could lead to expensive higher order
   * operations depending on the size of the array.
   */
  value_type norm() const {
    auto tile = builder_->operator()(idx_, range_);
    return tile.norm();
  }

  template <typename Archive>
  std::enable_if_t<madness::archive::is_output_archive<Archive>::value, void>
  serialize(Archive &ar) {
    ar &idx_;
    ar &range_;
    assert(builder_ != nullptr);
    ar & builder_->id();
  }

  template <typename Archive>
  std::enable_if_t<madness::archive::is_input_archive<Archive>::value, void>
  serialize(Archive &ar) {
    ar &idx_;
    ar &range_;
    madness::uniqueidT id;
    ar &id;

    assert(builder_ == nullptr);
    madness::World *world = madness::World::world_from_id(id.get_world_id());
    // dynamic cast
    builder_ = std::dynamic_pointer_cast<Builder>(
        world->template shared_ptr_from_id<Builder>(id));
    assert(builder_ != nullptr);
  }
};

/*! Class to hold a direct tile builder with its array. */
template <typename Tile, typename Policy, typename Builder>
class DirectArray {
 public:
  using Array = TA::DistArray<DirectTile<Tile, Builder>, Policy>;

 private:
  std::shared_ptr<Builder> builder_;
  Array array_;

 public:
  using tile_type = typename Array::value_type;

 public:
  DirectArray() = default;
  DirectArray(std::shared_ptr<Builder> b, Array a)
      : builder_(std::move(b)), array_(std::move(a)) {}

  DirectArray(std::shared_ptr<Builder> b) : builder_(std::move(b)), array_() {}

  bool is_initialized() const { return builder_ && array_.is_initialized(); }

  template <typename... Args>
  auto operator()(Args &&... args) {
    return array_(std::forward<Args>(args)...);
  }

  template <typename... Args>
  auto operator()(Args &&... args) const {
    return array_(std::forward<Args>(args)...);
  }

  void set_array(Array a) { array_ = std::move(a); }

  Array &array() { return array_; }

  Array const &array() const { return array_; }

  std::shared_ptr<Builder> builder() { return builder_; }

  const std::shared_ptr<Builder> builder() const { return builder_; }
};

/*
 * A Direct Integral Tile that computes Tile on the fly using Density Fitting
 * three center integral
 */
template <typename Tile, typename Policy>
class DirectDFTile {
 public:
  using Builder = DirectDFIntegralBuilder<Tile, Policy>;
  using value_type = typename Tile::numeric_type;
  using eval_type = Tile;
  using range_type = TA::Range;

  DirectDFTile() = default;
  DirectDFTile(DirectDFTile const &) = default;
  DirectDFTile(DirectDFTile &&) = default;
  DirectDFTile &operator=(DirectDFTile const &) = default;
  DirectDFTile &operator=(DirectDFTile &&) = default;

  DirectDFTile(std::vector<std::size_t> index, TA::Range range,
               std::shared_ptr<Builder> builder)
      : idx_(std::move(index)),
        range_(std::move(range)),
        builder_(std::move(builder)) {}

  /// compute and return Tile
  explicit operator TA::Future<eval_type>() const& { return builder_->operator()(idx_, range_); }
//  explicit operator eval_type () const& { return builder_->operator()(idx_, range_); }

  /// output serialize
  template <typename Archive>
  std::enable_if_t<madness::archive::is_output_archive<Archive>::value, void>
  serialize(Archive &ar) {
    ar &idx_;
    ar &range_;
    TA_ASSERT(builder_ != nullptr);
    ar & builder_->id();
  }

  /// input serialize
  template <typename Archive>
  std::enable_if_t<madness::archive::is_input_archive<Archive>::value, void>
  serialize(Archive &ar) {
    ar &idx_;
    ar &range_;
    madness::uniqueidT id;
    ar &id;

    TA_ASSERT(builder_ == nullptr);
    madness::World *world = madness::World::world_from_id(id.get_world_id());
    // dynamic cast
    builder_ = std::dynamic_pointer_cast<Builder>(
        world->template shared_ptr_from_id<Builder>(id));
    TA_ASSERT(builder_ != nullptr);
  }

 private:
  std::vector<std::size_t> idx_;
  TA::Range range_;
  std::shared_ptr<Builder> builder_;
};

}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_DIRECT_TILE_H_
