
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_DIRECT_TILE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_DIRECT_TILE_H_

#include "mpqc/chemistry/qc/lcao/integrals/task_integrals_common.h"
#include "mpqc/chemistry/qc/lcao/integrals/integral_builder.h"

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
template <typename Tile, typename Engine = libint2::Engine>
class DirectTile {
 public:
  using Builder = DirectIntegralBuilder<Tile, Engine>;
  using BaseBuilder = IntegralBuilder<Tile, Engine>;

 private:
  std::vector<std::size_t> idx_;
  TA::Range range_;
  std::shared_ptr<Builder> builder_;

 public:
  using value_type = double;
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
        world->template shared_ptr_from_id<BaseBuilder>(id));
    assert(builder_ != nullptr);
  }
};

/*! Class to hold a direct tile builder with its array. */
template <typename Tile, typename Policy, typename Engine = libint2::Engine>
class DirectArray {
 public:
  using Builder = DirectIntegralBuilder<Tile, Engine>;
  using Array = TA::DistArray<DirectTile<Tile, Engine>, Policy>;

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
  auto operator()(Args &&... args)
      -> decltype(array_(std::forward<Args>(args)...)) {
    return array_(std::forward<Args>(args)...);
  }

  template <typename... Args>
  auto operator()(Args &&... args) const
      -> decltype(array_(std::forward<Args>(args)...)) {
    return array_(std::forward<Args>(args)...);
  }

  void set_array(Array a) { array_ = std::move(a); }

  Array &array() { return array_; }

  Array const &array() const { return array_; }

  std::shared_ptr<Builder> builder() { return builder_; }

  const std::shared_ptr<Builder> builder() const { return builder_; }
};

}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_DIRECT_TILE_H_
