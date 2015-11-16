#pragma once
#ifndef MPQC_INTEGRALS_DIRECTTILE_H
#define MPQC_INTEGRALS_DIRECTTILE_H

#include "task_integrals_common.h"

#include "../tensor/tcc_tile.h"

#include "../include/tiledarray.h"

#include <madness/world/worldptr.h>

#include <memory>
#include <vector>
#include <functional>

namespace mpqc {
namespace integrals {

namespace detail {

template <typename T>
using ShrTile = tcc::tensor::Tile<T>;

} // namespace detail

/*! \brief A direct tile for integral construction
 *
 */
template <typename Builder>
class DirectTile {
  private:
    std::vector<std::size_t> idx_;
    TA::Range range_;
    std::shared_ptr<Builder> builder_;
    madness::detail::WorldPtr<madness::World> world_ptr_;
    madness::uniqueidT builder_id_;

    using TileType = decltype(std::declval<Builder>()(
          std::declval<std::vector<std::size_t>>(), std::declval<TA::Range>()));

  public:
    using value_type = double;
    using eval_type = TileType;
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
              builder_(std::move(builder)),
              world_ptr_(builder_->get_world(), &builder_->get_world()),
              builder_id_(builder_->id()) {}

    operator eval_type() const { return builder_->operator()(idx_, range_); }

    template <typename Archive>
    void serialize(Archive &) {
        assert(false);
    }

#if 0
    template <typename Archive>
    enable_if_t<madness::archive::is_output_archive<Archive>::value, void>
    serialize(Archive &ar) {
        ar &idx_;
        ar &range_;
        ar &world_ptr_;
        ar &builder_id_;
    }

    template <typename Archive>
    enable_if_t<madness::archive::is_input_archive<Archive>::value, void>
    serialize(Archive &ar) {
        ar &idx_;
        ar &range_;
        ar &world_ptr_;
        ar & builder_id_;

        assert(builder_ == nullptr);
        builder_ = world_ptr_.get_world().ptr_from_id<Builder>(builder_id_);
    }
#endif
};

} // namespace integrals
} // namespace mpqc

#endif // MPQC_INTEGRALS_DIRECTTILE_H
