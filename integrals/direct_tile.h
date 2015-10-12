#pragma once
#ifndef MPQC_INTEGRALS_DIRECTTILE_H
#define MPQC_INTEGRALS_DIRECTTILE_H

#include "task_integrals_common.h"
#include "integral_engine_pool.h"
#include "task_integrals_op_invoker.h"

#include "integral_screeners.h"

#include "../tensor/tcc_tile.h"

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
template <typename E, unsigned long N, typename Op>
class DirectTile {
  private:
    using TileType = detail::Ttype<Op>;
    using FnType = TileType(TA::Range const&);

    TA::Range range_;
    std::vector<std::size_t> idx_;
    ShrPool<E> engines_;
    detail::ShrShellVecArray<N> shell_vecs_;
    std::function<FnType> op_invoker_;

  public:
    using value_type = double;
    using eval_type = detail::Ttype<Op>;
    using range_type = TA::Range;

    DirectTile() = default;
    DirectTile(DirectTile const &) = default;
    DirectTile(DirectTile &&) = default;
    DirectTile &operator=(DirectTile const &) = default;
    DirectTile &operator=(DirectTile &&) = default;

    template <typename Invoker>
    DirectTile(TA::Range range, std::vector<std::size_t> index,
               ShrPool<E> engines, detail::ShrShellVecArray<N> shells,
               Invoker op_invoke)
            : range_(std::move(range)),
              idx_(std::move(index)),
              engines_(std::move(engines)),
              shell_vecs_(std::move(shells)),
              op_invoker_(std::function<FnType>(std::move(op_invoke))) {}

    operator eval_type() const { return op_invoker_(range_); }

    template <typename Archive>
    Archive &serialize(Archive &ar) {
        assert(false);
        return ar;
    }
};

template <typename ScreenOp = init_base_screen, typename E, unsigned long N,
          typename Op>
DirectTile<E, N, Op>
make_direct_tile(TA::Range range, std::vector<std::size_t> index,
                 ShrPool<E> engines, detail::ShrBases<N> bases, Op op) {

    auto shell_vecs = detail::get_shells(index, bases);
    auto screen = ScreenOp()(index, bases, engines);

    auto op_invoke = detail::make_op_invoke(index, engines, shell_vecs,
                                            std::move(op), std::move(screen));

    return DirectTile<E, N, Op>(std::move(range), std::move(index),
                                std::move(engines), shell_vecs,
                                std::move(op_invoke));
}

} // namespace integrals
} // namespace mpqc

#endif // MPQC_INTEGRALS_DIRECTTILE_H
