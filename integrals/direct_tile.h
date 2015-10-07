#pragma once
#ifndef MPQC_INTEGRALS_DIRECTTILE_H
#define MPQC_INTEGRALS_DIRECTTILE_H

#include "task_integrals_common.h"
#include "integral_engine_pool.h"
#include "task_integrals_op_invoker.h"

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
template <typename E, unsigned long N, typename Op, typename Screener>
class DirectTile {
  private:
    TA::Range range_;
    std::vector<std::size_t> idx_;
    ShrPool<E> engines_;
    detail::ShrShellVecArray<N> shell_vecs_;
    std::function<detail::Ttype<Op>(TA::TensorD &&)> op_;
    Screener screen_;

  public:
    using value_type = double;
    using eval_type = detail::Ttype<Op>;
    using range_type = TA::Range;

    DirectTile() = default;
    DirectTile(DirectTile const &) = default;
    DirectTile(DirectTile &&) = default;
    DirectTile &operator=(DirectTile const &) = default;
    DirectTile &operator=(DirectTile &&) = default;

    DirectTile(TA::Range range, std::vector<std::size_t> index,
               ShrPool<E> engines, detail::ShrShellVecArray<N> shells, Op op,
               Screener screen)
            : range_(std::move(range)),
              idx_(std::move(index)),
              engines_(std::move(engines)),
              shell_vecs_(std::move(shells)),
              op_(std::move(op)),
              screen_(std::move(screen)) {}

    operator eval_type() const {
        auto invoker
              = detail::make_screened_op_invoke(idx_, engines_, shell_vecs_,
                                                op_, screen_);
        return invoker(range_);
    }

    template <typename Archive>
    Archive &serialize(Archive &ar) {
        assert(false);
        return ar;
    }
};

template <typename E, unsigned long N, typename Op, typename Screener>
DirectTile<E, N, Op, Screener>
make_direct_tile(TA::Range range, std::vector<std::size_t> index,
                 ShrPool<E> engines, detail::ShrBases<N> const &bases, Op op,
                 Screener screen) {

    auto shell_vecs = detail::get_shells(index, bases);
    return DirectTile<E, N, Op, Screener>(std::move(range), std::move(index),
                                          std::move(engines), shell_vecs,
                                          std::move(op), std::move(screen));
}

} // namespace integrals
} // namespace mpqc

#endif // MPQC_INTEGRALS_DIRECTTILE_H
