#pragma once
#ifndef MPQC_INTEGRALS_INTEGRAL_BUILDER_H
#define MPQC_INTEGRALS_INTEGRAL_BUILDER_H

#include "../common/typedefs.h"
#include "../basis/basis.h"

#include "task_integrals_common.h"
#include "screening/screen_base.h"

#include "task_integral_kernels.h"
#include "integral_engine_pool.h"

#include "../include/tiledarray.h"

#include <array>

namespace mpqc {
namespace integrals {

/*! \brief Builds integrals from an  array of bases and an integral engine pool.
 *
 * \param Op is a function or functor that takes a TA::Tensor && and returns a
 *tile.
 * The simplest type of Op is simply the following:
 * ```
 * auto op = [](TA::Tensor<double> && t){ return std::move(t) };
 * ```
 */
template <unsigned long N, typename E, typename Op>
class IntegralBuilder : public madness::WorldObject<IntegralBuilder<N, E, Op>> {
  private:
    std::array<basis::Basis, N> bases_;
    Epool<E> engines_;
    std::unique_ptr<Screener> screen_;
    Op op_;

  public:
    using op_type = detail::Ttype<Op>;

    /*! \brief Only constructor for IntegralBuilder
     *
     * \param world Should be the same world as the one for the array which will
     * hold the tiles.
     *
     * \param screen is a screening type that inherited from Screener.
     */
    template <typename ScreenerType>
    IntegralBuilder(madness::World &world, E engine,
                    std::array<basis::Basis, N> const &bases,
                    ScreenerType screen, Op op)
            : madness::WorldObject<IntegralBuilder<N, E, Op>>(world),
              bases_(bases),
              engines_(Epool<E>(std::move(engine))),
              screen_(
                    std::unique_ptr<Screener>(new Screener(std::move(screen)))),
              op_(std::move(op)) {
                  // Must call for WorldObject Interface to be satisfied
                  this->process_pending();
              }


    op_type
    operator()(std::vector<std::size_t> const &idx, TA::Range range) {
        return op_(integrals(std::move(idx), std::move(range)));
    }

    TA::TensorD
    integrals(std::vector<std::size_t> const &idx, TA::Range range) {

        // Get integral shells
        detail::VecArray<N> shellvec_ptrs;
        for (auto i = 0ul; i < N; ++i) {
            shellvec_ptrs[i] = &bases_[i].cluster_shells()[idx[i]];
        }

        // Compute integrals over the selected shells.
        return detail::integral_kernel(engines_.local(), std::move(range),
                               shellvec_ptrs, *screen_);
    }

    op_type op(TA::TensorD &&tensor) { return op_(std::move(tensor)); }
};

/*!
 * \brief Function to make detection of template parameters easier, see
 * IntegralBuilder for details.
 */
template <typename E, typename Op, unsigned long N, typename ScreenerType>
IntegralBuilder<N, E, Op>
make_integral_builder(madness::World &world, E const &engine,
                      std::array<basis::Basis, N> const &bases, 
                      ScreenerType screen, Op op) {
    return IntegralBuilder<N, E, Op>(world, engine, bases, std::move(screen),
                                     std::move(op));
}

} // namespace integrals
} // namespace mpqc
#endif // MPQC_INTEGRALS_INTEGRAL_BUILDER_H
