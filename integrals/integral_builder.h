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
    detail::ShrBases<N> bases_;
    ShrPool<E> engines_;
    std::shared_ptr<Screener> screen_;
    Op op_;

  public:
    using op_type = detail::Ttype<Op>;

    /*! \brief Constructor which copies all shared_ptr members
     *
     * \param world Should be the same world as the one for the array which will
     * hold the tiles.
     * \param shr_epool is a shared pointer to an IntegralEnginePool
     * \param shr_bases is a shared pointer to an array of Basis
     * \param screen is a shared pointer to a Screener type
     * \param op should be a thread safe function or functor that takes a
     *  rvalue of a TA::TensorD and returns a valid TA::Array tile.
     */
    IntegralBuilder(madness::World &world, ShrPool<E> shr_epool,
                    detail::ShrBases<N> shr_bases,
                    std::shared_ptr<Screener> screen, Op op)
            : madness::WorldObject<IntegralBuilder<N, E, Op>>(world),
              bases_(std::move(shr_bases)),
              engines_(std::move(shr_epool)),
              screen_(std::move(screen)),
              op_(std::move(op)) {
        // Must call for WorldObject Interface to be satisfied
        this->process_pending();
    }


    op_type operator()(std::vector<std::size_t> const &idx, TA::Range range) {
        return op_(integrals(idx, std::move(range)));
    }

    TA::TensorD
    integrals(std::vector<std::size_t> const &idx, TA::Range range) {

        // Get integral shells
        detail::VecArray<N> shellvec_ptrs;
        for (auto i = 0ul; i < N; ++i) {
            auto const &basis_i = bases_->operator[](i);
            shellvec_ptrs[i] = &basis_i.cluster_shells()[idx[i]];
        }

        // Compute integrals over the selected shells.
        return detail::integral_kernel(engines_->local(), std::move(range),
                                       shellvec_ptrs, *screen_);
    }

    op_type op(TA::TensorD &&tensor) { return op_(std::move(tensor)); }
};

/*!
 * \brief Function to make detection of template parameters easier, see
 * IntegralBuilder for details.
 */
template <typename E, typename Op, unsigned long N>
IntegralBuilder<N, E, Op>
make_integral_builder(madness::World &world, ShrPool<E> shr_epool,
                      detail::ShrBases<N> shr_bases,
                      std::shared_ptr<Screener> shr_screen, Op op) {
    return IntegralBuilder<N, E, Op>(world, std::move(shr_epool),
                                     std::move(shr_bases),
                                     std::move(shr_screen), std::move(op));
}

} // namespace integrals
} // namespace mpqc
#endif // MPQC_INTEGRALS_INTEGRAL_BUILDER_H
