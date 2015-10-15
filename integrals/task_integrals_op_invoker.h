#pragma once
#ifndef MPQC_INTEGRALS_TASKINTEGRALSOPINVOKER_H
#define MPQC_INTEGRALS_TASKINTEGRALSOPINVOKER_H

#include "../common/typedefs.h"
#include "task_integrals_common.h"

#include "screening/screen_base.h"
#include "task_integral_kernels.h"

#include <array>
#include <vector>

namespace mpqc {
namespace integrals {
namespace detail {

template <unsigned long N>
using VecArray = std::array<ShellVec const *, N>;

template <typename E, unsigned long N, typename Op>
class op_invoke {
  private:
    std::vector<std::size_t> idx_;
    ShrPool<E> engines_;
    ShrShellVecArray<N> shell_vecs_;
    Op op_;
    std::shared_ptr<Screener> screener_;

  public:
    op_invoke() = default;
    op_invoke(std::vector<std::size_t> idx, ShrPool<E> es,
              ShrShellVecArray<N> sv, Op op, std::shared_ptr<Screener> screen)
            : idx_(std::move(idx)),
              engines_(std::move(es)),
              shell_vecs_(std::move(sv)),
              op_(std::move(op)),
              screener_(std::move(screen)) {}

    Ttype<Op> operator()(TA::Range rng) const {
        return op_(integrals(std::move(rng)));
    }

    TA::TensorD integrals(TA::Range &&rng) const {
        return integral_kernel(engines_->local(), std::move(rng),
                               shellvec_ptrs(), *(screener_.get()));
    }

    Ttype<Op> apply(TA::TensorD &&t) { return op_(std::move(t)); }

  private:
    VecArray<N> shellvec_ptrs() const {

        VecArray<N> shellvec_ptrs;
        for (auto i = 0ul; i < N; ++i) {
            shellvec_ptrs[i] = shell_vecs_[i].get();
        }

        return shellvec_ptrs;
    }
};

template <typename E, unsigned long N, typename Op, typename ScreenType>
op_invoke<E, N, Op> make_op_invoke(std::vector<std::size_t> idx, ShrPool<E> es,
                                   ShrShellVecArray<N> sh_vecs, Op op,
                                   ScreenType screen) {

    return op_invoke<E, N, Op>(
          std::move(idx), std::move(es), std::move(sh_vecs), std::move(op),
          std::make_shared<ScreenType>(std::move(screen)));
}

} // namespace detail
} // namespace integrals
} // namespace mpqc


#endif // MPQC_INTEGRALS_TASKINTEGRALSOPINVOKER_H
