#pragma once
#ifndef MPQC_INTEGRALS_TASKINTEGRALSOPINVOKER_H
#define MPQC_INTEGRALS_TASKINTEGRALSOPINVOKER_H

#include "task_integrals_common.h"
#include "../common/typedefs.h"

namespace mpqc {
namespace integrals {
namespace detail {

template <unsigned long N>
using VecArray = std::array<ShellVec const *, N>;

template <typename E, unsigned long N>
TA::TensorD no_screening(E &eng, TA::Range &&rng, VecArray<N> const &va) {
    return integral_kernel(eng, std::move(rng), va);
}

template <unsigned long N>
struct schwartz_screen;

template <>
struct schwartz_screen<3> {
    std::shared_ptr<const MatrixD> Qx_;
    std::shared_ptr<const MatrixD> Qab_;

    schwartz_screen(std::shared_ptr<MatrixD> x,
                    std::shared_ptr<const MatrixD> ab)
            : Qx_(std::move(x)), Qab_(std::move(ab)) {}

    template <typename E, unsigned long N>
    TA::TensorD
    operator()(E &eng, TA::Range &&rng, VecArray<N> const &va) const {
        return integral_kernel(eng, std::move(rng), va, *Qx_, *Qab_);
    }
};

template <typename E, unsigned long N, typename Op, typename Screener>
class op_invoke {
  private:
    std::vector<std::size_t> idx_;
    ShrPool<E> engines_;
    ShrShellVecArray<N> shell_vecs_;
    Op op_;
    Screener screen_ints_;

  public:
    op_invoke(std::vector<std::size_t> idx, ShrPool<E> es,
              ShrShellVecArray<N> sv, Op op, Screener si)
            : idx_(std::move(idx)),
              engines_(std::move(es)),
              shell_vecs_(std::move(sv)),
              op_(std::move(op)),
              screen_ints_(std::move(si)) {}

    Ttype<Op> operator()(TA::Range rng) const {
        return op_(integrals(std::move(rng)));
    }

    TA::TensorD integrals(TA::Range &&rng) const {
        return screen_ints_(engines_->local(), std::move(rng), shellvec_ptrs());
    }

    Ttype<Op> apply(TA::TensorD &&t) const { return op_(std::move(t)); }

  private:
    VecArray<N> shellvec_ptrs() const {

        VecArray<N> shellvec_ptrs;
        for (auto i = 0ul; i < N; ++i) {
            shellvec_ptrs[i] = shell_vecs_[i].get();
        }

        return shellvec_ptrs;
    }
};

template <typename E, unsigned long N, typename Op, typename Screener>
op_invoke<E, N, Op, Screener>
make_screened_op_invoke(std::vector<std::size_t> idx, ShrPool<E> es,
                        ShrShellVecArray<N> sh_vecs, Op op, Screener screen) {
    return op_invoke<E, N, Op, Screener>(std::move(idx), std::move(es),
                                         std::move(sh_vecs), std::move(op),
                                         std::move(screen));
}

template <typename E, unsigned long N, typename Op>
op_invoke<E, N, Op, schwartz_screen<N>>
make_schwartz_op_invoke(std::vector<std::size_t> idx, ShrPool<E> es,
                        ShrShellVecArray<N> sh_vecs, Op op,
                        std::shared_ptr<MatrixD> Qx,
                        std::shared_ptr<MatrixD> Qab) {
    auto screen = schwartz_screen<N>(std::move(Qx), std::move(Qab));
    return op_invoke<E, N, Op, schwartz_screen<N>>(
          std::move(idx), std::move(es), std::move(sh_vecs), std::move(op),
          std::move(screen));
}


template <typename E, unsigned long N>
using no_screen_t = decltype(&no_screening<E, N>);

template <typename E, unsigned long N, typename Op>
op_invoke<E, N, Op, no_screen_t<E, N>>
make_unscreened_op_invoke(std::vector<std::size_t> idx, ShrPool<E> es,
                          ShrShellVecArray<N> sh_vecs, Op op) {
    return op_invoke<E, N, Op, no_screen_t<E, N>>(std::move(idx), std::move(es),
                                                  std::move(sh_vecs),
                                                  std::move(op),
                                                  no_screening<E, N>);
}

} // namespace detail
} // namespace integrals
} // namespace mpqc


#endif // MPQC_INTEGRALS_TASKINTEGRALSOPINVOKER_H
