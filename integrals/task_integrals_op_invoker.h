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

template<unsigned long N>
struct schwartz_screen;

template<>
struct schwartz_screen<3> {
    MatrixD const *Qx_;
    MatrixD const *Qab_;

    schwartz_screen(MatrixD const *x, MatrixD const *ab) : Qx_(x), Qab_(ab) {}

    template <typename E, unsigned long N>
    TA::TensorD
    operator()(E &eng, TA::Range &&rng, VecArray<N> const &va) const {
        return integral_kernel(eng, std::move(rng), va, *Qx_, *Qab_);
    }
};

template <typename E, unsigned long N, typename Op, typename Screener>
class op_invoke {
  private:
    const TRange *const trange_ptr_;
    ShrPool<E> engines_;
    ShrBases<N> bases_;
    Op op_;
    Screener screen_ints_;

  public:
    op_invoke(const TRange *const tp, ShrPool<E> const &es,
              ShrBases<N> const &bs, Op op, Screener si)
            : trange_ptr_(tp),
              engines_(es),
              bases_(bs),
              op_(op),
              screen_ints_(si) {}

    Ttype<Op> operator()(int64_t tile_ord) const { return op_(integrals(tile_ord)); }

    TA::TensorD integrals(int64_t tile_ord) const {
        return screen_ints_(engines_->local(), make_range(tile_ord),
                            shellvec_ptrs(tile_ord));
    }

    Ttype<Op> apply(TA::TensorD &&t) const { return op_(std::move(t)); }

  private:
    VecArray<N> shellvec_ptrs(int64_t tile_ord) const {
        auto const &idx = trange_ptr_->tiles().idx(tile_ord);

        VecArray<N> shellvec_ptrs;
        for (auto i = 0ul; i < N; ++i) {
            shellvec_ptrs[i]
                  = &(bases_->operator[](i).cluster_shells()[idx[i]]);
        }

        return shellvec_ptrs;
    }

    TA::Range make_range(int64_t tile_ord) const {
        return trange_ptr_->make_tile_range(tile_ord);
    }
};


template <typename E, unsigned long N, typename Op>
op_invoke<E, N, Op, schwartz_screen<N>>
make_schwartz_op_invoke(const TRange *const tp, ShrPool<E> const &es,
                          ShrBases<N> const &bs, Op op, MatrixD const *Qx,
                          MatrixD const *Qab) {
    return op_invoke<E, N, Op, schwartz_screen<N>>(tp, es, bs, op,
                                                schwartz_screen<N>(Qx, Qab));
}

template <typename E, unsigned long N>
using no_screen_t = decltype(&no_screening<E, N>);

template <typename E, unsigned long N, typename Op>
op_invoke<E, N, Op, no_screen_t<E, N>>
make_unscreened_op_invoke(const TRange *const tp, ShrPool<E> const &es,
                        ShrBases<N> const &bs, Op op) {
    return op_invoke<E, N, Op, no_screen_t<E, N>>(tp, es, bs, op,
                                                  no_screening<E, N>);
}

} // namespace detail
} // namespace integrals
} // namespace mpqc


#endif // MPQC_INTEGRALS_TASKINTEGRALSOPINVOKER_H
