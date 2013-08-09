#ifndef MPQC_TENSOR_FUNCTIONAL_HPP
#define MPQC_TENSOR_FUNCTIONAL_HPP

#include "mpqc/tensor/forward.hpp"

namespace mpqc {
namespace detail {
namespace Tensor {

    template<class F, typename T, typename U, size_t N, class Order, class Seq>
    void apply(F f, TensorBase<T,N,Order> t, TensorBase<U,N,Order> u,
               const Seq &idx) {
        integral_tie<Seq> tie(idx);
        f(t(tie), u(tie));
    }

    template<class F, typename T, typename U, size_t N, class Order, class Seq>
    void apply(F f, TensorBase<T,N,Order> t, TensorBase<U,N,Order> u,
               const Seq &idx, boost::mpl::int_<0>) {
        assert(t.dims()[0] == u.dims()[0]);
        int n = t.dims()[0];
        for (int i = 0; i < n; ++i) {
            apply(f, t, u, boost::fusion::push_front(idx, boost::cref(i)));
        }
    }

    template<class F, typename T, typename U, size_t N, class Order,
             class Seq, int I>
    void apply(F f, TensorBase<T,N,Order> t, TensorBase<U,N,Order> u,
               const Seq &idx, boost::mpl::int_<I>) {
        assert(t.dims()[I] == u.dims()[I]);
        int n = t.dims()[I];
        for (int i = 0; i < n; ++i) {
            const auto &s = boost::fusion::push_front(idx, boost::cref(i));
            apply(f, t, u, s, boost::mpl::int_<I-1>());
        }
    }

    template<class F, typename T, typename U, size_t N, class Order>
    void apply(F f, TensorBase<T,N,Order> t, TensorBase<U,N,Order> u) {
        apply(f, t, u, boost::tuple<>(), boost::mpl::int_<N-1>());
    }

    /// @addtogroup Tensor
    /// @{

    struct assign {
        template <typename T, typename U>
        void operator()(T &t, const U &u) const {
            t = u;
        }
    };

    template<typename S>
    struct multiply_assign {
        multiply_assign(const S &s) : value_(s) {}
        template <typename T>
        void operator()(T &t) const {
            t *= value_;
        }
    private:
        S value_;
    };

    template<typename S>
    struct divide_assign {
        divide_assign(const S &s) : value_(s) {}
        template <typename T>
        void operator()(T &t) const {
            t /= value_;
        }
    private:
        S value_;
    };

    struct plus_assign {
        template <typename T, typename S>
        void operator()(T &t, const S &s) const {
            t += s;
        }
    };

    struct minus_assign {
        template <typename T, typename S>
        void operator()(T &t, const S &s) const {
            t -= s;
        }
    };

    /// @}

}
}
}

#endif /* MPQC_TENSOR_FUNCTIONAL_HPP */
