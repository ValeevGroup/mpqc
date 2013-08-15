#ifndef MPQC_MATH_TENSOR_FUNCTIONAL_HPP
#define MPQC_MATH_TENSOR_FUNCTIONAL_HPP

#include "mpqc/math/tensor/forward.hpp"
#include "mpqc/math/tensor/exception.hpp"

#include <boost/mpl/int.hpp>
#include <boost/fusion/include/vector.hpp>

namespace mpqc {
namespace detail {
namespace Tensor {

    template<class F, typename T, typename U, size_t N, class Order, class Seq>
    void apply(F f, TensorBase<T,N,Order> t, TensorBase<U,N,Order> u,
               const integral_tie<Seq> &idx) {
        f(t(idx), u(idx));
    }

    template<class F, typename T, typename U, size_t N, class Order, class Seq>
    void apply(F f, TensorBase<T,N,Order> t, TensorBase<U,N,Order> u,
               boost::mpl::int_<0>, const Seq &idx) {
        typedef integral_tie<Seq> tie;
        f(t(tie(idx)), u(tie(idx)));
    }

    template<class F, typename T, typename U, size_t N, class Order,
             class Seq, int I>
    void apply(F f, TensorBase<T,N,Order> t, TensorBase<U,N,Order> u,
               boost::mpl::int_<I>, const Seq &idx) {
#ifndef NDEBUG
        if (t.dims()[I-1] != u.dims()[I-1]) {
            throw TensorDimensionsException(I-1, t.dims()[I-1], u.dims()[I-1]);
        }
#endif
        int n = t.dims()[I-1];
        for (int i = 0; i < n; ++i) {
            apply(f, t, u, boost::mpl::int_<I-1>(),
                  boost::fusion::push_front(idx, boost::cref(i)));
        }
    }

    template<class F, typename T, typename U, size_t N, class Order>
    void apply(F f, TensorBase<T,N,Order> t, TensorBase<U,N,Order> u) {
        apply(f, t, u, boost::mpl::int_<N>(), boost::fusion::vector<>());
    }

    /// @addtogroup MathTensor
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

#endif /* MPQC_MATH_TENSOR_FUNCTIONAL_HPP */
