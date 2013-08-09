#ifndef MPQC_TENSOR_FUNCTIONAL_HPP
#define MPQC_TENSOR_FUNCTIONAL_HPP

#include "mpqc/tensor/base.hpp"

namespace mpqc {
namespace detail {
namespace Tensor {

    /// @addtogroup Tensor
    /// @{

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
