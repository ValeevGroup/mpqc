#ifndef MPQC_MATH_TENSOR_FORWARD_HPP
#define MPQC_MATH_TENSOR_FORWARD_HPP

#include "mpqc/math/tensor/order.hpp"

#include <boost/fusion/include/boost_tuple.hpp>
#include <boost/mpl/accumulate.hpp>
#include <boost/mpl/print.hpp>

#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/remove_reference.hpp>

namespace mpqc {

    template <typename T, size_t N, class Order = TensorColumnMajor >
    struct TensorBase;

}


namespace mpqc {
namespace detail {
namespace Tensor {

    template<class T>
    struct is_integral {
        typedef typename boost::is_integral<
            typename boost::remove_reference<T>::type
            >::type type;
    };

    /// "returns" true if every element of tuple T is an integral type
    template<class T>
    struct is_integral_tuple {
        static const bool value =
            boost::mpl::accumulate<
            T,
            boost::mpl::bool_<true>,
            boost::mpl::and_< is_integral<boost::mpl::_2>, boost::mpl::_1 >
            >::type::value;
    };

    /// index tie wrapper
    template<class Tie>
    struct integral_tie : Tie {
        integral_tie(const Tie &t) : Tie(t) {}
    };

    /// range tie wrapper
    template<class Tie>
    struct range_tie : Tie {
        range_tie(const Tie &t) : Tie(t) {}
    };

    /// returns integral_tie if every type in Tie is an integral type,
    /// disabled otherwise
    template<class Tie>
    typename boost::enable_if<
        is_integral_tuple<Tie>,
        integral_tie<Tie>
        >::type
    tie(const Tie &t) {
        return integral_tie<Tie>(t);
    }

    /// returns range_tie if NOT every type in Tie is an integral type,
    /// disabled otherwise
    template<class Tie>
    typename boost::disable_if<
        is_integral_tuple<Tie>,
        range_tie<Tie>
        >::type
    tie(const Tie &t) {
        return range_tie<Tie>(t);
    }

} // Tensor
} // detail
} // mpqc


#endif /* MPQC_MATH_TENSOR_FORWARD_HPP */
