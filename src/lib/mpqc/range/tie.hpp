#ifndef MPQC_RANGE_TIE_HPP
#define MPQC_RANGE_TIE_HPP

#include <boost/array.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

#include "mpqc/range.hpp"

namespace mpqc {

    /// boost::tuple tie wrapper
    /// @ingroup MathRange
    template<class Tuple>
    struct range::tie : Tuple {
        static const int N = boost::tuples::length<Tuple>::value;
        tie(const Tuple &s) : Tuple(s) {}
        /// cast range tuple as boost::array
        operator boost::array<range,N>() const {
            boost::array<range,N> a;
            cast(a, int_<N>());
            return a;
        }
        /// cast range tuple as std::vector
        operator std::vector<range>() const {
            boost::array<range,N> a = *this;
            return std::vector<range>(a.begin(), a.end());
        }
    private:
        template<int I>
        struct int_ {};
        template<class V, int I>
        void cast(V &v, int_<I>) const {
            v[I-1] = range_cast(boost::tuples::get<I-1>(*this));
            //std::cout << "cast " << I-1 << "=" << v[I-1] << std::endl;
            cast(v, int_<I-1>());
        }
        template<class V>
        void cast(V &v, int_<0>) const {}
    };

    template<class Tuple>
    std::ostream& operator<<(std::ostream& os, const range::tie<Tuple> &tie) {
        os << (const Tuple&)tie;
        return os;
    }

}


namespace mpqc {
namespace detail {

    // wrap boost::tuple T in range::tie
    template<class T>
    range::tie<T> range_tie(const T &t) {
        return range::tie<T>(t);
    }

}
}

#endif // MPQC_RANGE_TIE_HPP
