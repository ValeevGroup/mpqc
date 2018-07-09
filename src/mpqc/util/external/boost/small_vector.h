//
// Created by Eduard Valeyev on 7/9/18.
//

#ifndef MPQC4_SRC_MPQC_UTIL_EXTERNAL_BOOST_SMALL_VECTOR_H
#define MPQC4_SRC_MPQC_UTIL_EXTERNAL_BOOST_SMALL_VECTOR_H

// this is already defined in TiledArray/math/btas.h if defined(BTAS_HAS_BOOST_CONTAINER)
#if ! (defined(TILEDARRAY_MATH_BTAS_H__INCLUDED) && defined(BTAS_HAS_BOOST_CONTAINER))
#if __has_include(<boost/container/small_vector.hpp>)
#include <boost/container/small_vector.hpp>

namespace madness {
namespace archive {

template <class Archive, typename T, std::size_t N, typename A>
struct ArchiveLoadImpl<Archive, boost::container::small_vector<T, N, A>> {
  static inline void load(const Archive& ar,
                          boost::container::small_vector<T, N, A>& x) {
    std::size_t n{};
    ar& n;
    x.resize(n);
    for (auto& xi : x) ar& xi;
  }
};

template <class Archive, typename T, std::size_t N, typename A>
struct ArchiveStoreImpl<Archive, boost::container::small_vector<T, N, A>> {
  static inline void store(const Archive& ar,
                           const boost::container::small_vector<T, N, A>& x) {
    ar& x.size();
    for (const auto& xi : x) ar& xi;
  }
};

}  // namespace archive
}  // namespace madness

#endif  // __has_include(<boost/container/small_vector.hpp>)
#endif  // ! (defined(TILEDARRAY_MATH_BTAS_H__INCLUDED) && defined(BTAS_HAS_BOOST_CONTAINER))

#endif  // MPQC4_SRC_MPQC_UTIL_EXTERNAL_BOOST_SMALL_VECTOR_H
