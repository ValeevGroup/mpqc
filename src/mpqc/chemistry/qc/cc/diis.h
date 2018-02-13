//
// Created by Chong Peng on 7/7/15.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_DIIS_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_DIIS_H_

#include <cmath>
#include <utility>

#include <mpqc/util/misc/assert.h>

namespace mpqc {
namespace cc {

/// a vector of objects that can be used for DIIS, e.g. {T1, T2, T3}

/// @tparam T the type of an amplitude tensor
template <typename T>
class TPack : public std::vector<T> {
public:
  typedef typename T::element_type element_type;
  typedef decltype(norm2(std::declval<T>())) scalar_type;

  // constructor
  TPack(T &_t1, T &_t2) : std::vector<T>{_t1,_t2} {}
  TPack(T &_t1, T &_t2, T &_t3) : std::vector<T>{_t1,_t2,_t3} {}

  TPack() = default;
  TPack(std::size_t i) : std::vector<T>(i) {}
};  // TPack<T>

template <typename T>
inline auto norm2(const TPack<T> &a) -> decltype(norm2(std::declval<T>())) {
  using scalar_type = typename TPack<T>::scalar_type;
  scalar_type norm2_squared = 0;
  for(auto& t: a) {
    auto t_norm = norm2(t);
    norm2_squared += t_norm * t_norm;
  }
  return double(std::sqrt(norm2_squared));
}

template <typename T>
inline void zero(TPack<T> &a) {
  for(auto& t: a) {
    zero(t);
  }
}

template <typename T>
inline auto dot_product(const TPack<T> &a, const TPack<T> &b) {
  MPQC_ASSERT(a.size() == b.size());
  typename TPack<T>::scalar_type result = 0;
  for(auto i=0; i!=a.size(); ++i) {
    result += dot_product(a[i], b[i]);
  }
  return result;
}

template <typename T, typename Scalar>
inline void axpy(TPack<T> &y, Scalar a, const TPack<T> &x) {
  MPQC_ASSERT(x.size() == y.size());
  for(auto i=0; i!=x.size(); ++i) {
    axpy(y[i], a, x[i]);
  }
}

template <typename T>
inline TPack<T> copy(TPack<T> &a) {
  return a;
}

template <typename T, typename Scalar>
inline void scale(TPack<T> &y, Scalar a) {
  for(auto& t: y)
    scale(t, a);
}

}  // namespace cc
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_DIIS_H_
