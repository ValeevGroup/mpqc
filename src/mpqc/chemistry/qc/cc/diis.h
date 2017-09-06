//
// Created by Chong Peng on 7/7/15.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_DIIS_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_DIIS_H_

#include <cmath>
#include <utility>

namespace mpqc {
namespace cc {

/// the {T1,T2} amplitude pair

/// @tparam T1 the type representing the set of 1-body amplitudes
/// @tparam T2 the type representing the set of 2-body amplitudes
template <typename T1, typename T2>
struct T1T2 {
  typedef typename T1::element_type element_type;
  typedef decltype(norm2(std::declval<T1>())) scalar_type;

  // constructor
  T1T2(T1 &_t1, T2 &_t2) : t1(_t1), t2(_t2) {}

  // default constructor
  T1T2() : t1(), t2() {}

  T1 t1;
  T2 t2;

  auto norm() const{
    auto t1_norm = norm2(t1);
    auto t2_norm = norm2(t2);
    return double(std::sqrt(t1_norm * t1_norm + t2_norm * t2_norm));
  }
};

template <typename T1, typename T2>
auto norm2(const T1T2<T1, T2> &a) -> decltype(norm2(std::declval<T1>())) {
  return a.norm();
}

template <typename T1, typename T2>
inline void zero(T1T2<T1, T2> &a) {
  zero(a.t1);
  zero(a.t2);
}

template <typename T1, typename T2>
inline auto dot_product(const T1T2<T1, T2> &a, const T1T2<T1, T2> &b) {
  return dot_product(a.t1, b.t1) + dot_product(a.t2, b.t2);
}

template <typename T1, typename T2, typename Scalar>
inline void axpy(T1T2<T1, T2> &y, Scalar a, const T1T2<T1, T2> &x) {
  axpy(y.t1, a, x.t1);
  axpy(y.t2, a, x.t2);
};

template <typename T1, typename T2>
inline T1T2<T1,T2> copy(T1T2<T1, T2> &a) {
  return a;
}

template <typename T1, typename T2, typename Scalar>
inline void scale(T1T2<T1, T2> &y, Scalar a) {
  scale(y.t1, a);
  scale(y.t2, a);
};

template <typename T1, typename T2>
inline std::size_t size(T1T2<T1,T2>& x) {
  return size(x.t1) + size(x.t2);
};


}  // namespace cc
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_DIIS_H_
