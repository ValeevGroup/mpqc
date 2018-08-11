//
// Created by Eduard Valeyev on 7/18/18.
//

#ifndef MPQC4_SRC_MPQC_MATH_EXTERNAL_EIGEN_RANDOM_H_
#define MPQC4_SRC_MPQC_MATH_EXTERNAL_EIGEN_RANDOM_H_

#include <random>

#include "mpqc/math/external/eigen/eigen.h"

#include <Eigen/QR>

namespace mpqc {
namespace math {

template <typename T>
RowMatrix<T> random_unitary(int size, int seed = 42) {
  static_assert(std::is_floating_point<T>::value, "random_unitary is only implemented for real types so far");

  auto engine = std::mt19937{}; engine.seed(seed);
  auto dist = std::uniform_real_distribution<T>(0, 1);

  RowMatrix<T> rmat(size, size); // random real
  for(int r=0; r!=size; ++r)
    for(int c=0; c!=size; ++c)
      rmat(r, c) = dist(engine);

  Eigen::HouseholderQR<RowMatrix<T>> qr(rmat);
  RowMatrix<T> Q = qr.householderQ();
  // fix phase of Q: see p. 9 of http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20[paper].pdf
  RowMatrix<T> R = qr.matrixQR().template triangularView<Eigen::Upper>();
  for(int c=0; c!=size; ++c) {
    if (R(c,c) < 0) {
      Q.col(c) *= -1;
    }
  }
  return Q;
}

}
}
#endif //MPQC_RANDOM_H
