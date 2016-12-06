#ifndef MPQC4_SRC_MPQC_MATH_LINALG_CHOLESKY_INVERSE_H_
#define MPQC4_SRC_MPQC_MATH_LINALG_CHOLESKY_INVERSE_H_

#include <tiledarray.h>

namespace mpqc {
namespace array_ops {

TA::DistArray<TA::TensorD, TA::SparsePolicy> cholesky_inverse(
    TA::DistArray<TA::TensorD, TA::SparsePolicy> const &A);

}  // namespace ta_routines
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_MATH_LINALG_CHOLESKY_INVERSE_H_
