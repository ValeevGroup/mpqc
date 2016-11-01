#ifndef MPQC_TAROUTINES_CHOLEKSYINVERSE_H
#define MPQC_TAROUTINES_CHOLEKSYINVERSE_H

#include <tiledarray.h>

namespace mpqc {
namespace array_ops {

TA::DistArray<TA::TensorD, TA::SparsePolicy>
cholesky_inverse(TA::DistArray<TA::TensorD, TA::SparsePolicy> const &A);

} // namespace ta_routines
} // namespace mpqc

#endif // MPQC_TAROUTINES_CHOLEKSYINVERSE_H
