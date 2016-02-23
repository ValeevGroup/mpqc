#pragma once
#ifndef MPQC_TAROUTINES_CHOLEKSYINVERSE_H
#define MPQC_TAROUTINES_CHOLEKSYINVERSE_H

#include "../include/tiledarray.h"
#include "../common/namespaces.h"
#include "../common/typedefs.h"

namespace mpqc {
namespace array_ops {

TA::DistArray<TA::TensorD, SpPolicy>
cholesky_inverse(TA::DistArray<TA::TensorD, SpPolicy> const &A);

} // namespace ta_routines
} // namespace mpqc

#endif // MPQC_TAROUTINES_CHOLEKSYINVERSE_H
