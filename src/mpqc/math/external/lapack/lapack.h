//
// Created by Eduard Valeyev on 5/31/18.
//

#ifndef MPQC4_SRC_MPQC_MATH_EXTERNAL_LAPACK_LAPACK_H_
#define MPQC4_SRC_MPQC_MATH_EXTERNAL_LAPACK_LAPACK_H_

#include <TiledArray/madness.h>

#if HAVE_INTEL_MKL || __has_include(<lapacke.h>)
#ifndef MADNESS_LINALG_USE_LAPACKE
#  define MADNESS_LINALG_USE_LAPACKE
#endif
#endif

#include "madness/tensor/clapack.h"

#endif //MPQC4_SRC_MPQC_MATH_EXTERNAL_LAPACK_LAPACK_H_
