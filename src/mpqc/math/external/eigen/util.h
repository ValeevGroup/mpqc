#ifndef MPQC4_SRC_MPQC_MATH_EXTERNAL_EIGEN_UTIL_H_
#define MPQC4_SRC_MPQC_MATH_EXTERNAL_EIGEN_UTIL_H_

#include "mpqc/math/external/eigen/eigen.h"

namespace mpqc {
namespace detail {

/*!
 * \brief This sorts eigenvalues and eigenvectors
 * in ascending order of the real parts of eigenvalues
 *
 * \param eigVal the vector of complex eigenvalues
 * \param eigVec the complex matrix consisting of complex eigenvectors
 */
void sort_eigen(VectorZ &eigVal, MatrixZ &eigVec);

}  // namespace detail
}  // namespace mpqc

#endif // MPQC4_SRC_MPQC_MATH_EXTERNAL_EIGEN_UTIL_H_
