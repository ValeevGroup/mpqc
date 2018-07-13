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

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
auto format_commainit(const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& matrix,
                      int precision = -1, bool align_cols = true) {
  return matrix.format(Eigen::IOFormat(precision > 0 ? precision : Eigen::FullPrecision,
                                       align_cols ? 0 : Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";"));
}

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void write_commainit(std::ostream& os,
                     const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& matrix,
                     const std::string name = "matrix",
                     int precision = -1, bool align_cols = true) {
  os << name << "(" << matrix.rows() << "," << matrix.cols() << "); "
     << name << format_commainit(matrix, precision, align_cols);
}

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
auto format_cpp(const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& matrix) {
  return matrix.format(Eigen::IOFormat(Eigen::FullPrecision, 0, ", ", ",\n", "{", "}", "{", "}"));
}

}  // namespace detail
}  // namespace mpqc

#endif // MPQC4_SRC_MPQC_MATH_EXTERNAL_EIGEN_UTIL_H_
