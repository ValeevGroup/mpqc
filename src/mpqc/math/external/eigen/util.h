#ifndef MPQC4_SRC_MPQC_MATH_EXTERNAL_EIGEN_UTIL_H_
#define MPQC4_SRC_MPQC_MATH_EXTERNAL_EIGEN_UTIL_H_

#include "mpqc/math/external/eigen/eigen.h"

namespace mpqc {
namespace detail {

/*!
 * \brief This sorts eigenvalues and eigenvectors
 * in ascending order of the real parts of eigenvalues
 *
 * @param[in,out] eigVal the vector of complex eigenvalues
 * @param[in,out] eigVec the complex matrix consisting of complex eigenvectors
 */
void sort_eigen(VectorZ &eigVal, MatrixZ &eigVec);

/** Creates an Eigen::IOFormat object for printing Eigen matrices in format understood by comma initialization operator.
 *
 * @param[in] matrix an Eigen::Matrix object
 * @param[in] precision number of digits to print; default (-1) is to use full precision
 * @param[in] align_cols a flag that controls whether columns are aligned
 * @return an Eigen::IOFormat object; to use a format object @c fmt to print matrix @c M do @code std::cout << M.format(fmt) << std::endl; @endcode
 * @sa mpqc::detail::write_commainit
 * @note Eigen::IOFormat is described at https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
 * @note comma initialization is described at https://eigen.tuxfamily.org/dox/group__TutorialAdvancedInitialization.html
 * */
template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
auto format_commainit(const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& matrix,
                      int precision = -1, bool align_cols = true) {
  return matrix.format(Eigen::IOFormat(precision > 0 ? precision : Eigen::FullPrecision,
                                       align_cols ? 0 : Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";"));
}

/** Writes an Eigen matrix to a stream in format understood by comma initialization operator.
 *
 * @param[in,out] os a std::ostream object reference
 * @param[in] matrix an Eigen::Matrix object
 * @param[in] name a string identifying this matrix
 * @param[in] precision number of digits to print; default (-1) is to use full precision
 * @param[in] align_cols a flag that controls whether columns are aligned
 * @note comma initialization is described at https://eigen.tuxfamily.org/dox/group__TutorialAdvancedInitialization.html
 * */
template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void write_commainit(std::ostream& os,
                     const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& matrix,
                     const std::string name = "matrix",
                     int precision = -1, bool align_cols = true) {
  os << name << "(" << matrix.rows() << "," << matrix.cols() << "); "
     << name << format_commainit(matrix, precision, align_cols);
}

/** Creates an Eigen::IOFormat object for printing Eigen matrices in format useful for static initialization of C++ multidimentional arrays.
 *
 * @param[in] matrix an Eigen::Matrix object
 * @return an Eigen::IOFormat object; to use a format object @c fmt to print matrix @c M do @code std::cout << M.format(fmt) << std::endl; @endcode
 * @note Eigen::IOFormat is described at https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
 * */
template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
auto format_cpp(const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& matrix) {
  return matrix.format(Eigen::IOFormat(Eigen::FullPrecision, 0, ", ", ",\n", "{", "}", "{", "}"));
}

}  // namespace detail
}  // namespace mpqc

#endif // MPQC4_SRC_MPQC_MATH_EXTERNAL_EIGEN_UTIL_H_
